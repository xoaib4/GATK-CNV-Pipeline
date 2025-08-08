#!/usr/bin/env bash

# Usage:
# ./gatk_cnv_pipeline.sh -r path/to/reference.fasta -b path/to/targets.bed -l bam_list.txt

# -------------
# Input Parsing
# -------------

while getopts "r:b:l:" opt; do
  case $opt in
    r) REFERENCE=$OPTARG ;;
    b) TARGETS=$OPTARG ;;
    l) BAM_LIST_FILE=$OPTARG ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Check required files
if [[ ! -f $REFERENCE || ! -f $TARGETS || ! -f $BAM_LIST_FILE ]]; then
  echo "Missing required input(s). Please check reference, targets, or BAM list." >&2
  exit 1
fi

# ------------------------
# Set up derived variables
# ------------------------

WORKDIR="gatk_cnv_work"
mkdir -p "$WORKDIR/cvg" "$WORKDIR/ploidy"

PREPROCESSED_INTERVALS="$WORKDIR/targets.preprocessed.interval_list"
ANNOTATED_INTERVALS="$WORKDIR/targets.annotated.tsv"
FILTERED_INTERVALS="$WORKDIR/targets.filtered.interval_list"
CONTIG_PRIORS="$WORKDIR/contig_ploidy_priors.tsv"

# Dummy priors file
echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2" > $CONTIG_PRIORS
echo -e "chr1\t0.01\t0.1\t0.89" >> $CONTIG_PRIORS

# -----------------------------------
# Step 1: Preprocess target intervals
# -----------------------------------
gatk PreprocessIntervals \
  -R "$REFERENCE" \
  -L "$TARGETS" \
  --bin-length 0 \
  -imr OVERLAPPING_ONLY \
  -O "$PREPROCESSED_INTERVALS"

# --------------------------
# Step 2: CollectReadCounts
# --------------------------

BAM_COUNT=0
READCOUNT_INPUTS=()
while read -r BAM; do
  BAM_NAME=$(basename "$BAM" .bam)
  OUTPUT_TSV="$WORKDIR/cvg/${BAM_NAME}.tsv"

  gatk CollectReadCounts \
    -L "$PREPROCESSED_INTERVALS" \
    -R "$REFERENCE" \
    -I "$BAM" \
    --format TSV \
    -imr OVERLAPPING_ONLY \
    -O "$OUTPUT_TSV"

  READCOUNT_INPUTS+=("-I" "$OUTPUT_TSV")
  ((BAM_COUNT++))
done < "$BAM_LIST_FILE"

# --------------------------
# Step 3: Annotate intervals
# --------------------------
gatk AnnotateIntervals \
  -L "$PREPROCESSED_INTERVALS" \
  -R "$REFERENCE" \
  -imr OVERLAPPING_ONLY \
  -O "$ANNOTATED_INTERVALS"

# --------------------------
# Step 4: Filter intervals
# --------------------------
gatk FilterIntervals \
  -L "$PREPROCESSED_INTERVALS" \
  --annotated-intervals "$ANNOTATED_INTERVALS" \
  "${READCOUNT_INPUTS[@]}" \
  -imr OVERLAPPING_ONLY \
  -O "$FILTERED_INTERVALS"

# ---------------------------------------
# Step 5: Determine Germline Contig Ploidy
# ---------------------------------------
gatk DetermineGermlineContigPloidy \
  -L "$FILTERED_INTERVALS" \
  "${READCOUNT_INPUTS[@]}" \
  --interval-merging-rule OVERLAPPING_ONLY \
  --contig-ploidy-priors "$CONTIG_PRIORS" \
  --output "$WORKDIR/ploidy" \
  --output-prefix ploidy

# --------------------------
# Step 6: GermlineCNVCaller
# --------------------------

gatk GermlineCNVCaller \
  --run-mode COHORT \
  -L "$FILTERED_INTERVALS" \
  "${READCOUNT_INPUTS[@]}" \
  --contig-ploidy-calls "$WORKDIR/ploidy" \
  --annotated-intervals "$ANNOTATED_INTERVALS" \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output "$WORKDIR/cohort_cnv" \
  --output-prefix cohort_all

# -----------------------------
# Step 7: Post-process example
# -----------------------------

# You must loop this per-sample to get VCFs
# Replace "NA19017" with each sample name
gatk PostprocessGermlineCNVCalls \
  --model-shard-path "$WORKDIR/cohort_cnv/cohort_all-model" \
  --calls-shard-path "$WORKDIR/cohort_cnv/cohort_all-calls" \
  --allosomal-contig chrX --allosomal-contig chrY \
  --contig-ploidy-calls "$WORKDIR/ploidy" \
  --sample-index 0 \
  --output-genotyped-intervals "$WORKDIR/NA19017.genotyped.intervals.vcf.gz" \
  --output-genotyped-segments "$WORKDIR/NA19017.genotyped.segments.vcf.gz" \
  --sequence-dictionary "${REFERENCE%.fasta}.dict"
