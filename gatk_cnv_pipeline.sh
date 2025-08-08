#!/usr/bin/env bash

# Usage:
# ./gatk_cnv_pipeline.sh -r reference.fasta -b targets.bed -l bam_list.txt -p ploidy_priors.tsv

set -euo pipefail

# Parse arguments
while getopts ":r:b:l:p:" opt; do
  case ${opt} in
    r ) REF=${OPTARG};;
    b ) BED=${OPTARG};;
    l ) BAM_LIST=${OPTARG};;
    p ) PLOIDY_PRIORS=${OPTARG};;
    \? ) echo "Usage: cmd [-r reference] [-b bed_file] [-l bam_list] [-p ploidy_priors]"; exit 1;;
  esac
done

if [[ -z "${REF:-}" || -z "${BED:-}" || -z "${BAM_LIST:-}" || -z "${PLOIDY_PRIORS:-}" ]]; then
  echo "Missing required arguments. Use -r, -b, -l, and -p."
  exit 1
fi

WORKDIR="gatk_cnv_work"
COVERAGE_DIR="$WORKDIR/cvg"
PLOIDY_DIR="$WORKDIR/ploidy"
CNV_DIR="$WORKDIR/cohort_cnv"

mkdir -p "$WORKDIR" "$COVERAGE_DIR" "$PLOIDY_DIR" "$CNV_DIR"

# Step 1: Preprocess Intervals
gatk PreprocessIntervals \
  -R "$REF" \
  -L "$BED" \
  --bin-length 0 \
  -imr OVERLAPPING_ONLY \
  -O "$WORKDIR/targets.preprocessed.interval_list"

# Step 2: Collect Read Counts
while IFS= read -r BAM; do
  SAMPLE=$(basename "$BAM" .bam)
  gatk CollectReadCounts \
    -L "$WORKDIR/targets.preprocessed.interval_list" \
    -R "$REF" \
    -imr OVERLAPPING_ONLY \
    -I "$BAM" \
    --format TSV \
    -O "$COVERAGE_DIR/$SAMPLE.tsv"
done < "$BAM_LIST"

# Step 3: Annotate Intervals
gatk AnnotateIntervals \
  -L "$WORKDIR/targets.preprocessed.interval_list" \
  -R "$REF" \
  -imr OVERLAPPING_ONLY \
  -O "$WORKDIR/targets.annotated.tsv"

# Step 4: Filter Intervals
INPUTS=()
for TSV in "$COVERAGE_DIR"/*.tsv; do
  INPUTS+=("-I" "$TSV")
done

gatk FilterIntervals \
  -L "$WORKDIR/targets.preprocessed.interval_list" \
  --annotated-intervals "$WORKDIR/targets.annotated.tsv" \
  "${INPUTS[@]}" \
  -imr OVERLAPPING_ONLY \
  -O "$WORKDIR/targets.filtered.interval_list"

# Step 5: Determine Contig Ploidy
gatk DetermineGermlineContigPloidy \
  -L "$WORKDIR/targets.filtered.interval_list" \
  --interval-merging-rule OVERLAPPING_ONLY \
  "${INPUTS[@]}" \
  --contig-ploidy-priors "$PLOIDY_PRIORS" \
  --output "$PLOIDY_DIR" \
  --output-prefix ploidy

# Step 6: Germline CNV Calling
gatk GermlineCNVCaller \
  --run-mode COHORT \
  -L "$WORKDIR/targets.filtered.interval_list" \
  "${INPUTS[@]}" \
  --contig-ploidy-calls "$PLOIDY_DIR/ploidy-calls" \
  --annotated-intervals "$WORKDIR/targets.annotated.tsv" \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output "$CNV_DIR" \
  --output-prefix cohort_all
