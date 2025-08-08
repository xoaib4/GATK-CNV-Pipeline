# Germline CNV Discovery Pipeline (GATK 4)

This pipeline performs **multi-sample germline CNV discovery** using [GATK GermlineCNVCaller](https://gatk.broadinstitute.org/) on a cohort of samples with aligned BAM files.

## Requirements

- GATK 4.x installed and available in `PATH`
- Bash
- Samtools (for `faidx` if used)
- Input files:
  - Reference genome FASTA (indexed with `.fai` and `.dict`)
  - Interval list (or BED) file of target regions
  - List of BAM files (one path per line)

## Input

- `-r` → Path to reference genome (e.g. `Homo_sapiens_assembly38.fasta`)
- `-b` → Target intervals in BED or `.interval_list` format
- `-l` → List of BAM file paths (text file, one per line)

## Usage

```
chmod +x gatk_cnv_pipeline.sh

./gatk_cnv_pipeline.sh \
  -r ref/Homo_sapiens_assembly38.fasta \
  -b targets.bed \
  -l bam_list.txt
```
## Output Structure
After running, the script will create:

```
gatk_cnv_work/
├── targets.preprocessed.interval_list
├── targets.annotated.tsv
├── targets.filtered.interval_list
├── contig_ploidy_priors.tsv
├── cvg/
│   └── *.tsv                # Read count files per sample
├── ploidy/
│   ├── ploidy-calls/       # Contig ploidy results
├── cohort_cnv/
│   ├── cohort_all-model/
│   ├── cohort_all-calls/
├── NA19017.genotyped.*.vcf.gz
```
Replace NA19017 with your sample name(s) and loop through to generate individual VCFs.

## Example BAM List Format
```
/path/to/sample1.bam
/path/to/sample2.bam
/path/to/sample3.bam
...
```
## Notes
The CNV calling is done in cohort mode, and assumes all samples share similar target intervals.

The script creates a dummy contig_ploidy_priors.tsv file for autosomal chromosomes. You can adjust this manually.

You can further split intervals for scatter-gather parallelism using IntervalListTools if needed.
