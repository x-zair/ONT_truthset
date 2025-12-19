# Sub-consensus haploid variant calling in Long-read sequencing technology

[![Preprint](https://img.shields.io/badge/Preprint-v1-orange.svg?logo=researchgate)](https://www.researchsquare.com/article/rs-6226988/v1)
[![BioProject](https://img.shields.io/badge/BioProject-PRJNA1245633-blue.svg)](https://www.ncbi.nlm.nih.gov/bioproject/1245633)

Contact: **Xavier Zair** · ephxz@nus.edu.sg  
Saw Swee Hock School of Public Health, National University of Singapore

---

> This repository provides scripts, reference sequences, precise workflow and templates for benchmarking haploid variant detection accuracy on Nanopore MinION/GridION platforms (R9.4.1 and R10.4.1 flow cells).  
> Raw reads are available at NCBI BioProject [PRJNA1245633](https://www.ncbi.nlm.nih.gov/bioproject/1245633).

The dataset includes a simple, low-complexity chimeric plasmid truth set designed with short 261 bp inserts, and longer bacterial genome sequences from two closely related *E. coli* strains (∼95% similarity) that contain insertion and deletion variants. This resource allows researchers to benchmark rare variant calling performance specifically on Nanopore sequencing data, reflecting both controlled and complex genomic scenarios.This repository also provides detailed used for the Phred calibration tool [QUAD](https://github.com/andreas-wilm/quad), and accompanying callibration files for Phred score adjusments on `.fastq` files.

<details>
 
<summary>Truthset Composition and details</summary>

## Plasmid Dataset

This dataset contains chimeric plasmid mixtures designed within a controlled truth-set framework.

- The Delta spike plasmid (pUNO-01, 7179 bp) is used as the base donor.  
- A 261 bp insert region (positions 2554–2785) was selected for modification.  
- Two chimeric plasmids were created by replacing this insert with sequences from Wild-type and Omicron variants.  

Plasmid mixtures were prepared at varying ratios to simulate ground-truth complexity. Nine mixture sets, along with three clonal plasmids, are available at the [repository](https://www.ncbi.nlm.nih.gov/bioproject/1245633).

<img src="/img/vector.png" alt="plasmid design" style="width:40%; height:auto;"/>


**Dataset highlights:**

| Feature                     | Details                                         |
|-----------------------------|-------------------------------------------------|
| Insert length               | 261 base pairs                                  |
| Variant positions           | 85                                              |
| Unique mutants              | 100                                             |
| Single nucleotide variants  | 60 positions with single nucleotide differences across Omicron, Delta, and Wild-type plasmids |
| Multiple nucleotide variants| 15 positions with distinct nucleotides in all three plasmids; variant callers must detect two alternate alleles per site, doubling mutants for these positions (15 × 2 = 30 mutants)  |
 
![dot_aln](img/dot_aln.png)

List of variants provided at within the [plasmid variant sheet](/plasmids/list_vector.xlsx). Reference fasta's provided in the same folder. 

| Barcode    | Omicron | Delta | Wild-type |
|------------|---------|-------|-----------|
| Barcode 1  | 1       | 0     | 0         |
| Barcode 2  | 0       | 1     | 0         |
| Barcode 3  | 0       | 0     | 1         |
| Barcode 4  | 0.8     | 0.1   | 0.1       |
| Barcode 5  | 0.8     | 0.15  | 0.05      |
| Barcode 6  | 0.8     | 0.19  | 0.01      |
| Barcode 7  | 0.1     | 0.8   | 0.1       |
| Barcode 8  | 0.15    | 0.8   | 0.05      |
| Barcode 9  | 0.19    | 0.8   | 0.01      |
| Barcode 10 | 0.1     | 0.1   | 0.8       |
| Barcode 11 | 0.05    | 0.15  | 0.8       |
| Barcode 12 | 0.01    | 0.19  | 0.8       |

### Bacterial Genomes

Two closely related *E. coli* strains (~95% similarity) were selected for sequencing: **DH5-α** and **STBL3 (HB101)**. 

<img src="img/genome_pair.png" alt="dot" style="width:40%; height:auto;"/>

Each strain was sequenced on separate flow cells. Additionally, a third flow cell contained a mixture of the two strains at varying ratios:

| Barcode        | DH5-α  | STBL3  |
|----------------|--------|--------|
|Barcode 1       | 0.990  | 0.010  |
|Barcode 3       | 0.950  | 0.050  |
|Barcode 4       | 0.900  | 0.100  |

Reference genomes and variant list enclosed in the [genomes](/genomes) folder. Variant list is seperated into subsitutional(list_s), deletional(list_d) and insertional(list_i) variants. Clonal genomes of both strains, and the mixes are available at the [repository](https://www.ncbi.nlm.nih.gov/bioproject/1245633).

</details>

<details>

<summary>Workflow Overview</summary>

## ONT Variant Calling Benchmarking with Phred Score Calibration
                                                                 
```
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                                        STUDY DESIGN                                                                 │
└──────────────────────────────────────────────────────┬──────────────────────────────────────────────────────────────────────────────┘
                                                       │
                                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                                 TRUTH SET CONSTRUCTION                                                               │
├──────────────────────────────────────────────────────────────────┬───────────────────────────────────────────────────────────────────┤
│                  PLASMID APPROACH                                │                    GENOMIC APPROACH                               │
│                                                                  │                                                                   │
│  3 SARS-CoV-2 Spike Variants (Delta, Omicron, Wild-type)         │     2 E. coli K-12 Strains (HB101 & DH5α)                         │
│  • 100 SNVs across 7,179 bp                                      │     • 1,496 SNVs, 13,770 Insertions, 47,990 Deletions             │
│                                                                  │     • 95.2% pairwise genome similarity                            │
│                                                                  │                                                                   │
│  Mix at AF: 0.01, 0.05, 0.1, 0.15, 0.19                          │     Mix at AF: 0.1                                                │
└──────────────────────────────────────────────────────────────────┴───────────────────────────────────────────────────────────────────┘
                                                       │
                                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                          LIBRARY PREPARATION & SEQUENCING                                                            │
├──────────────────────────────────────────────────────────────────┬───────────────────────────────────────────────────────────────────┤
│                  R9.4.1 Chemistry                                │                  R10.4.1 Chemistry                                │
│         SQT-LSK109 Kit + EXP-NBD104 Barcodes                     │              SQK-NBD114.24 Native Barcoding Kit                   │
└──────────────────────────────────────────────────────────────────┴───────────────────────────────────────────────────────────────────┘
                                                       │
                                                       ▼
                                          GridION Mk1 Platform (72-hour run)
                                                       │
                                                       ▼
                                             Raw Sequencing Data (FAST5/POD5)
                                                       │
                                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                               BIOINFORMATIC PROCESSING PIPELINE                                                      │
├──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│  1. Basecalling → Dorado v0.5.3 (Super Accurate Model v4.2.0)          4. Alignment → Minimap2 v2.17 (-ax map-ont)                   │
│  2. QC & Filtering → NanoFilt v2.8.0 (Q ≥ 10)                          5. Depth Normalization → Rasusa v0.6.1(12 replicates/plasmids)│
│  3. Adapter Trimming → Porechop v0.2.3 (auto mode)                                            Jvarkit v4b65b20b2 (genomes)           │
└───────────────────────────────────────────────┬──────────────────────────────────────────────┬───────────────────────────────────────┘
                                                │                                              │
                                                ▼                                              ▼
                    ┌───────────────────────────────────────┐              ┌───────────────────────────────────────────────────┐
                    │       STANDARD WORKFLOW               │              │         PHRED CALIBRATION WORKFLOW                │
                    ├───────────────────────────────────────┤              ├───────────────────────────────────────────────────┤
                    │                                       │              │  1. Extract 1M bases (12 iterations, Rasusa)      │
                    │  LoFreq Variant Calling v2.1.3.1      │              │  2. QvsQ Plot Generation (BBmap v38.18 qahist)    │
                    │  • -B mode (disable BAQ)              │              │                                                   │
                    │  • Standard Phred scores              │              │  3. Calibration Approaches Tested:                │
                    │                                       │              │     a) Fitted Curves (R² > 0.95):                 │
                    │                                       │              │        • Quadratic: y = ax² + bx + c              │
                    │                                       │              │        • Polynomial: up to n=4                    │
                    │                                       │              │        • Gaussian: y = Ae^(-(x-μ)²/2σ²)           │
                    │                                       │              │     b) Scalar Reduction: -1 to -10                │
                    │                                       │              │     c) Control: No adjustment                     │
                    │                                       │              │                                                   │
                    │                                       │              │  4. Apply Calibration (QUAD v0.2)                 │
                    │                                       │              │  5. LoFreq Variant Calling v2.1.3.1               │
                    └───────────────────┬───────────────────┘              └──────────────────────┬────────────────────────────┘
                                        │                                                         │
                                        └─────────────────────┬───────────────────────────────────┘
                                                              ▼
                                    ┌─────────────────────────────────────────────────────────┐
                                    │         BENCHMARK AGAINST TRUTH SET                     │
                                    ├─────────────────────────────────────────────────────────┤
                                    │  • Calculate Recall (Sensitivity)                       │
                                    │  • Calculate False Discovery Rate (FDR)                 │
                                    │  • Compare: Chemistries, Depths, Calibration Methods    │
                                    └──────────────────────────┬──────────────────────────────┘
                                                               ▼
                                    ┌─────────────────────────────────────────────────────────┐
                                    │              STATISTICAL ANALYSIS                       │
                                    ├─────────────────────────────────────────────────────────┤
                                    │  • Unpaired t-tests (chemistry comparisons)             │
                                    │  • Linear regression (Phred reduction effect, R²=0.73)  │
                                    │  • Cohen's d effect sizes (d > 2.9, very large)         │
                                    │  • ANOVA with Tukey's HSD (depth comparisons)           │
                                    │                                                         │
                                    │  Tools: GraphPad Prism v8.0.1, R v4.2.0 (effsize)       │
                                    └──────────────────────────┬──────────────────────────────┘
                                                               ▼
                                              ┌────────────────────────────────┐
                                              │         FINAL RESULTS          │
                                              ├────────────────────────────────┤
                                              │  • Recall Metrics              │
                                              │  • False Discovery Rates       │
                                              │  • Effect Sizes (Cohen's d)    │
                                              │  • Optimal Calibration Method  │
                                              └────────────────────────────────┘
```
</details>

<details>

<summary>Codes for VCF extraction and analysis</summary> 

## Bioinfromatic Processing

> 1) Concatenate raw fastq files by barcode.
> 2) Trim adapters using [Porechop](https://github.com/rrwick/Porechop).
> 3)  Filter reads by quality (min PHRED score 10) with [NanoFilt](https://github.com/wdecoster/nanofilt);/
>     for plasmids, addition filter by length (max length 7180).
> 5)  Align reads to the reference genome using [minimap2](https://github.com/lh3/minimap2) in ONT mode for plasmids;\
>     Denovo assembly of genomes with [Flye](https://github.com/mikolmogorov/Flye) asssembler.
> 6)  Convert, sort, and index SAM files to BAM format using [samtools](https://www.htslib.org/doc/samtools.html).\
      [Qualimap](https://docs.seqera.io/multiqc/modules/qualimap) for QC information
> 8)  Perform variant calling using [LoFreq](https://csb5.github.io/lofreq/) with multi-threading.

```bash

cat *.fastq > [barcode].fastq

porechop -i [barcode].fastq -o trim_[barcode].fastq

NanoFilt -q 10 --maxlength [length] trim_[barcode].fastq > filt_trim_[barcode].fastq

minimap2 -ax map-ont [reference.fasta] filt_trim_[barcode].fastq > [barcode].sam
# OR
flye --nano-raw [in.fastq] --out-dir [path/to/dir] --threads [int]

samtools sort --threads 50 [barcode].sam > [barcode].bam
samtools index [barcode].bam

qualimap bamqc -nr 250 -bam [barcode].bam

lofreq call-parallel \
  --pp-threads [num_threads] -B -f [reference.fasta] -o [barcode]_lf.vcf [barcode].bam

```
> [!Note]
> Replace all placeholders enclosed in square brackets (e.g., `[barcode]`, `[reference.fasta]`, `[num_threads]`) with the actual values for your samples, reference genome files, and desired thread counts before running the commands.
## Additional variant caller benchmarks
[Longshot](https://github.com/pjedge/longshot):
```bash
longshot -n –-bam [file] –ref [file] –out [file].vcf
```
[Clair3](https://github.com/HKU-BAL/Clair3)
```bash
run_clair3.sh \
--bam_fn=[path/to/bam] \
--ref_fn=[path/to/fasta] \
--threads=[int] \
--platform="ont" \
--model_path="[path/to/model]" \
--haploid_sensitive \
--output="[out]" \
--include_all_ctgs \
--no_phasing_for_fa

```
> [!Note]
> The reference has to be the same seq that was used to do the alignment. Also important for the alignment toolkit to fit the platform- ergo ONT uses minimap2, realign if necessary. The –include_all_ctgs is vital.
### Custom scripts for VCF extraction
R scripts in respective [plasmids](/plasmids) and [genomes](/genomes) folder will extract variant information from variant files and summarise true positives, false positives, and false negatives based on provided variant lists. Scripts within genome folder will further segregate variants into subsiutial, deletional and insertional variants. 
Sample VCFs as well as expected outputs provided within [samples](/samples) folder.
Respective reference fasta's can also be found within respective folders. 
> [!Important]
> Scripts are designed to parse VCF files generated by [LoFreq](https://csb5.github.io/lofreq/); compatibility with VCFs from other variant callers may require parser modifications.
## Toolkits for detailed analysis
- [Rasusa](https://github.com/mbhall88/rasusa) for plasmid subsampling
- [Jvarkit/Biostars154220](https://lindenb.github.io/jvarkit/Biostar154220.html) for genome subsampling
- [BBmap](https://github.com/BioInfoTools/BBMap) Phred score analysis:
   - raw Phred score extaction with reformat.sh -qchist
   - QvQ plot on clonal datasets with reformat.sh -qahist
-  
   
```bash
# Subsampling 
rasusa --input [input.fq] --coverage [Target depth] --genome-size [genome size]
java -jar [/path/to/jvarkit.jar] biostar154220 -d [Target Depth]  -o [output.sam] [input.sam]

# BBmap install:
conda install -c bioconda bbmap
# Phred score extraction:
reformat.sh in=[fastq] qchist=[output.csv] maxcalledquality=[int; 90 for ONT]
# QvQ:
reformat.sh in=[bam] qahist=[output.csv] maxcalledquality=[int; 90 for ONT] ref=[reference.fasta]
```
> [!Tip]
> If `qahist` fails, filtering adding MD tags to bam alignment might be required:
```bash
samtools view -h -F 3844 yourinput.bam | samtools calmd - yourref.fasta | samtools view - -o youroutput.bam
```

</details>

