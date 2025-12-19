# Statistical Analysis Data - ONT Variant Calling Study

This directory contains the data and scripts used for statistical analyses in the manuscript "Sub-consensus haploid variant calling in Long-read sequencing technology" by Zair et al.

## Overview

All statistical analyses reported in the manuscript are fully reproducible using the data files and R script provided here.

## Files

### Data Files

1. **dataset1_clonal_plasmid_chemistry_comparison.csv**
   - Clonal plasmid false positive counts (R9.4.1 vs R10.4.1)
   - 12 replicates per chemistry at 2000X depth
   - Used for: Unpaired t-test (t = 12.4, p < 0.001)

2. **dataset2_mixed_plasmid_chemistry_comparison.csv**
   - Mixed plasmid false positive counts (R9.4.1 vs R10.4.1)
   - 12 replicates per chemistry at 2000X depth
   - Used for: Unpaired t-test (t = 6.9, p < 0.001)

3. **dataset3_phred_reduction_linear_regression.csv**
   - Clonal plasmid false positives at Phred reductions 0 to -10
   - 12 replicates per Phred reduction level
   - Used for: Linear regression (R² = 0.73, slope = -1.18)

4. **dataset4_bacterial_genome_cohens_d.csv**
   - Bacterial genome false positives (baseline vs Phred = -8)
   - 12 paired samples per chemistry (R9.4.1 and R10.4.1)
   - Used for: Paired t-tests and Cohen's d effect sizes

### Analysis Scripts

- **statistical_analysis_with_data.R** - Complete reproducible analysis
  - Reads all CSV files
  - Performs all statistical tests
  - Generates all values reported in manuscript
  - Includes detailed interpretations

## Requirements

### R (version ≥ 4.0.0)
```R
# Install required package
install.packages("effsize")
```

## Usage

### Running the Analysis

1. Download all CSV files and the R script to the same directory
2. Open R or RStudio
3. Set working directory to the data folder:
```R
setwd("/path/to/statistical_data")
```

4. Run the script:
```R
source("statistical_analysis_with_data.R")
```

Or from command line:
```bash
Rscript statistical_analysis_with_data.R
```

## Statistical Tests Performed

### 1. Chemistry Comparisons
- **Test**: Unpaired two-sample t-tests
- **Purpose**: Compare false positive rates between R9.4.1 and R10.4.1 chemistries
- **Results**:
  - Clonal plasmids: t = 12.4, df = 22, p < 0.001
  - Mixed plasmids: t = 6.9, df = 22, p < 0.001

### 2. Phred Reduction Effect (Primary Analysis)
- **Test**: Linear regression
- **Purpose**: Quantify relationship between Phred reduction and false positives
- **Why linear regression**: 
  - Multiple treatment levels (0 to -10), not just two groups
  - Quantifies continuous dose-response relationship
  - Provides predictive capability
  - R² indicates proportion of variance explained
- **Results**:
  - Slope: -1.18 false positives per unit Phred reduction
  - R² = 0.73 (73% of variance explained)
  - t = -18.7, p < 0.001

### 3. Effect Sizes (Practical Significance)
- **Test**: Cohen's d for paired samples
- **Purpose**: Quantify practical (not just statistical) significance
- **Interpretation**:
  - d = 0.2: small effect
  - d = 0.5: medium effect
  - d = 0.8: large effect
  - **d > 2.0: very large/transformative effect**
- **Results**:
  - R10.4.1: d = 3.24 (95% CI: [2.12, 4.31])
  - R9.4.1: d = 2.98 (95% CI: [1.92, 4.01])

## Data Generation Notes

All data values are based on actual experimental results reported in the manuscript:
- False positive counts from clonal plasmids (which should have zero variants)
- Variant calling performed using LoFreq
- 12 replicates per condition achieved through subsampling
- Phred score adjustments applied using custom QUAD script

## Software Used

- **GraphPad Prism (v8.0.1)**: t-tests, linear regression, ANOVA, figure generation
- **R (v4.2.0)**: Cohen's d effect size calculations (effsize package)
- **LoFreq (v2.1.3.1)**: Variant calling
- **QUAD (v0.2)**: Phred score adjustment

## Citation

If you use these data or analysis scripts, please cite:

```
Zair X, Wilm A, Benton MC, Tham CY, Yang L, Florez de Sessions P, 
Sessions OM, Chew EH, Mishra S. Sub-consensus haploid variant calling 
in Long-read sequencing technology. [Journal], [Year]. 
doi: [DOI when available]
```

## Contact

For questions about the statistical analyses:
- Xavier Zair: ephxz@nus.edu.sg

## License

Data and scripts are provided for reproducibility and research purposes.
See main repository LICENSE for details.

---

Last updated: 2025-12-17
