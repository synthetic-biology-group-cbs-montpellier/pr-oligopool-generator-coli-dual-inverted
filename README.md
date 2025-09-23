# pr-oligopool-generator-coli-dual-inverted

<img width="940" height="372" alt="image" src="https://github.com/user-attachments/assets/4cb986d5-21aa-4709-a1e8-7181cda1b147" />


Python tool for generating dual inverted promoter sequences from flow-seq data, with automated sequence design for synthesis.


## Overview

This script processes flow-seq data from Kosuri et al. 2012 to select two distinct sets of genetic constructs for generating dual inverted promoter-RBS architectures separated by a bidirectional terminator. It uses PCA clustering and 2D grid selection to minimize recombination risk while maximizing expression diversity.

## Key Features

- **Sequence Selection**: Selects maximally different sequence sets using PCA and k-mer analysis
- **Architecture Design**: Generates synthesis-ready sequences with primers and BsaI sites
- **Junction Analysis**: Comprehensive restriction site detection at sequence boundaries
- **Quality Control**: Length verification and sequence validation
- **Visualization**: Multi-page PDF reports with detailed analysis plots

## Architecture

**Final Synthesis Construct:**
```
PRIMER_5_PRIME - BSAI_SITE_FWD - RevComp(Set1) - Terminator - Set2 - BSAI_SITE_REV - PRIMER_3_PRIME
```

## Input Files

- `sd01.xlsx`: Promoter sequences and data
- `sd02.xlsx`: RBS sequences and data  
- `sd03.xlsx`: Flow-seq construct data with fluorescence measurements

## Output Files

### CSV Files
- `05_concatenated_dual_inverted_promoters.csv`: Complete sequence data with all components
- `05_synthesis_ready_sequences.csv`: All sequences ready for synthesis
- `05_synthesis_ready_sequences_CLEAN.csv`: Only sequences passing all quality checks
- `05_selected_variants_SET1.csv` / `05_selected_variants_SET2.csv`: Individual set data

### Reports
- `dual_inverted_promoter_analysis.pdf`: Comprehensive multi-page analysis report
- `overview_plot.png`: Key results summary visualization

## Key Parameters

```python
SEQUENCES_PER_SET = 44        # Sequences in each set
MIN_SEQUENCE_LENGTH = 20      # Minimum sequence length (bp)
MAX_SEQUENCE_LENGTH = 70      # Maximum sequence length (bp)
N_BINS = 8                    # Grid dimensions for 2D selection
```

## Quality Control Features

- **Junction Analysis**: Detects unwanted restriction sites at component boundaries
- **Length Verification**: Ensures component lengths sum correctly
- **Sequence Validation**: Removes sequences with problematic sites
- **Architecture Compliance**: Verifies complete synthesis construct integrity

## Usage

1. Place input Excel files (`sd01.xlsx`, `sd02.xlsx`, `sd03.xlsx`) in the script directory
2. Adjust parameters as needed at the top of the script
3. Run the script:
   ```bash
   python Dual_inverted2.py
   ```

## Dependencies

- pandas
- numpy
- matplotlib
- sklearn
- re, os, shutil, datetime

## Output Interpretation

- **PASS**: Sequences with no restriction site issues - ready for synthesis
- **WARN**: Sequences with minor issues - evaluate case by case
- **FAIL**: Sequences with junction restriction sites - avoid synthesis

## Notes

- Script creates timestamped output directories for each run
- All sequences are sanitized to remove whitespace
- Length verification ensures data integrity
- Junction analysis prevents unwanted cutting during cloning

---



