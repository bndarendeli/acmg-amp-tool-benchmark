# Data Requirements for Figure 2

## Required Data Files

### 1. merged_results.tsv (PRIMARY DATA SOURCE)
**Location**: `../../analysis_last/data/merged_results.tsv`  
**Size**: ~6.6 MB  
**Description**: Main dataset containing all tool classifications across all variants

#### Content:
- All variant classifications from 14 tools
- Multiple datasets: ClinGen, FOXL2, HGMD+ClinVar combinations
- Columns include:
  - Variant identifiers (Chr, Pos, Ref, Alt)
  - Dataset name
  - Ground truth classification
  - Tool classifications (InterVar_2018, InterVar_2025, BIAS, CharGer, etc.)
  - ACMG criteria for each tool

#### How it's used:
- **Panel A**: Filtered to calculate call rates (via `corrected_call_rates.csv`)
- **Panel B**: Analyzes unresolved variants (Unknown + NULL)
- **Panel C**: Compares InterVar 2018 vs 2025 classifications

#### Source:
This file is a **pre-processed master dataset** that combines:
- Tool output files from all 14 classification tools
- Ground truth annotations
- Multiple evaluation cohorts

#### How to Generate merged_results.tsv:
✅ **Script is NOW INCLUDED in this directory!**

```bash
python3 parse_all_results.py
```

This script:
1. Loads ground truth data from all datasets
2. Parses tool output files from 14 classification tools
3. Merges all tool results with ground truth
4. Saves to `merged_results.tsv`

**Dependencies**:
- `parsers/` directory (included) - Tool-specific parsers for 14 tools
- `utils_criteria_parser.py` (included) - ACMG criteria parsing utilities
- Ground truth files (must be in `../../analysis_2/ground_truth/`)
- Tool result files (must be in `../../Results/`)

---

### 2. corrected_call_rates.csv (DERIVED DATA)
**Location**: `./corrected_call_rates.csv` (included in this directory)  
**Size**: ~2.6 KB  
**Description**: Pre-calculated call rates for all tools across datasets

#### Content:
```csv
Tool,Dataset,Call_Rate
Genebe,clingen_28012026,99.8
BIAS,clingen_28012026,99.7
...
```

#### How it's generated:
```bash
python3 recalculate_call_rates.py
```

This script:
1. Reads `merged_results.tsv`
2. Calculates call rates for each tool on each dataset
3. Saves to `corrected_call_rates.csv`

#### How it's used:
- **Panel A only**: Provides call rate data for the horizontal stacked bar chart

---

## Data Flow Diagram

```
merged_results.tsv (6.6 MB)
    ↓
    ├─→ recalculate_call_rates.py
    │       ↓
    │   corrected_call_rates.csv (2.6 KB)
    │       ↓
    │   [Panel A: Call Rate Spectrum]
    │
    ├─→ create_figure1_comprehensive.py
    │       ↓
    │   [Panel B: Unresolved Variants Analysis]
    │       ↓
    │   [Panel C: InterVar Comparison]
    │
    ↓
Figure2_Comprehensive.png/pdf
```

---

## Checking Data Availability

### Quick Check
```bash
# Check if merged_results.tsv exists
ls -lh ../../analysis_last/data/merged_results.tsv

# Check if corrected_call_rates.csv exists
ls -lh corrected_call_rates.csv
```

### Expected Output
```
-rwxrwxrwx 1 user user 6.6M Jul 26 15:35 ../../analysis_last/data/merged_results.tsv
-rwxrwxrwx 1 user user 2.6K Sep  5 18:23 corrected_call_rates.csv
```

---

## Generating merged_results.tsv from Scratch

✅ **The script is NOW INCLUDED!**

If `merged_results.tsv` is not available, you can generate it:

```bash
python3 parse_all_results.py
```

### Prerequisites:
1. **Ground truth files** must exist in:
   - `../../analysis_2/ground_truth/`
   - Files: `clingen_28012026.tsv`, `FOXL2.tsv`, etc.

2. **Tool result files** must exist in:
   - `../../Results/`
   - Subdirectories for each tool: `Intervar_20180118/`, `Genebe/`, `BIAS/`, etc.

### What the script does:
1. Loads ground truth for all 5 datasets
2. Parses output files from 14 tools using specialized parsers
3. Merges all tool classifications with ground truth
4. Saves to `../../analysis_2/results/all_variants/merged_results.tsv`
5. Also saves per-dataset files

### Included Components:
- ✅ `parse_all_results.py` - Main parsing script
- ✅ `parsers/` directory - 14 tool-specific parsers
- ✅ `utils_criteria_parser.py` - ACMG criteria utilities

---

## Regenerating corrected_call_rates.csv

If you need to regenerate `corrected_call_rates.csv`:

```bash
# Make sure merged_results.tsv exists first
python3 recalculate_call_rates.py
```

This will create a new `corrected_call_rates.csv` in:
- `../../analysis_last/results/corrected_metrics/corrected_call_rates.csv`

You can then copy it to this directory:
```bash
cp ../../analysis_last/results/corrected_metrics/corrected_call_rates.csv ./
```

---

## Summary

**Self-Contained**: ✅ `corrected_call_rates.csv` (included)  
**External Dependency**: ⚠️ `merged_results.tsv` (must exist in `../../analysis_last/data/`)

**To run Figure 2 generation**, you MUST have:
1. ✅ `corrected_call_rates.csv` (already included)
2. ⚠️ `../../analysis_last/data/merged_results.tsv` (external dependency)

If both files are available, simply run:
```bash
python3 create_figure2_comprehensive.py
```
