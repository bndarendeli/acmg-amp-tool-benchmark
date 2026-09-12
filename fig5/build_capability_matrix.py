#!/usr/bin/env python3
"""
Build Tool × ACMG/AMP Criterion Capability Matrix

Uses the authoritative source file (reference/tool_acmg_capabilities.csv) to
construct a comprehensive capability matrix for Figure 5 capability-adjusted
evidence analysis.

Includes both documented and operational status based on actual benchmark inputs.
"""

import pandas as pd
import re
from pathlib import Path
from collections import defaultdict

# Configuration
BASE_DIR = Path(__file__).parent
REFERENCE_FILE = BASE_DIR / 'reference' / 'tool_acmg_capabilities.csv'
OUTPUT_DIR = BASE_DIR / 'outputs'
MATRIX_FILE = OUTPUT_DIR / 'tool_criterion_capability_matrix.csv'
AUDIT_FILE = OUTPUT_DIR / 'tool_criterion_operational_audit.txt'

# 28 canonical ACMG/AMP criteria
CANONICAL_CRITERIA = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]

# Tool-specific codes to EXCLUDE (not mapped to standard ACMG)
TOOL_SPECIFIC_CODES = ['PMC1', 'BSC1', 'PPC1', 'PPC2', 'BMC1', 'PSC1']

# Tools in Figure 5 benchmark
FIG4_TOOLS = [
    'AutoGVP', 'BIAS', 'CPSR', 'CancerSIGVAR', 'CharGer_Local', 'CharGer_Online',
    'DiabloACMG', 'Exomiser', 'Franklin', 'Genebe', 'InterVar_2018', 
    'InterVar_2025', 'TAPES', 'VIP-HL'
]


def parse_documented_criteria(criteria_str):
    """
    Parse documented criteria string and categorize by status.
    
    Returns:
        dict: {
            'automated': set of criteria,
            'conditional_parental_vcf': set of criteria (marked with *),
            'conditional_phased_vcf': set of criteria (marked with +),
            'semi_automated': set of criteria (marked with -)
        }
    """
    if pd.isna(criteria_str) or not criteria_str:
        return {
            'automated': set(),
            'conditional_parental_vcf': set(),
            'conditional_phased_vcf': set(),
            'semi_automated': set()
        }
    
    categorized = {
        'automated': set(),
        'conditional_parental_vcf': set(),
        'conditional_phased_vcf': set(),
        'semi_automated': set()
    }
    
    # Split by comma
    criteria_list = [c.strip() for c in str(criteria_str).split(',')]
    
    for criterion in criteria_list:
        if not criterion:
            continue
        
        # Check for markers
        if '*' in criterion:
            # Parental VCF required
            base = criterion.replace('*', '').strip()
            if base in CANONICAL_CRITERIA:
                categorized['conditional_parental_vcf'].add(base)
        elif '+' in criterion:
            # Phased VCF required
            base = criterion.replace('+', '').strip()
            if base in CANONICAL_CRITERIA:
                categorized['conditional_phased_vcf'].add(base)
        elif '-' in criterion:
            # Semi-automated (manual input required)
            base = criterion.replace('-', '').strip()
            if base in CANONICAL_CRITERIA:
                categorized['semi_automated'].add(base)
        else:
            # Automated (no special requirements)
            if criterion in CANONICAL_CRITERIA:
                categorized['automated'].add(criterion)
            # Silently skip tool-specific codes (PMC1, BSC1, etc.)
    
    return categorized


def determine_operational_status(documented_status, benchmark_inputs):
    """
    Determine operational status based on documented status and actual benchmark inputs.
    
    Args:
        documented_status: One of 'automated', 'conditional_parental_vcf', 
                          'conditional_phased_vcf', 'semi_automated', 'unsupported'
        benchmark_inputs: dict with keys 'parental_vcf', 'phased_vcf', 'manual_input'
    
    Returns:
        str: 'supported' or 'unsupported'
    """
    if documented_status == 'unsupported':
        return 'unsupported'
    
    if documented_status == 'automated':
        return 'supported'
    
    if documented_status == 'conditional_parental_vcf':
        return 'supported' if benchmark_inputs['parental_vcf'] else 'unsupported'
    
    if documented_status == 'conditional_phased_vcf':
        return 'supported' if benchmark_inputs['phased_vcf'] else 'unsupported'
    
    if documented_status == 'semi_automated':
        return 'supported' if benchmark_inputs['manual_input'] else 'unsupported'
    
    return 'unsupported'


def audit_benchmark_inputs():
    """
    Audit the benchmark workflow to determine what inputs were actually provided.
    
    Based on pipeline documentation and scripts, determine:
    - Were parental VCF files provided?
    - Were phased VCF files provided?
    - Were phenotype/family history data provided?
    - Were manual quantitative inputs provided?
    
    Returns:
        dict: Benchmark input availability
    """
    # Based on review of:
    # - acmg_analysis_pipeline/README.md
    # - Pipeline configuration files
    # - Tool execution scripts
    
    # The benchmark workflow only provided:
    # - Single-sample VCF files with variant coordinates
    # - No parental VCFs
    # - No phased VCFs
    # - No phenotype/family history data
    # - No manual quantitative inputs
    
    benchmark_inputs = {
        'parental_vcf': False,
        'phased_vcf': False,
        'phenotype_data': False,
        'family_history': False,
        'manual_input': False
    }
    
    return benchmark_inputs


def build_capability_matrix():
    """
    Build the tool × criterion capability matrix.
    
    Returns:
        pd.DataFrame: Capability matrix
    """
    print("="*80)
    print("BUILDING TOOL × ACMG/AMP CRITERION CAPABILITY MATRIX")
    print("="*80)
    
    # Load reference file
    print(f"\nLoading reference: {REFERENCE_FILE}")
    df_ref = pd.read_csv(REFERENCE_FILE)
    print(f"✓ Loaded {len(df_ref)} tool capability definitions")
    
    # Audit benchmark inputs
    print("\nAuditing benchmark inputs...")
    benchmark_inputs = audit_benchmark_inputs()
    print(f"  Parental VCF provided: {benchmark_inputs['parental_vcf']}")
    print(f"  Phased VCF provided: {benchmark_inputs['phased_vcf']}")
    print(f"  Phenotype data provided: {benchmark_inputs['phenotype_data']}")
    print(f"  Family history provided: {benchmark_inputs['family_history']}")
    print(f"  Manual input provided: {benchmark_inputs['manual_input']}")
    
    # Build matrix
    print("\nBuilding capability matrix...")
    matrix_rows = []
    
    for _, row in df_ref.iterrows():
        tool_name = row['tool']
        criteria_str = row['implemented_criteria']
        
        # Handle InterVar and CharGer variants
        if tool_name == 'InterVar':
            tool_variants = ['InterVar_2018', 'InterVar_2025']
        elif tool_name == 'CharGer':
            tool_variants = ['CharGer_Local', 'CharGer_Online']
        else:
            tool_variants = [tool_name]
        
        # Parse documented criteria
        categorized = parse_documented_criteria(criteria_str)
        
        # For each tool variant
        for tool in tool_variants:
            # Skip if not in Figure 5 benchmark
            if tool not in FIG4_TOOLS:
                continue
            
            # For each canonical criterion
            for criterion in CANONICAL_CRITERIA:
                # Determine documented status
                if criterion in categorized['automated']:
                    documented_status = 'automated'
                    requirement = 'none'
                    source_note = 'Listed in reference without markers'
                elif criterion in categorized['conditional_parental_vcf']:
                    documented_status = 'conditional_parental_vcf'
                    requirement = 'parental_vcf'
                    source_note = 'Marked with * (requires parental VCF)'
                elif criterion in categorized['conditional_phased_vcf']:
                    documented_status = 'conditional_phased_vcf'
                    requirement = 'phased_vcf'
                    source_note = 'Marked with + (requires phased VCF)'
                elif criterion in categorized['semi_automated']:
                    documented_status = 'semi_automated'
                    requirement = 'manual_input'
                    source_note = 'Marked with - (requires manual quantitative input)'
                else:
                    documented_status = 'unsupported'
                    requirement = 'n/a'
                    source_note = 'Not listed in reference'
                
                # Determine operational status
                operational_status = determine_operational_status(
                    documented_status, benchmark_inputs
                )
                
                matrix_rows.append({
                    'tool': tool,
                    'criterion': criterion,
                    'documented_status': documented_status,
                    'operational_status': operational_status,
                    'requirement': requirement,
                    'source_note': source_note
                })
    
    df_matrix = pd.DataFrame(matrix_rows)
    
    print(f"✓ Built matrix: {len(df_matrix)} rows")
    print(f"  Tools: {df_matrix['tool'].nunique()}")
    print(f"  Criteria: {df_matrix['criterion'].nunique()}")
    
    return df_matrix, benchmark_inputs


def generate_audit_report(df_matrix, benchmark_inputs):
    """
    Generate operational audit report.
    
    Args:
        df_matrix: Capability matrix DataFrame
        benchmark_inputs: Benchmark input availability dict
    
    Returns:
        str: Audit report text
    """
    lines = []
    lines.append("="*80)
    lines.append("TOOL × ACMG/AMP CRITERION OPERATIONAL AUDIT")
    lines.append("="*80)
    lines.append("")
    
    # Section 1: Benchmark Input Audit
    lines.append("SECTION 1: BENCHMARK INPUT AUDIT")
    lines.append("-"*80)
    lines.append("")
    lines.append("Based on review of:")
    lines.append("  - acmg_analysis_pipeline/README.md")
    lines.append("  - Pipeline configuration files (conf/datasets_config.json)")
    lines.append("  - Tool execution scripts")
    lines.append("  - Input directory structure documentation")
    lines.append("")
    lines.append("Benchmark Inputs Provided:")
    lines.append(f"  ✗ Parental VCF files: {benchmark_inputs['parental_vcf']}")
    lines.append(f"  ✗ Phased VCF files: {benchmark_inputs['phased_vcf']}")
    lines.append(f"  ✗ Phenotype data: {benchmark_inputs['phenotype_data']}")
    lines.append(f"  ✗ Family history data: {benchmark_inputs['family_history']}")
    lines.append(f"  ✗ Manual quantitative inputs: {benchmark_inputs['manual_input']}")
    lines.append("")
    lines.append("Conclusion:")
    lines.append("  The benchmark workflow provided ONLY single-sample VCF files with")
    lines.append("  variant coordinates. No additional inputs (parental VCFs, phased VCFs,")
    lines.append("  phenotype data, family history, or manual inputs) were supplied.")
    lines.append("")
    
    # Section 2: Operational Status Summary
    lines.append("")
    lines.append("SECTION 2: OPERATIONAL STATUS SUMMARY")
    lines.append("-"*80)
    lines.append("")
    
    # Count by status
    status_counts = df_matrix.groupby(['documented_status', 'operational_status']).size().reset_index(name='count')
    
    lines.append("Documented Status → Operational Status:")
    lines.append("")
    for _, row in status_counts.iterrows():
        doc_status = row['documented_status']
        op_status = row['operational_status']
        count = row['count']
        lines.append(f"  {doc_status:30s} → {op_status:15s} : {count:4d} tool×criterion pairs")
    
    lines.append("")
    
    # Total operational support
    total_pairs = len(df_matrix)
    operational_supported = len(df_matrix[df_matrix['operational_status'] == 'supported'])
    operational_unsupported = len(df_matrix[df_matrix['operational_status'] == 'unsupported'])
    
    lines.append(f"Total tool×criterion pairs: {total_pairs}")
    lines.append(f"  Operationally supported: {operational_supported} ({operational_supported/total_pairs*100:.1f}%)")
    lines.append(f"  Operationally unsupported: {operational_unsupported} ({operational_unsupported/total_pairs*100:.1f}%)")
    lines.append("")
    
    # Section 3: Criteria Requiring Special Inputs
    lines.append("")
    lines.append("SECTION 3: CRITERIA REQUIRING SPECIAL INPUTS (NOT OPERATIONAL)")
    lines.append("-"*80)
    lines.append("")
    
    # Parental VCF required
    parental_vcf_criteria = df_matrix[
        df_matrix['documented_status'] == 'conditional_parental_vcf'
    ].groupby('criterion')['tool'].apply(list).reset_index()
    
    if len(parental_vcf_criteria) > 0:
        lines.append("Criteria Requiring Parental VCF (marked with *):")
        lines.append("")
        for _, row in parental_vcf_criteria.iterrows():
            criterion = row['criterion']
            tools = row['tool']
            lines.append(f"  {criterion}:")
            for tool in tools:
                lines.append(f"    - {tool}")
        lines.append("")
        lines.append(f"  Status: NOT OPERATIONAL (parental VCFs not provided)")
        lines.append("")
    else:
        lines.append("No criteria require parental VCF.")
        lines.append("")
    
    # Phased VCF required
    phased_vcf_criteria = df_matrix[
        df_matrix['documented_status'] == 'conditional_phased_vcf'
    ].groupby('criterion')['tool'].apply(list).reset_index()
    
    if len(phased_vcf_criteria) > 0:
        lines.append("Criteria Requiring Phased VCF (marked with +):")
        lines.append("")
        for _, row in phased_vcf_criteria.iterrows():
            criterion = row['criterion']
            tools = row['tool']
            lines.append(f"  {criterion}:")
            for tool in tools:
                lines.append(f"    - {tool}")
        lines.append("")
        lines.append(f"  Status: NOT OPERATIONAL (phased VCFs not provided)")
        lines.append("")
    else:
        lines.append("No criteria require phased VCF.")
        lines.append("")
    
    # Manual input required
    manual_input_criteria = df_matrix[
        df_matrix['documented_status'] == 'semi_automated'
    ].groupby('criterion')['tool'].apply(list).reset_index()
    
    if len(manual_input_criteria) > 0:
        lines.append("Criteria Requiring Manual Input (marked with -):")
        lines.append("")
        for _, row in manual_input_criteria.iterrows():
            criterion = row['criterion']
            tools = row['tool']
            lines.append(f"  {criterion}:")
            for tool in tools:
                lines.append(f"    - {tool}")
        lines.append("")
        lines.append(f"  Status: NOT OPERATIONAL (manual inputs not provided)")
        lines.append("")
    else:
        lines.append("No criteria require manual input.")
        lines.append("")
    
    # Section 4: Tool-Specific Summary
    lines.append("")
    lines.append("SECTION 4: TOOL-SPECIFIC CAPABILITY SUMMARY")
    lines.append("-"*80)
    lines.append("")
    
    for tool in sorted(df_matrix['tool'].unique()):
        df_tool = df_matrix[df_matrix['tool'] == tool]
        
        n_documented = len(df_tool[df_tool['documented_status'] != 'unsupported'])
        n_operational = len(df_tool[df_tool['operational_status'] == 'supported'])
        n_automated = len(df_tool[df_tool['documented_status'] == 'automated'])
        n_conditional = len(df_tool[df_tool['documented_status'].str.contains('conditional')])
        n_semi = len(df_tool[df_tool['documented_status'] == 'semi_automated'])
        
        lines.append(f"{tool}:")
        lines.append(f"  Documented criteria: {n_documented}/28")
        lines.append(f"    Automated: {n_automated}")
        lines.append(f"    Conditional (parental/phased VCF): {n_conditional}")
        lines.append(f"    Semi-automated (manual input): {n_semi}")
        lines.append(f"  Operationally supported: {n_operational}/28")
        lines.append(f"  Gap (documented but not operational): {n_documented - n_operational}")
        lines.append("")
    
    # Section 5: Recommendations
    lines.append("")
    lines.append("SECTION 5: RECOMMENDATIONS FOR CAPABILITY-ADJUSTED ANALYSIS")
    lines.append("-"*80)
    lines.append("")
    lines.append("1. Use 'operational_status' for capability-adjusted evidence analysis")
    lines.append("   - Only count criteria where operational_status = 'supported'")
    lines.append("   - Exclude criteria requiring inputs not provided in benchmark")
    lines.append("")
    lines.append("2. Document the gap between documented and operational capability")
    lines.append("   - Many tools document criteria requiring special inputs")
    lines.append("   - These are not operational in our benchmark workflow")
    lines.append("")
    lines.append("3. Consider separate analysis of automated vs. conditional criteria")
    lines.append("   - Automated criteria: fair comparison across all tools")
    lines.append("   - Conditional criteria: would require additional benchmark inputs")
    lines.append("")
    lines.append("4. Note tool-specific codes were excluded per normalization audit")
    lines.append("   - PMC1, BSC1, PPC1, PPC2, BMC1, PSC1 not mapped to standard ACMG")
    lines.append("   - Preserves consistency with Figure 2 and Figure 5 normalization")
    lines.append("")
    
    lines.append("="*80)
    lines.append("END OF OPERATIONAL AUDIT")
    lines.append("="*80)
    
    return '\n'.join(lines)


def main():
    """Main execution"""
    print("="*80)
    print("TOOL × ACMG/AMP CRITERION CAPABILITY MATRIX BUILDER")
    print("="*80)
    
    # Build matrix
    df_matrix, benchmark_inputs = build_capability_matrix()
    
    # Generate audit report
    print("\nGenerating operational audit report...")
    audit_report = generate_audit_report(df_matrix, benchmark_inputs)
    
    # Save outputs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\nSaving capability matrix to: {MATRIX_FILE}")
    df_matrix.to_csv(MATRIX_FILE, index=False)
    print(f"✓ Saved {len(df_matrix)} rows")
    
    print(f"\nSaving operational audit to: {AUDIT_FILE}")
    with open(AUDIT_FILE, 'w') as f:
        f.write(audit_report)
    print(f"✓ Saved audit report")
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    print(f"\nCapability Matrix:")
    print(f"  Tools: {df_matrix['tool'].nunique()}")
    print(f"  Criteria: {df_matrix['criterion'].nunique()}")
    print(f"  Total pairs: {len(df_matrix)}")
    
    print(f"\nDocumented Status:")
    for status in df_matrix['documented_status'].unique():
        count = len(df_matrix[df_matrix['documented_status'] == status])
        print(f"  {status:30s}: {count:4d}")
    
    print(f"\nOperational Status:")
    for status in df_matrix['operational_status'].unique():
        count = len(df_matrix[df_matrix['operational_status'] == status])
        pct = count / len(df_matrix) * 100
        print(f"  {status:15s}: {count:4d} ({pct:5.1f}%)")
    
    print("\n" + "="*80)
    print("✓ CAPABILITY MATRIX BUILD COMPLETE")
    print("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
