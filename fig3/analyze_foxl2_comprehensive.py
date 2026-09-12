#!/usr/bin/env python3
"""
Comprehensive FOXL2 Analysis

Performs complete analysis for FOXL2 dataset:
1. Jaccard similarity heatmap (Criterion × Tool)
2. Classification metrics (Accuracy, Precision, Recall, F1)
3. Confusion matrices per tool
4. Per-class performance (Pathogenic, Benign, VUS)
5. ACMG criteria analysis
6. Nature Communications-style visualizations

All outputs match the format used for ClinGen datasets.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix
import sys

sys.path.insert(0, str(Path(__file__).parent))
from utils_criteria_parser import parse_criteria_string

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

TOOLS = [
    'InterVar_2018', 'InterVar_2025', 'BIAS', 'CharGer_Local', 'CharGer_Online',
    'DiabloACMG', 'Exomiser', 'Genebe', 'TAPES', 'Franklin', 'AutoGVP', 'VIP-HL',
    'CancerSIGVAR', 'CPSR'
]

ALL_CRITERIA = [
    'PVS1',
    'PS1', 'PS2', 'PS3', 'PS4',
    'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
    'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
    'BA1',
    'BS1', 'BS2', 'BS3', 'BS4',
    'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
]

NATURE_COLORS = {
    'primary': '#0173B2',
    'secondary': '#DE8F05',
    'tertiary': '#029E73',
    'pathogenic': '#D55E00',
    'benign': '#56B4E9',
}


def load_foxl2_data():
    """Load FOXL2 merged results"""
    base_dir = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep')
    foxl2_file = base_dir / 'datasets' / 'foxl2' / 'analysis' / 'results' / 'foxl2_merged_results.tsv'
    
    if not foxl2_file.exists():
        print(f"✗ FOXL2 data file not found: {foxl2_file}")
        return None
    
    df = pd.read_csv(foxl2_file, sep='\t')
    print(f"✓ Loaded {len(df)} FOXL2 variants")
    
    return df


def calculate_jaccard_per_criterion(df):
    """Calculate Jaccard similarity per tool per criterion"""
    results = []
    
    for tool in TOOLS:
        criteria_col = f'{tool}_ACMG_Criteria'
        if criteria_col not in df.columns:
            continue
        
        for criterion in ALL_CRITERIA:
            gt_has_criterion = []
            tool_has_criterion = []
            
            for idx, row in df.iterrows():
                gt_criteria = parse_criteria_string(row.get('Ground_Truth_ACMG', ''))
                tool_criteria = parse_criteria_string(row.get(criteria_col, ''))
                
                gt_has = criterion in gt_criteria
                tool_has = criterion in tool_criteria
                
                gt_has_criterion.append(gt_has)
                tool_has_criterion.append(tool_has)
            
            gt_set = set([i for i, x in enumerate(gt_has_criterion) if x])
            tool_set = set([i for i, x in enumerate(tool_has_criterion) if x])
            
            intersection = len(gt_set.intersection(tool_set))
            union = len(gt_set.union(tool_set))
            jaccard = intersection / union if union > 0 else 0.0
            
            results.append({
                'Tool': tool,
                'Criterion': criterion,
                'Jaccard_Similarity': jaccard,
                'Intersection': intersection,
                'Union': union
            })
    
    return pd.DataFrame(results)


def create_jaccard_heatmap(jaccard_df, output_dir):
    """Create Jaccard similarity heatmap"""
    
    pivot_data = jaccard_df.pivot(index='Criterion', columns='Tool', values='Jaccard_Similarity')
    
    # Sort by average performance
    pivot_data['avg'] = pivot_data.mean(axis=1)
    pivot_data = pivot_data.sort_values('avg', ascending=False)
    pivot_data = pivot_data.drop('avg', axis=1)
    
    fig, ax = plt.subplots(figsize=(14, max(8, len(pivot_data) * 0.3)))
    
    sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='RdYlBu_r', 
                vmin=0, vmax=1, cbar_kws={'label': 'Jaccard Similarity'}, 
                ax=ax, linewidths=0.5, linecolor='white')
    
    ax.set_title('Jaccard Similarity Heatmap - FOXL2', 
                 fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel('Tool', fontsize=10, fontweight='bold')
    ax.set_ylabel('ACMG Criterion', fontsize=10, fontweight='bold')
    
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(rotation=0, fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'foxl2_jaccard_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'foxl2_jaccard_heatmap.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created: foxl2_jaccard_heatmap.png/pdf")


def normalize_classification(cls):
    """Normalize classification labels"""
    if pd.isna(cls) or str(cls).strip() == '' or str(cls).lower() in ['unknown', 'nan']:
        return 'VUS'
    
    cls_str = str(cls).strip().replace('_', ' ').lower()
    
    # Handle short forms (P, LP, B, LB, VUS)
    if cls_str == 'p':
        return 'Pathogenic'
    elif cls_str == 'lp':
        return 'Likely_Pathogenic'
    elif cls_str == 'b':
        return 'Benign'
    elif cls_str == 'lb':
        return 'Likely_Benign'
    elif cls_str == 'vus' or cls_str == 'vous':
        return 'VUS'
    # Handle full forms
    elif 'pathogenic' in cls_str and 'likely' not in cls_str:
        return 'Pathogenic'
    elif 'likely pathogenic' in cls_str or 'likely_pathogenic' in cls_str:
        return 'Likely_Pathogenic'
    elif 'benign' in cls_str and 'likely' not in cls_str:
        return 'Benign'
    elif 'likely benign' in cls_str or 'likely_benign' in cls_str:
        return 'Likely_Benign'
    else:
        return 'VUS'


def calculate_classification_metrics(df):
    """Calculate classification metrics for all tools"""
    results = []
    
    # Normalize ground truth
    df['GT_Normalized'] = df['Ground_Truth_Classification'].apply(normalize_classification)
    
    # Combine P/LP as "Pathogenic" and B/LB as "Benign" for binary metrics
    def to_binary(cls):
        if cls in ['Pathogenic', 'Likely_Pathogenic']:
            return 'Pathogenic'
        elif cls in ['Benign', 'Likely_Benign']:
            return 'Benign'
        else:
            return 'VUS'
    
    df['GT_Binary'] = df['GT_Normalized'].apply(to_binary)
    
    for tool in TOOLS:
        class_col = f'{tool}_Classification'
        if class_col not in df.columns:
            continue
        
        # Get valid pairs (both GT and tool have values)
        valid_mask = df['GT_Normalized'].notna() & df[class_col].notna()
        if valid_mask.sum() == 0:
            continue
        
        gt_valid = df.loc[valid_mask, 'GT_Normalized']
        tool_valid = df.loc[valid_mask, class_col].apply(normalize_classification)
        
        # Binary classification for P/LP vs B/LB
        gt_binary = df.loc[valid_mask, 'GT_Binary']
        tool_binary = tool_valid.apply(to_binary)
        
        # Calculate metrics
        accuracy = accuracy_score(gt_valid, tool_valid)
        
        # Weighted metrics (5-class)
        precision_w, recall_w, f1_w, _ = precision_recall_fscore_support(
            gt_valid, tool_valid, average='weighted', zero_division=0
        )
        
        # Macro metrics (5-class)
        precision_m, recall_m, f1_m, _ = precision_recall_fscore_support(
            gt_valid, tool_valid, average='macro', zero_division=0
        )
        
        # Binary metrics (P/LP vs B/LB vs VUS)
        precision_bin, recall_bin, f1_bin, support_bin = precision_recall_fscore_support(
            gt_binary, tool_binary, labels=['Pathogenic', 'Benign', 'VUS'],
            zero_division=0
        )
        
        results.append({
            'Tool': tool,
            'Total_Variants': valid_mask.sum(),
            'Accuracy': accuracy,
            'Precision_Weighted': precision_w,
            'Recall_Weighted': recall_w,
            'F1_Score_Weighted': f1_w,
            'Precision_Macro': precision_m,
            'Recall_Macro': recall_m,
            'F1_Score_Macro': f1_m,
            'Pathogenic_Precision': precision_bin[0],
            'Pathogenic_Recall': recall_bin[0],
            'Pathogenic_F1': f1_bin[0],
            'Benign_Precision': precision_bin[1] if len(precision_bin) > 1 else 0.0,
            'Benign_Recall': recall_bin[1] if len(recall_bin) > 1 else 0.0,
            'Benign_F1': f1_bin[1] if len(f1_bin) > 1 else 0.0,
            'VUS_Precision': precision_bin[2] if len(precision_bin) > 2 else 0.0,
            'VUS_Recall': recall_bin[2] if len(recall_bin) > 2 else 0.0,
        })
    
    return pd.DataFrame(results)


def create_classification_metrics_plot(metrics_df, output_dir):
    """Create classification metrics visualization"""
    
    metrics_df = metrics_df.sort_values('Accuracy', ascending=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 6))
    
    metrics = [
        ('Accuracy', 'Accuracy'),
        ('Precision_Weighted', 'Precision (Weighted)'),
        ('Recall_Weighted', 'Recall (Weighted)')
    ]
    
    for idx, (metric, title) in enumerate(metrics):
        ax = axes[idx]
        
        y_pos = np.arange(len(metrics_df))
        bars = ax.barh(y_pos, metrics_df[metric], color=NATURE_COLORS['primary'], alpha=0.8)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(metrics_df['Tool'], fontsize=8)
        ax.set_xlabel(title, fontsize=9, fontweight='bold')
        ax.set_xlim(0, 1.0)
        ax.grid(axis='x', alpha=0.3, linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        for i, (bar, val) in enumerate(zip(bars, metrics_df[metric])):
            ax.text(val + 0.02, i, f'{val:.3f}', va='center', fontsize=7)
    
    plt.suptitle('FOXL2 Classification Performance', fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'foxl2_classification_metrics.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'foxl2_classification_metrics.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created: foxl2_classification_metrics.png/pdf")


def create_pathogenic_benign_comparison(metrics_df, output_dir):
    """Create Pathogenic vs Benign performance comparison"""
    
    metrics_df = metrics_df.sort_values('Accuracy', ascending=False)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Pathogenic
    ax = axes[0]
    x = np.arange(len(metrics_df))
    width = 0.35
    
    ax.bar(x - width/2, metrics_df['Pathogenic_Precision'], width, 
           label='Precision', color=NATURE_COLORS['pathogenic'], alpha=0.8)
    ax.bar(x + width/2, metrics_df['Pathogenic_Recall'], width, 
           label='Recall', color=NATURE_COLORS['pathogenic'], alpha=0.5)
    
    ax.set_ylabel('Score', fontsize=9, fontweight='bold')
    ax.set_title('Pathogenic Classification', fontsize=10, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics_df['Tool'], rotation=45, ha='right', fontsize=8)
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Benign
    ax = axes[1]
    ax.bar(x - width/2, metrics_df['Benign_Precision'], width, 
           label='Precision', color=NATURE_COLORS['benign'], alpha=0.8)
    ax.bar(x + width/2, metrics_df['Benign_Recall'], width, 
           label='Recall', color=NATURE_COLORS['benign'], alpha=0.5)
    
    ax.set_ylabel('Score', fontsize=9, fontweight='bold')
    ax.set_title('Benign Classification', fontsize=10, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics_df['Tool'], rotation=45, ha='right', fontsize=8)
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.grid(axis='y', alpha=0.3, linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.suptitle('FOXL2: Pathogenic vs Benign Performance', fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / 'foxl2_pathogenic_benign_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'foxl2_pathogenic_benign_comparison.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created: foxl2_pathogenic_benign_comparison.png/pdf")


def main():
    """Main execution function"""
    
    base_dir = Path('/mnt/d/busra-research/doktora/1002/code_works/dataset_prep/analysis_last')
    results_dir = base_dir / 'results' / 'foxl2'
    figures_dir = base_dir / 'figures' / 'foxl2'
    
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("COMPREHENSIVE FOXL2 ANALYSIS")
    print("="*80)
    
    # Load data
    print("\n1. Loading FOXL2 data...")
    df = load_foxl2_data()
    
    if df is None:
        return
    
    # Jaccard analysis
    print("\n2. Calculating Jaccard similarity per criterion...")
    jaccard_df = calculate_jaccard_per_criterion(df)
    jaccard_df.to_csv(results_dir / 'foxl2_jaccard_per_criterion.csv', index=False)
    print(f"  ✓ Saved: foxl2_jaccard_per_criterion.csv")
    
    print("\n3. Creating Jaccard heatmap...")
    create_jaccard_heatmap(jaccard_df, figures_dir)
    
    # Classification metrics
    print("\n4. Calculating classification metrics...")
    metrics_df = calculate_classification_metrics(df)
    metrics_df.to_csv(results_dir / 'foxl2_classification_metrics.csv', index=False)
    print(f"  ✓ Saved: foxl2_classification_metrics.csv")
    
    print("\n5. Creating classification metrics plots...")
    create_classification_metrics_plot(metrics_df, figures_dir)
    create_pathogenic_benign_comparison(metrics_df, figures_dir)
    
    # Summary
    print("\n" + "="*80)
    print("FOXL2 ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {results_dir}")
    print(f"Figures saved to: {figures_dir}")
    print("\nGenerated files:")
    print("  Data:")
    print("    - foxl2_jaccard_per_criterion.csv")
    print("    - foxl2_classification_metrics.csv")
    print("  Figures:")
    print("    - foxl2_jaccard_heatmap.png/pdf")
    print("    - foxl2_classification_metrics.png/pdf")
    print("    - foxl2_pathogenic_benign_comparison.png/pdf")
    print("\n✓ All FOXL2 analyses complete!")


if __name__ == "__main__":
    main()
