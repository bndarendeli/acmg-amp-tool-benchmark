from base_parser import BaseParser
import pandas as pd

class DiabloACMGParser(BaseParser):
    """Parser for DiabloACMG TSV output"""
    
    def __init__(self, tool_name: str = "DiabloACMG"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse DiabloACMG TSV file"""
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            # DiabloACMG has very long lines, find the essential columns
            if 'chrom' not in df.columns or 'pos' not in df.columns:
                print(f"Error: Required columns not found in DiabloACMG file")
                return pd.DataFrame()
            
            # Define ACMG criterion names to look for
            acmg_criterion_names = [
                'PVS1', 'PS1', 'PS2', 'PS3', 'PS4',
                'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
                'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
                'BA1', 'BS1', 'BS2', 'BS3', 'BS4',
                'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
            ]
            
            results = []
            for idx, row in df.iterrows():
                try:
                    chrom = str(row['chrom']).replace('chr', '')
                    pos = int(row['pos'])
                    ref = str(row['ref_base'])
                    alt = str(row['alt_base'])
                    
                    # DiabloACMG stores classification in ACMG column
                    classification = None
                    if 'ACMG' in df.columns and pd.notna(row['ACMG']):
                        classification = str(row['ACMG'])
                    elif 'extra_vcf_info.Classification' in df.columns and pd.notna(row['extra_vcf_info.Classification']):
                        classification = str(row['extra_vcf_info.Classification'])
                    elif 'extra_vcf_info.CLASSIFICATION' in df.columns and pd.notna(row['extra_vcf_info.CLASSIFICATION']):
                        classification = str(row['extra_vcf_info.CLASSIFICATION'])
                    
                    # Extract ACMG criteria from individual criterion columns
                    # DiabloACMG has columns like PVS1, PS1, PS3, PM1, PM2, etc. with values 0 or 1
                    acmg_criteria = []
                    for criterion in acmg_criterion_names:
                        if criterion in df.columns and pd.notna(row[criterion]):
                            try:
                                # DiabloACMG uses 0 for not assigned, and non-zero values (1, 2, etc.) for assigned
                                # Different values may indicate strength levels
                                if int(row[criterion]) > 0:
                                    acmg_criteria.append(criterion)
                            except (ValueError, TypeError):
                                continue
                    
                    variant_key = self.make_variant_key(chrom, pos, ref, alt)
                    
                    results.append({
                        'Chr': chrom,
                        'Pos': pos,
                        'Ref': ref,
                        'Alt': alt,
                        'Variant_Key': variant_key,
                        'Classification': self.standardize_classification(classification) if classification else None,
                        'ACMG_Criteria': ','.join(acmg_criteria) if acmg_criteria else None,
                        'Tool': self.tool_name
                    })
                except Exception as e:
                    continue
            
            return pd.DataFrame(results)
            
        except Exception as e:
            print(f"Error parsing DiabloACMG file {file_path}: {e}")
            return pd.DataFrame()
