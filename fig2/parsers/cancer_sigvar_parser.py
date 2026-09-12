from base_parser import BaseParser
import pandas as pd

class CancerSIGVARParser(BaseParser):
    """Parser for Cancer SIGVAR TSV output"""
    
    def __init__(self, tool_name: str = "Cancer_SIGVAR"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse Cancer SIGVAR TSV file"""
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            if 'CHROM' not in df.columns or 'POS' not in df.columns:
                print(f"Error: Required columns not found in Cancer SIGVAR file")
                return pd.DataFrame()
            
            results = []
            for idx, row in df.iterrows():
                try:
                    chrom = str(row['CHROM']).replace('chr', '')
                    pos = int(row['POS'])
                    ref = str(row['REF'])
                    alt = str(row['ALT'])
                    
                    classification = None
                    if 'CancerSIGVAR_Automated_Interpretation' in row:
                        classification = str(row['CancerSIGVAR_Automated_Interpretation'])
                    
                    acmg_criteria = []
                    if 'Calculated_Evidences' in row and pd.notna(row['Calculated_Evidences']):
                        evidences = str(row['Calculated_Evidences'])
                        if evidences != '.' and evidences != 'NA':
                            acmg_criteria = [c.strip() for c in evidences.split('|') if c.strip()]
                    
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
            print(f"Error parsing Cancer SIGVAR file {file_path}: {e}")
            return pd.DataFrame()
