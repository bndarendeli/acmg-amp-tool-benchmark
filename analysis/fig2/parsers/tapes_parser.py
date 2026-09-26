from base_parser import BaseParser
import pandas as pd

class TAPESParser(BaseParser):
    """Parser for TAPES CSV output"""
    
    def __init__(self, tool_name: str = "TAPES"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse TAPES CSV file"""
        try:
            df = pd.read_csv(file_path, low_memory=False)
            
            if 'Chr' not in df.columns or 'Start' not in df.columns:
                print(f"Error: Required columns not found in TAPES file")
                return pd.DataFrame()
            
            results = []
            for idx, row in df.iterrows():
                try:
                    chrom = str(row['Chr']).replace('chr', '')
                    pos = int(row['Start'])
                    ref = str(row['Ref'])
                    alt = str(row['Alt'])
                    
                    classification = None
                    if 'Prediction_ACMG_tapes' in row and pd.notna(row['Prediction_ACMG_tapes']):
                        classification = str(row['Prediction_ACMG_tapes'])
                    
                    # Extract ACMG criteria from columns ending with _contrib
                    acmg_criteria = []
                    for col in df.columns:
                        if '_contrib' in col and pd.notna(row[col]) and row[col] == 1:
                            criterion = col.replace('_contrib', '').upper()
                            acmg_criteria.append(criterion)
                    
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
            print(f"Error parsing TAPES file {file_path}: {e}")
            return pd.DataFrame()
