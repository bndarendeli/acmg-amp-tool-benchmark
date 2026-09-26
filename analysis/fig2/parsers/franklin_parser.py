from base_parser import BaseParser
import pandas as pd

class FranklinParser(BaseParser):
    """Parser for Franklin CSV output"""
    
    def __init__(self, tool_name: str = "Franklin"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse Franklin CSV file"""
        try:
            df = pd.read_csv(file_path, low_memory=False)
            
            if 'Chr' not in df.columns or 'Start Position' not in df.columns:
                print(f"Error: Required columns not found in Franklin file")
                return pd.DataFrame()
            
            results = []
            for idx, row in df.iterrows():
                try:
                    chrom = str(row['Chr']).replace('chr', '').replace('Chr', '')
                    pos = int(row['Start Position'])
                    ref = str(row['Ref'])
                    alt = str(row['Alt'])
                    
                    classification = None
                    if 'Classification' in row and pd.notna(row['Classification']):
                        classification = str(row['Classification'])
                    elif 'Genoox Classification' in row and pd.notna(row['Genoox Classification']):
                        classification = str(row['Genoox Classification'])
                    
                    acmg_criteria = []
                    if 'ACMG' in row and pd.notna(row['ACMG']):
                        acmg = str(row['ACMG'])
                        if acmg != '.' and acmg != 'NA' and acmg != '':
                            acmg_criteria = [c.strip() for c in acmg.split(',') if c.strip()]
                    
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
            print(f"Error parsing Franklin file {file_path}: {e}")
            return pd.DataFrame()
