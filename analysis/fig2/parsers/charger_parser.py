from base_parser import BaseParser
import pandas as pd

class CharGerParser(BaseParser):
    """Parser for CharGer (local and online) TSV output"""
    
    def __init__(self, tool_name: str = "CharGer"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse CharGer TSV file"""
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            if 'Chromosome' not in df.columns or 'Start' not in df.columns:
                print(f"Error: Required columns not found in CharGer file")
                return pd.DataFrame()
            
            results = []
            for idx, row in df.iterrows():
                try:
                    chrom = str(row['Chromosome']).replace('chr', '')
                    pos = int(row['Start'])
                    ref = str(row['Reference'])
                    alt = str(row['Alternate'])
                    
                    classification = None
                    if 'CharGer_Classification' in row and pd.notna(row['CharGer_Classification']):
                        classification = str(row['CharGer_Classification'])
                    elif 'ACMG_Classification' in row and pd.notna(row['ACMG_Classification']):
                        classification = str(row['ACMG_Classification'])
                    
                    acmg_criteria = []
                    if 'Positive_Evidence' in row and pd.notna(row['Positive_Evidence']):
                        pos_ev = str(row['Positive_Evidence'])
                        if pos_ev != 'NA' and pos_ev != '.':
                            acmg_criteria.extend([c.strip() for c in pos_ev.split(',') if c.strip()])
                    
                    if 'Negative_Evidence' in row and pd.notna(row['Negative_Evidence']):
                        neg_ev = str(row['Negative_Evidence'])
                        if neg_ev != 'NA' and neg_ev != '.':
                            acmg_criteria.extend([c.strip() for c in neg_ev.split(',') if c.strip()])
                    
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
            print(f"Error parsing CharGer file {file_path}: {e}")
            return pd.DataFrame()
