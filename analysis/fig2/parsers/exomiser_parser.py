from base_parser import BaseParser
import pandas as pd

class ExomiserParser(BaseParser):
    """Parser for Exomiser TSV output"""
    
    def __init__(self, tool_name: str = "Exomiser"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse Exomiser TSV file"""
        try:
            # Exomiser files have # in header, so don't use comment parameter
            # Read the file and strip # from column names
            df = pd.read_csv(file_path, sep='\t', comment=None, low_memory=False)
            df.columns = df.columns.str.lstrip('#')
            
            if 'CONTIG' not in df.columns or 'START' not in df.columns:
                print(f"Error: Required columns not found in Exomiser file")
                print(f"Available columns: {df.columns.tolist()[:10]}")
                return pd.DataFrame()
            
            results = []
            for idx, row in df.iterrows():
                try:
                    chrom = str(row['CONTIG']).replace('chr', '')
                    pos = int(row['START'])
                    ref = str(row['REF'])
                    alt = str(row['ALT'])
                    
                    classification = None
                    if 'EXOMISER_ACMG_CLASSIFICATION' in row and pd.notna(row['EXOMISER_ACMG_CLASSIFICATION']):
                        classification = str(row['EXOMISER_ACMG_CLASSIFICATION'])
                    
                    acmg_criteria = []
                    if 'EXOMISER_ACMG_EVIDENCE' in row and pd.notna(row['EXOMISER_ACMG_EVIDENCE']):
                        evidence = str(row['EXOMISER_ACMG_EVIDENCE'])
                        if evidence != '.' and evidence != 'NA' and evidence != '':
                            acmg_criteria = [c.strip() for c in evidence.split(',') if c.strip()]
                    
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
            print(f"Error parsing Exomiser file {file_path}: {e}")
            return pd.DataFrame()
