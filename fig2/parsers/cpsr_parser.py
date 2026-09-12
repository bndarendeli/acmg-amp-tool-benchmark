from base_parser import BaseParser
import pandas as pd

class CPSRParser(BaseParser):
    """Parser for CPSR TSV output"""
    
    def __init__(self, tool_name: str = "CPSR"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse CPSR TSV file"""
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            if 'GENOMIC_CHANGE' not in df.columns:
                print(f"Error: GENOMIC_CHANGE column not found in CPSR file")
                return pd.DataFrame()
            
            results = []
            for idx, row in df.iterrows():
                try:
                    # Parse genomic change format: Chr:g.Pos_Ref>Alt
                    genomic_change = str(row['GENOMIC_CHANGE'])
                    # Format example: 21:g.34859477G>A
                    parts = genomic_change.split(':g.')
                    if len(parts) != 2:
                        continue
                    
                    chrom = parts[0].replace('chr', '')
                    pos_alleles = parts[1]
                    
                    # Extract position and alleles
                    for i, char in enumerate(pos_alleles):
                        if char.isalpha():
                            pos = int(pos_alleles[:i])
                            alleles = pos_alleles[i:]
                            break
                    
                    if '>' in alleles:
                        ref, alt = alleles.split('>')
                    else:
                        continue
                    
                    classification = None
                    if 'CPSR_CLASSIFICATION' in row and pd.notna(row['CPSR_CLASSIFICATION']):
                        classification = str(row['CPSR_CLASSIFICATION'])
                    elif 'FINAL_CLASSIFICATION' in row and pd.notna(row['FINAL_CLASSIFICATION']):
                        classification = str(row['FINAL_CLASSIFICATION'])
                    
                    acmg_criteria = []
                    if 'CPSR_CLASSIFICATION_CODE' in row and pd.notna(row['CPSR_CLASSIFICATION_CODE']):
                        codes = str(row['CPSR_CLASSIFICATION_CODE'])
                        if codes != '.' and codes != 'NA':
                            acmg_criteria = [c.strip() for c in codes.split('|') if c.strip() and c.strip().startswith('ACMG_')]
                            acmg_criteria = [c.replace('ACMG_', '').split('_')[0] for c in acmg_criteria]
                    
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
            print(f"Error parsing CPSR file {file_path}: {e}")
            return pd.DataFrame()
