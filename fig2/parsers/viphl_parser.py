from base_parser import BaseParser
import pandas as pd

class VIPHLParser(BaseParser):
    """Parser for VIP-HL TSV output"""
    
    def __init__(self, tool_name: str = "VIP-HL"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse VIP-HL TSV file
        
        Format: Variant_Input contains chr-pos-ref-alt (hg19 coordinates)
        Example: 1-100000-A-G
        """
        try:
            df = pd.read_csv(file_path, sep='\t')
            
            variants = []
            for idx, row in df.iterrows():
                # Parse Variant_Input: format is "Chr-Pos-Ref-Alt" (hg19)
                variant_input = str(row.get('Variant_Input', ''))
                if not variant_input or variant_input == 'nan':
                    continue
                
                # Split by dash
                parts = variant_input.split('-')
                if len(parts) < 4:
                    continue
                
                # Handle cases where ref/alt might contain dashes
                chrom = parts[0]
                pos = parts[1]
                # Join remaining parts and split again to get ref and alt
                remaining = '-'.join(parts[2:])
                # Try to split ref and alt (usually last dash)
                ref_alt_parts = remaining.rsplit('-', 1)
                if len(ref_alt_parts) == 2:
                    ref, alt = ref_alt_parts
                else:
                    # If can't split, use parts directly
                    ref = parts[2] if len(parts) > 2 else ''
                    alt = parts[3] if len(parts) > 3 else ''
                
                # Clean chromosome
                chrom = chrom.replace('chr', '').replace('Chr', '')
                
                # Create variant key (this is hg19)
                variant_key = self.make_variant_key(chrom, pos, ref, alt)
                
                # Get classification and criteria
                classification = self.standardize_classification(row.get('Classification', ''))
                criteria = row.get('Criteria', '')
                
                # Parse criteria
                if pd.notna(criteria) and criteria != '':
                    criteria_list = [c.strip() for c in str(criteria).split(',')]
                    criteria = ','.join([c for c in criteria_list if c and c != 'nan' and c != ''])
                else:
                    criteria = ''
                
                variants.append({
                    'Chr': chrom,
                    'Pos': int(pos),
                    'Ref': ref,
                    'Alt': alt,
                    'Variant_Key': variant_key,
                    'Classification': classification,
                    'ACMG_Criteria': criteria,
                    'Tool': self.tool_name
                })
            
            return pd.DataFrame(variants)
            
        except Exception as e:
            print(f"Error parsing VIP-HL file {file_path}: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
