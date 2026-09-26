import pandas as pd
from base_parser import BaseParser

class GenebeParser(BaseParser):
    """Parser for GeneBe output files (.vcf format)"""
    
    def __init__(self, tool_name: str = "GeneBe"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse GeneBe VCF file"""
        try:
            variants = []
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue
                    
                    chrom = fields[0].replace('chr', '')
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    info = fields[7]
                    
                    classification = self._extract_info_field(info, 'acmg_classification_base')
                    met_codes = self._extract_info_field(info, 'acmg_criteria_base')
                    
                    variant_key = self.make_variant_key(chrom, pos, ref, alt)
                    
                    variants.append({
                        'Chr': chrom,
                        'Pos': pos,
                        'Ref': ref,
                        'Alt': alt,
                        'Variant_Key': variant_key,
                        'Classification': self.standardize_classification(classification),
                        'ACMG_Criteria': met_codes if met_codes else '',
                        'Tool': self.tool_name
                    })
            
            return pd.DataFrame(variants)
            
        except Exception as e:
            print(f"Error parsing GeneBe file {file_path}: {e}")
            return pd.DataFrame()
    
    def _extract_info_field(self, info_string: str, field_name: str) -> str:
        """Extract specific field from VCF INFO column"""
        for item in info_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                if key == field_name:
                    return value
        return ''
