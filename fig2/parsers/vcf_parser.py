import pandas as pd
from base_parser import BaseParser

class VCFParser(BaseParser):
    """Generic VCF parser for tools that output VCF format"""
    
    def __init__(self, tool_name: str = "VCF_Tool"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse generic VCF file"""
        try:
            variants = []
            
            with open(file_path, 'r') as f:
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
                    
                    classification = self._extract_classification(info)
                    criteria = self._extract_criteria(info)
                    
                    variants.append({
                        'Chr': chrom,
                        'Pos': pos,
                        'Ref': ref,
                        'Alt': alt,
                        'Classification': self.standardize_classification(classification),
                        'ACMG_Criteria': criteria,
                        'Tool': self.tool_name
                    })
            
            df = pd.DataFrame(variants)
            return self.create_variant_key(df)
            
        except Exception as e:
            print(f"Error parsing VCF file {file_path}: {e}")
            return pd.DataFrame()
    
    def _extract_classification(self, info_string: str) -> str:
        """Extract classification from INFO field"""
        for field in ['CLASSIFICATION', 'CLNSIG', 'acmg_classification_base']:
            value = self._extract_info_field(info_string, field)
            if value:
                return value
        return ''
    
    def _extract_criteria(self, info_string: str) -> str:
        """Extract ACMG criteria from INFO field"""
        for field in ['MET_CODES', 'ACMG_CRITERIA', 'acmg_criteria_base']:
            value = self._extract_info_field(info_string, field)
            if value:
                return value
        return ''
    
    def _extract_info_field(self, info_string: str, field_name: str) -> str:
        """Extract specific field from VCF INFO column"""
        for item in info_string.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                if key == field_name:
                    return value
        return ''
