import pandas as pd
import json
from base_parser import BaseParser

class BIASParser(BaseParser):
    """Parser for BIAS output files (.tsv format)"""
    
    def __init__(self, tool_name: str = "BIAS"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse BIAS TSV file with JSON rationale"""
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            standardized = pd.DataFrame({
                'Chr': df['chromosome'].astype(str).str.replace('chr', ''),
                'Pos': df['position'].astype(int),
                'Ref': df['refAllele'],
                'Alt': df['altAllele'],
                'Classification': df['acmgClassification'].apply(self.standardize_classification),
                'Tool': self.tool_name
            })
            
            criteria_list = []
            for idx, row in df.iterrows():
                criteria = self._extract_criteria_from_rationale(row['rationale'])
                criteria_list.append(','.join(sorted(criteria)) if criteria else '')
            
            standardized['ACMG_Criteria'] = criteria_list
            
            return self.create_variant_key(standardized)
            
        except Exception as e:
            print(f"Error parsing BIAS file {file_path}: {e}")
            return pd.DataFrame()
    
    def _extract_criteria_from_rationale(self, rationale: str) -> list:
        """Extract ACMG criteria from JSON rationale field"""
        if pd.isna(rationale) or rationale == '':
            return []
        
        try:
            rationale_dict = json.loads(rationale)
            criteria = []
            
            for category, codes in rationale_dict.items():
                if isinstance(codes, dict):
                    for code, value in codes.items():
                        if isinstance(value, list) and len(value) > 0:
                            if value[0] > 0:
                                criteria.append(code.upper())
            
            return criteria
        except:
            return self.extract_acmg_criteria(rationale)
