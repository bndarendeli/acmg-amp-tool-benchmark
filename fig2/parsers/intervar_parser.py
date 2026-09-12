import pandas as pd
from base_parser import BaseParser

class InterVarParser(BaseParser):
    """Parser for InterVar output files (.txt format)"""
    
    def __init__(self, tool_name: str = "InterVar"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse InterVar multianno.txt file"""
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            standardized = pd.DataFrame({
                'Chr': df['Chr'].astype(str).str.replace('chr', ''),
                'Pos': df['Start'].astype(int),
                'Ref': df['Ref'],
                'Alt': df['Alt'],
                'Classification': df['InterVar_automated'].apply(self.standardize_classification),
                'Tool': self.tool_name
            })
            
            acmg_columns = ['PVS1', 'PS1', 'PS2', 'PS3', 'PS4',
                           'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
                           'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
                           'BA1', 'BS1', 'BS2', 'BS3', 'BS4',
                           'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7']
            
            criteria_list = []
            for idx, row in df.iterrows():
                criteria = []
                for col in acmg_columns:
                    if col in df.columns and str(row[col]) == '1':
                        criteria.append(col)
                criteria_list.append(','.join(sorted(criteria)) if criteria else '')
            
            standardized['ACMG_Criteria'] = criteria_list
            
            return self.create_variant_key(standardized)
            
        except Exception as e:
            print(f"Error parsing InterVar file {file_path}: {e}")
            return pd.DataFrame()
