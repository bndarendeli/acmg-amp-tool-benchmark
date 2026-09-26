from abc import ABC, abstractmethod
import pandas as pd
from typing import Dict, List, Tuple

class BaseParser(ABC):
    """
    Abstract base class for tool-specific parsers using Strategy Pattern.
    Each tool parser must implement the parse method.
    """
    
    def __init__(self, tool_name: str):
        self.tool_name = tool_name
    
    @abstractmethod
    def parse(self, file_path: str) -> pd.DataFrame:
        """
        Parse tool output file and return standardized DataFrame.
        
        Returns:
            DataFrame with columns: Chr, Pos, Ref, Alt, Classification, ACMG_Criteria
        """
        pass
    
    def create_variant_key(self, df: pd.DataFrame) -> pd.DataFrame:
        """Create unique variant identifier: Chr-Pos-Ref-Alt"""
        df['Variant_Key'] = (
            df['Chr'].astype(str) + '-' + 
            df['Pos'].astype(str) + '-' + 
            df['Ref'].astype(str) + '-' + 
            df['Alt'].astype(str)
        )
        return df
    
    def make_variant_key(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """Create variant key from individual components"""
        return f"{chrom}-{pos}-{ref}-{alt}"
    
    def standardize_classification(self, classification: str) -> str:
        """Standardize classification labels"""
        if pd.isna(classification) or classification in ['.', '', 'nan']:
            return 'Unknown'
        
        classification = str(classification).lower().replace('_', ' ').strip()
        
        if 'pathogenic' in classification and 'likely' in classification:
            return 'Likely_Pathogenic'
        elif 'pathogenic' in classification:
            return 'Pathogenic'
        elif 'benign' in classification and 'likely' in classification:
            return 'Likely_Benign'
        elif 'benign' in classification:
            return 'Benign'
        elif 'uncertain' in classification or 'vus' in classification:
            return 'Uncertain_Significance'
        else:
            return 'Unknown'
    
    def extract_acmg_criteria(self, criteria_string: str) -> List[str]:
        """Extract list of ACMG criteria from various formats"""
        if pd.isna(criteria_string) or criteria_string in ['.', '', 'nan']:
            return []
        
        criteria_string = str(criteria_string)
        criteria = []
        
        acmg_codes = [
            'PVS1', 'PS1', 'PS2', 'PS3', 'PS4',
            'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
            'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
            'BA1', 'BS1', 'BS2', 'BS3', 'BS4',
            'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
        ]
        
        for code in acmg_codes:
            if code in criteria_string.upper():
                criteria.append(code)
        
        return sorted(list(set(criteria)))
