import pandas as pd
from base_parser import BaseParser

class AutoGVPParser(BaseParser):
    """Parser for AutoGVP preprocessed output files (*_acmg_criteria_parsed.tsv format)
    
    NOTE: AutoGVP ACMG Criteria Handling
    =====================================
    This parser works with PREPROCESSED AutoGVP files that have been parsed to extract
    specific ACMG criteria from the InterVar evidence column.
    
    The preprocessed files contain:
    - chr, pos, ref, alt: Variant coordinates
    - autogvp_call: AutoGVP classification (Pathogenic, Likely_pathogenic, etc.)
    - autogvp_criteria: Comma-separated list of specific ACMG criteria (e.g., "PM1, PM2, PP5")
    - autogvp_criteria_count: Number of criteria applied
    
    The original AutoGVP output uses simplified criteria groups (PM, BP, etc.) in binary
    columns, but the preprocessed files extract the actual specific criteria from the
    InterVar evidence field, making them suitable for detailed ACMG criteria analysis.
    
    Therefore, this parser:
    - DOES extract variant classifications (Pathogenic, Likely_Pathogenic, etc.)
    - DOES extract specific ACMG criteria from the autogvp_criteria column
    - AutoGVP can be used for both classification and criteria analysis
    """
    
    def __init__(self, tool_name: str = "AutoGVP"):
        super().__init__(tool_name)
    
    def parse(self, file_path: str) -> pd.DataFrame:
        """Parse AutoGVP parsed TSV file (*_acmg_criteria_parsed.tsv)"""
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            
            # Parse HGVSg to extract chr, pos, ref, alt
            # Format: chr1:g.9710459C>A
            parsed_variants = df['HGVSg'].apply(self._parse_hgvsg)
            
            standardized = pd.DataFrame({
                'Chr': [v[0] for v in parsed_variants],
                'Pos': [v[1] for v in parsed_variants],
                'Ref': [v[2] for v in parsed_variants],
                'Alt': [v[3] for v in parsed_variants],
                'Classification': df['autogvp_call'].apply(self.standardize_classification),
                'Tool': self.tool_name
            })
            
            # Extract ACMG criteria from active_criteria column
            # Format: "PM1, PM2, PP5" or "None" or empty
            criteria_list = []
            for idx, row in df.iterrows():
                criteria_str = str(row['active_criteria']) if 'active_criteria' in df.columns else ''
                if pd.isna(criteria_str) or criteria_str in ['None', 'nan', '', 'na']:
                    criteria_list.append('')
                else:
                    # Clean up the criteria string: remove spaces, split by comma
                    criteria = [c.strip() for c in criteria_str.split(',') if c.strip() and c.strip() not in ['None', 'na']]
                    criteria_list.append(','.join(sorted(criteria)) if criteria else '')
            
            standardized['ACMG_Criteria'] = criteria_list
            
            return self.create_variant_key(standardized)
            
        except Exception as e:
            print(f"Error parsing AutoGVP file {file_path}: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
    
    def _parse_hgvsg(self, hgvsg_str: str) -> tuple:
        """Parse HGVSg string to extract chr, pos, ref, alt
        
        Example: chr1:g.9710459C>A -> ('1', 9710459, 'C', 'A')
        """
        import re
        
        if pd.isna(hgvsg_str) or not hgvsg_str:
            return ('', 0, '', '')
        
        try:
            # Pattern: chr1:g.9710459C>A or chr1:g.9710459del or chr1:g.9710459_9710460insA
            match = re.match(r'chr(\w+):g\.(\d+)([A-Z]+)>([A-Z]+)', str(hgvsg_str))
            if match:
                chrom = match.group(1)
                pos = int(match.group(2))
                ref = match.group(3)
                alt = match.group(4)
                return (chrom, pos, ref, alt)
            
            # Handle deletions: chr1:g.9710459del
            match = re.match(r'chr(\w+):g\.(\d+)del([A-Z]*)', str(hgvsg_str))
            if match:
                chrom = match.group(1)
                pos = int(match.group(2))
                ref = match.group(3) if match.group(3) else 'N'
                alt = '-'
                return (chrom, pos, ref, alt)
            
            # If no pattern matches, return empty
            return ('', 0, '', '')
            
        except Exception:
            return ('', 0, '', '')
