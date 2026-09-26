#!/usr/bin/env python3
"""
Utility functions for parsing ACMG criteria.

Ground truth criteria are often combined with underscores (e.g., PVS1_PS2_PM6).
These need to be split into individual criteria for proper analysis.
"""

import pandas as pd
import re


def parse_criteria_string(criteria_str):
    """
    Parse criteria string into a set of individual criteria.
    
    Handles combined criteria with multiple separators:
    - "PVS1_PS2_PM6" → {"PVS1", "PS2", "PM6"}
    - "PM2, PP4_PM, PVS1_PS1" → {"PM2", "PP4", "PVS1", "PS1"} (PM excluded)
    - "PVS1&PM2&PP5" → {"PVS1", "PM2", "PP5"} (Genebe format with &)
    
    Excludes incomplete criteria: PS, PM, PP (without numbers)
    
    Args:
        criteria_str: String containing ACMG criteria
        
    Returns:
        Set of individual ACMG criteria codes
    """
    if pd.isna(criteria_str) or criteria_str == '' or criteria_str == 'Unknown':
        return set()
    
    # Split by comma first
    criteria_parts = [c.strip() for c in str(criteria_str).split(',')]
    
    individual_criteria = set()
    
    for part in criteria_parts:
        if not part or part in ['None', 'nan', 'na', 'Unknown']:
            continue
        
        # Split combined criteria by underscore AND ampersand
        # e.g., "PVS1_PS2_PM6" → ["PVS1", "PS2", "PM6"]
        # e.g., "PVS1&PM2&PP5" → ["PVS1", "PM2", "PP5"]
        # Replace & with _ first for uniform splitting
        part = part.replace('&', '_')
        sub_criteria = part.split('_')
        
        for criterion in sub_criteria:
            criterion = criterion.strip()
            
            # Validate it looks like an ACMG criterion (starts with P, B, or M)
            # Must have a number after the prefix (e.g., PVS1, PM2, PS1, not just PS, PM, PP)
            # Pattern: P/B/M + one or more letters + one or more digits
            if criterion and re.match(r'^[PBM][PVSMABP]+\d+', criterion):
                individual_criteria.add(criterion)
    
    return individual_criteria


def parse_criteria_to_list(criteria_str):
    """Parse criteria string and return as sorted list"""
    return sorted(list(parse_criteria_string(criteria_str)))


def parse_criteria_to_string(criteria_str, separator=','):
    """Parse criteria string and return as comma-separated string"""
    criteria_list = parse_criteria_to_list(criteria_str)
    return separator.join(criteria_list) if criteria_list else ''


# Test the parser
if __name__ == "__main__":
    test_cases = [
        "PM2, PP4_PM, PVS1_PS1",
        "PVS1_PS2_PM6",
        "PVS1&PM2&PP5",  # Genebe format
        "PM1&PM2&PP2&PP3",  # Genebe format
        "PM2",
        "PP4_PM",  # Should exclude PM
        "PS, PM, PP",  # Should exclude all (no numbers)
        "",
        "None"
    ]
    
    print("Testing criteria parser:")
    print("="*60)
    for test in test_cases:
        result = parse_criteria_string(test)
        print(f"Input:  '{test}'")
        print(f"Output: {result}")
        print()
