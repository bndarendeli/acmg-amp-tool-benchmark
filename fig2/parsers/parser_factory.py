from intervar_parser import InterVarParser
from bias_parser import BIASParser
from genebe_parser import GenebeParser
from vcf_parser import VCFParser
from charger_parser import CharGerParser
from cancer_sigvar_parser import CancerSIGVARParser
from diabloacmg_parser import DiabloACMGParser
from tapes_parser import TAPESParser
from viphl_parser import VIPHLParser
from cpsr_parser import CPSRParser
from exomiser_parser import ExomiserParser
from franklin_parser import FranklinParser
from autogvp_parser import AutoGVPParser
import os

class ParserFactory:
    """Factory class to create appropriate parser based on tool name and file format"""
    
    @staticmethod
    def get_parser(tool_name: str, file_path: str):
        """
        Return appropriate parser based on tool name and file extension.
        
        Args:
            tool_name: Name of the ACMG tool
            file_path: Path to the file to be parsed
            
        Returns:
            Parser instance
        """
        tool_lower = tool_name.lower()
        file_lower = file_path.lower()
        
        # Tool-specific parsers
        if 'autogvp' in tool_lower:
            return AutoGVPParser(tool_name)
        elif 'intervar' in tool_lower:
            return InterVarParser(tool_name)
        elif 'bias' in tool_lower:
            return BIASParser(tool_name)
        elif 'genebe' in tool_lower or 'gene.be' in tool_lower:
            return GenebeParser(tool_name)
        elif 'charger' in tool_lower or 'charger_local' in tool_lower or 'charger_online' in tool_lower:
            return CharGerParser(tool_name)
        elif 'cancer' in tool_lower and 'sigvar' in tool_lower:
            return CancerSIGVARParser(tool_name)
        elif 'diablo' in tool_lower or 'diabloacmg' in tool_lower:
            return DiabloACMGParser(tool_name)
        elif 'tapes' in tool_lower:
            return TAPESParser(tool_name)
        elif 'vip' in tool_lower and 'hl' in tool_lower:
            return VIPHLParser(tool_name)
        elif 'cpsr' in tool_lower:
            return CPSRParser(tool_name)
        elif 'exomiser' in tool_lower:
            return ExomiserParser(tool_name)
        elif 'franklin' in tool_lower:
            return FranklinParser(tool_name)
        # File extension-based fallback
        elif file_path.endswith('.vcf'):
            return VCFParser(tool_name)
        elif file_path.endswith(('.tsv', '.txt')):
            if 'intervar' in file_path.lower():
                return InterVarParser(tool_name)
            elif 'bias' in file_path.lower():
                return BIASParser(tool_name)
            else:
                return VCFParser(tool_name)
        else:
            return VCFParser(tool_name)
