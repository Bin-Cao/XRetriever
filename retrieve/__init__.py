"""
XRetriever - X-ray Diffraction Pattern Matching and Crystal Structure Retrieval
================================================================================

A powerful Python package for XRD pattern matching with state-of-the-art algorithms.

Author: Bin Cao
Affiliation: Hong Kong University of Science and Technology (Guangzhou)
Email: bcao686@connect.hkust-gz.edu.cn
GitHub: https://github.com/Bin-Cao/XRetriever

Components:
- xrd_reader: Read experimental XRD CSV/TXT files
- peak_detector: Advanced peak detection with preprocessing
- matcher: Pattern matching with dual scoring (Weighted + FOM)
- retriever: Main API for XRD matching and retrieval

Features:
- Robust peak detection with Savitzky-Golay filtering
- Hungarian algorithm-based optimal peak assignment
- Dual scoring metrics: Weighted Score + FOM (Figure of Merit)
- Element-based filtering for fast candidate screening
- Support for combined scoring mode
"""

import sys
import logging

# Package metadata
__version__ = '0.1.0'
__author__ = 'Bin Cao'
__email__ = 'bcao686@connect.hkust-gz.edu.cn'
__affiliation__ = 'Hong Kong University of Science and Technology (Guangzhou)'
__github__ = 'https://github.com/Bin-Cao/XRetriever'
__license__ = 'MIT'

# Import main classes
from .xrd_reader import XRDReader
from .peak_detector import PeakDetector
from .matcher import XRDMatcher
from .retriever import XRDRetriever

# Define public API
__all__ = [
    'XRDReader',
    'PeakDetector',
    'XRDMatcher',
    'XRDRetriever',
    '__version__',
    '__author__',
    '__email__',
    '__affiliation__',
    '__github__',
]

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Display package information when imported
def _print_package_info():
    """Print package information on import."""
    info_message = f"""
{'='*80}
XRetriever v{__version__} - X-ray Diffraction Pattern Matching
{'='*80}
Author:       {__author__}
Affiliation:  {__affiliation__}
Email:        {__email__}
GitHub:       {__github__}
License:      {__license__}
"""
    print(info_message, file=sys.stderr)

# Print info on import (can be disabled by setting environment variable)
import os
if os.environ.get('XRETRIEVER_QUIET', '').lower() not in ('1', 'true', 'yes'):
    _print_package_info()

