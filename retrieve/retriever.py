"""
XRD Retriever Module
====================

Main API for XRD pattern matching and crystal structure retrieval.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union
import numpy as np

from .xrd_reader import XRDReader
from .peak_detector import PeakDetector
from .matcher import XRDMatcher

# Set logger to WARNING level by default (suppress INFO messages)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class XRDRetriever:
    """
    Main interface for XRD pattern matching and retrieval.
    
    This class provides a high-level API that combines XRD file reading,
    peak detection, and pattern matching to find crystal structures that
    match experimental XRD data.
    """
    
    def __init__(self,
                 database_path: str,
                 position_tolerance: float = 0.2,
                 min_peak_height: float = 5.0,
                 min_peak_prominence: float = 3.0,
                 n_peaks: int = 5,
                 scoring_method: str = 'combined'):
        """
        Initialize the XRD retriever.

        Args:
            database_path: Path to the XRD database file (.pkl)
            position_tolerance: Maximum allowed difference in 2θ (degrees)
            min_peak_height: Minimum peak height for detection (% of max)
            min_peak_prominence: Minimum peak prominence (% of max)
            n_peaks: Number of top peaks to use for matching
            scoring_method: Scoring method ('weighted', 'fom', or 'combined')
        """
        self.database_path = database_path
        self.n_peaks = n_peaks
        self.scoring_method = scoring_method

        # Load database
        logger.info(f"Loading XRD database from {database_path}")
        self.database = self._load_database(database_path)

        # Initialize components
        self.reader = XRDReader()
        self.peak_detector = PeakDetector(
            min_peak_height=min_peak_height,
            min_peak_prominence=min_peak_prominence
        )
        self.matcher = XRDMatcher(
            position_tolerance=position_tolerance,
            scoring_method=scoring_method
        )
        
        logger.info("XRDRetriever initialized successfully")
    
    def _load_database(self, database_path: str) -> Dict:
        """Load the XRD database."""
        import pickle
        
        with open(database_path, 'rb') as f:
            database = pickle.load(f)
        
        n_entries = len(database['xrd_database'])
        logger.info(f"Loaded database with {n_entries} entries")
        
        return database
    
    def retrieve_from_file(self,
                          xrd_file: Union[str, Path],
                          elements: List[str],
                          top_n: int = 10,
                          auto_detect_format: bool = True) -> List[Dict]:
        """
        Retrieve matching crystal structures from an XRD CSV file.
        
        Args:
            xrd_file: Path to experimental XRD CSV file
            elements: List of chemical elements (e.g., ['Al', 'O'])
            top_n: Number of top matches to return
            auto_detect_format: Whether to auto-detect file format
            
        Returns:
            List of matching entries with scores
        """
        logger.info(f"Processing XRD file: {xrd_file}")
        logger.info(f"Required elements: {elements}")
        
        # Read XRD file
        if auto_detect_format:
            xrd_data = self.reader.read_auto(xrd_file)
        else:
            xrd_data = self.reader.read_csv(xrd_file)
        
        two_theta = xrd_data['two_theta']
        intensity = xrd_data['intensity']
        
        # Detect peaks
        logger.info("Detecting peaks...")
        peaks = self.peak_detector.get_top_peaks(
            two_theta,
            intensity,
            n_peaks=self.n_peaks,
            preprocess=True
        )
        
        if len(peaks) == 0:
            logger.warning("No peaks detected in experimental data")
            return []
        
        # Extract peak positions and intensities
        exp_positions, exp_intensities = self.peak_detector.extract_peak_positions_and_intensities(peaks)
        
        logger.info(f"Detected {len(peaks)} peaks:")
        for i, peak in enumerate(peaks[:5], 1):
            logger.info(f"  Peak {i}: 2θ={peak['two_theta']:.2f}°, I={peak['intensity']:.1f}")
        
        # Match against database
        logger.info("Matching against database...")
        results = self.matcher.match_pattern(
            exp_positions,
            exp_intensities,
            self.database,
            elements=elements,
            top_n=top_n
        )
        
        return results
    
    def retrieve_from_peaks(self,
                           peak_positions: List[float],
                           peak_intensities: List[float],
                           elements: List[str],
                           top_n: int = 10) -> List[Dict]:
        """
        Retrieve matching crystal structures from peak data.
        
        Args:
            peak_positions: List of peak positions (2θ angles)
            peak_intensities: List of peak intensities
            elements: List of chemical elements
            top_n: Number of top matches to return
            
        Returns:
            List of matching entries with scores
        """
        logger.info(f"Matching {len(peak_positions)} peaks")
        logger.info(f"Required elements: {elements}")
        
        exp_positions = np.array(peak_positions)
        exp_intensities = np.array(peak_intensities)
        
        # Match against database
        results = self.matcher.match_pattern(
            exp_positions,
            exp_intensities,
            self.database,
            elements=elements,
            top_n=top_n
        )
        
        return results
    
    def get_entry_details(self, entry_id: int) -> Optional[Dict]:
        """
        Get detailed information about a database entry.
        
        Args:
            entry_id: Database entry ID
            
        Returns:
            Entry details or None if not found
        """
        xrd_db = self.database['xrd_database']
        
        if entry_id in xrd_db:
            return xrd_db[entry_id]
        else:
            logger.warning(f"Entry {entry_id} not found in database")
            return None
    
    def search_by_formula(self, formula: str) -> List[Dict]:
        """
        Search database by chemical formula.
        
        Args:
            formula: Chemical formula (e.g., 'Al2O3')
            
        Returns:
            List of matching entries
        """
        xrd_db = self.database['xrd_database']
        matches = []
        
        for entry_id, entry in xrd_db.items():
            if entry['formula'] == formula:
                matches.append(entry)
        
        logger.info(f"Found {len(matches)} entries with formula {formula}")
        
        return matches
    
    def search_by_elements(self, elements: List[str]) -> List[Dict]:
        """
        Search database by chemical elements.
        
        Args:
            elements: List of element symbols
            
        Returns:
            List of matching entries
        """
        candidate_ids = self.matcher.filter_by_elements(self.database, elements)
        
        xrd_db = self.database['xrd_database']
        matches = [xrd_db[entry_id] for entry_id in candidate_ids]
        
        return matches
    
    def get_database_statistics(self) -> Dict:
        """
        Get statistics about the loaded database.
        
        Returns:
            Dictionary with database statistics
        """
        xrd_db = self.database['xrd_database']
        metadata = self.database['metadata']
        
        # Collect element statistics
        all_elements = set()
        formulas = []
        
        for entry in xrd_db.values():
            all_elements.update(entry['elements'])
            formulas.append(entry['formula'])
        
        stats = {
            'total_entries': len(xrd_db),
            'unique_elements': sorted(list(all_elements)),
            'n_unique_elements': len(all_elements),
            'unique_formulas': len(set(formulas)),
            'wavelength': metadata['wavelength'],
            'n_peaks_per_entry': metadata['n_peaks'],
            'two_theta_range': metadata['two_theta_range']
        }
        
        return stats
    
    def print_results(self, results: List[Dict], max_results: int = 10):
        """
        Print matching results in a formatted way.

        Args:
            results: List of match results
            max_results: Maximum number of results to print
        """
        if len(results) == 0:
            print("No matches found.")
            return

        print(f"\nTop {min(len(results), max_results)} Matches:")
        print("=" * 80)

        for i, result in enumerate(results[:max_results], 1):
            print(f"\n{i}. Score: {result['score']:.2f}", end="")

            # Display FOM if available (combined scoring mode)
            if 'fom' in result:
                print(f" | FOM: {result['fom']:.2f}", end="")

            # Display matched peaks if available
            if 'n_matched_peaks' in result:
                print(f" | Matched: {result['n_matched_peaks']} peaks", end="")

            print()  # New line
            print(f"   MPID: {result['mpid']}")
            print(f"   Formula: {result['formula']}")
            print(f"   Elements: {', '.join(result['elements'])}")
            print(f"   Space Group: {result['spacegroup']}")

        print("\n" + "=" * 80)

