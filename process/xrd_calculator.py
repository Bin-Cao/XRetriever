"""
XRD Spectrum Calculator Module
===============================

This module calculates XRD spectra from crystal structures using pymatgen
and extracts the most intense peaks.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.diffraction.xrd import XRDCalculator as PymatgenXRDCalculator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class XRDCalculator:
    """
    Calculate XRD spectra from crystal structures.
    
    This class uses pymatgen's XRD calculator to compute powder diffraction patterns
    and extract the most intense peaks.
    """
    
    def __init__(self, 
                 wavelength: str = 'CuKa',
                 two_theta_range: Tuple[float, float] = (10, 90),
                 min_d_spacing: float = 0.5):
        """
        Initialize the XRD calculator.
        
        Args:
            wavelength: X-ray wavelength. Common options:
                       'CuKa' (1.54184 Å), 'MoKa' (0.71073 Å), 
                       'CrKa' (2.29100 Å), 'FeKa' (1.93604 Å),
                       or a float value in Angstroms
            two_theta_range: Range of 2θ angles to calculate (min, max) in degrees
            min_d_spacing: Minimum d-spacing to consider (Angstroms)
        """
        self.wavelength = wavelength
        self.two_theta_range = two_theta_range
        self.min_d_spacing = min_d_spacing

        # Initialize pymatgen XRD calculator
        self.xrd_calc = PymatgenXRDCalculator(wavelength=wavelength)

        logger.info(f"XRD Calculator initialized with wavelength={wavelength}, "
                   f"2θ range={two_theta_range}, min d-spacing={min_d_spacing} Å")
    
    def calculate_pattern(self, atoms: Atoms) -> Optional[Dict]:
        """
        Calculate XRD pattern for a crystal structure.
        
        Args:
            atoms: ASE Atoms object representing the crystal structure
            
        Returns:
            Dictionary containing XRD pattern data or None if calculation failed
        """
        try:
            # Convert ASE Atoms to pymatgen Structure
            adaptor = AseAtomsAdaptor()
            structure = adaptor.get_structure(atoms)
            
            # Calculate XRD pattern
            pattern = self.xrd_calc.get_pattern(structure, 
                                                two_theta_range=self.two_theta_range)
            
            # Extract data
            two_theta = np.array(pattern.x)  # 2θ angles
            intensities = np.array(pattern.y)  # Intensities
            hkls = pattern.hkls  # Miller indices (list of tuples)
            d_hkls = np.array(pattern.d_hkls)  # d-spacings

            # Filter by minimum d-spacing
            valid_indices = d_hkls >= self.min_d_spacing

            result = {
                'two_theta': two_theta[valid_indices],
                'intensities': intensities[valid_indices],
                'hkls': [hkls[i] for i in range(len(hkls)) if valid_indices[i]],
                'd_spacings': d_hkls[valid_indices]
            }

            return result
            
        except Exception as e:
            logger.error(f"Failed to calculate XRD pattern: {e}")
            return None
    
    def extract_top_peaks(self, 
                         pattern: Dict, 
                         n_peaks: int = 5,
                         normalize: bool = True) -> List[Dict]:
        """
        Extract the top N most intense peaks from an XRD pattern.
        
        Args:
            pattern: XRD pattern dictionary from calculate_pattern()
            n_peaks: Number of top peaks to extract
            normalize: Whether to normalize intensities to 100
            
        Returns:
            List of peak dictionaries, sorted by intensity (descending)
        """
        if pattern is None or len(pattern['two_theta']) == 0:
            logger.warning("Empty or invalid XRD pattern")
            return []
        
        two_theta = pattern['two_theta']
        intensities = pattern['intensities']
        hkls = pattern['hkls']
        d_spacings = pattern['d_spacings']
        
        # Normalize intensities if requested
        if normalize and len(intensities) > 0:
            max_intensity = np.max(intensities)
            if max_intensity > 0:
                intensities = 100 * intensities / max_intensity
        
        # Create list of peaks
        peaks = []
        for i in range(len(two_theta)):
            peak = {
                'two_theta': float(two_theta[i]),
                'intensity': float(intensities[i]),
                'hkl': hkls[i],
                'd_spacing': float(d_spacings[i])
            }
            peaks.append(peak)
        
        # Sort by intensity (descending)
        peaks.sort(key=lambda x: x['intensity'], reverse=True)
        
        # Return top N peaks
        top_peaks = peaks[:n_peaks]
        
        logger.debug(f"Extracted {len(top_peaks)} top peaks from {len(peaks)} total peaks")
        
        return top_peaks
    
    def calculate_and_extract_peaks(self, 
                                   atoms: Atoms, 
                                   n_peaks: int = 5,
                                   normalize: bool = True) -> Optional[List[Dict]]:
        """
        Calculate XRD pattern and extract top peaks in one step.
        
        Args:
            atoms: ASE Atoms object
            n_peaks: Number of top peaks to extract
            normalize: Whether to normalize intensities
            
        Returns:
            List of top peak dictionaries or None if calculation failed
        """
        pattern = self.calculate_pattern(atoms)
        
        if pattern is None:
            return None
        
        top_peaks = self.extract_top_peaks(pattern, n_peaks=n_peaks, normalize=normalize)
        
        return top_peaks
    
    def get_peak_positions_and_intensities(self, peaks: List[Dict]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract peak positions and intensities as numpy arrays.
        
        Args:
            peaks: List of peak dictionaries
            
        Returns:
            Tuple of (positions, intensities) as numpy arrays
        """
        if not peaks:
            return np.array([]), np.array([])
        
        positions = np.array([p['two_theta'] for p in peaks])
        intensities = np.array([p['intensity'] for p in peaks])
        
        return positions, intensities
    
    def format_peaks_for_storage(self, peaks: List[Dict]) -> Dict:
        """
        Format peaks for compact storage in database.
        
        Args:
            peaks: List of peak dictionaries
            
        Returns:
            Dictionary with compact peak representation
        """
        if not peaks:
            return {
                'positions': [],
                'intensities': [],
                'hkls': [],
                'd_spacings': []
            }
        
        return {
            'positions': [p['two_theta'] for p in peaks],
            'intensities': [p['intensity'] for p in peaks],
            'hkls': [p['hkl'] for p in peaks],
            'd_spacings': [p['d_spacing'] for p in peaks]
        }

