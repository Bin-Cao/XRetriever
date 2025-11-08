"""
Peak Detection Module
=====================

This module implements robust peak detection algorithms for XRD spectra.
"""

import logging
import numpy as np
from typing import Dict, List, Optional, Tuple
from scipy.signal import find_peaks, peak_widths, savgol_filter
from scipy.ndimage import gaussian_filter1d

# Set logger to WARNING level by default (suppress INFO messages)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class PeakDetector:
    """
    Detect peaks in XRD spectra using robust algorithms.
    
    This class implements state-of-the-art peak detection methods that handle
    noise, baseline variations, and overlapping peaks in experimental XRD data.
    """
    
    def __init__(self,
                 min_peak_height: float = 5.0,
                 min_peak_prominence: float = 3.0,
                 min_peak_distance: float = 0.1,
                 smooth_window: int = 5):
        """
        Initialize the peak detector.
        
        Args:
            min_peak_height: Minimum peak height (% of max intensity)
            min_peak_prominence: Minimum peak prominence (% of max intensity)
            min_peak_distance: Minimum distance between peaks (degrees in 2θ)
            smooth_window: Window size for smoothing (must be odd)
        """
        self.min_peak_height = min_peak_height
        self.min_peak_prominence = min_peak_prominence
        self.min_peak_distance = min_peak_distance
        self.smooth_window = smooth_window if smooth_window % 2 == 1 else smooth_window + 1
        
        logger.info(f"PeakDetector initialized with height={min_peak_height}%, "
                   f"prominence={min_peak_prominence}%, distance={min_peak_distance}°")
    
    def preprocess_spectrum(self,
                           two_theta: np.ndarray,
                           intensity: np.ndarray,
                           normalize: bool = True,
                           remove_baseline: bool = True,
                           smooth: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Preprocess XRD spectrum before peak detection.
        
        Args:
            two_theta: 2θ angles
            intensity: Intensity values
            normalize: Whether to normalize to 100
            remove_baseline: Whether to remove baseline
            smooth: Whether to apply smoothing
            
        Returns:
            Tuple of (two_theta, processed_intensity)
        """
        processed_intensity = intensity.copy()
        
        # Normalize to 100
        if normalize:
            max_val = np.max(processed_intensity)
            if max_val > 0:
                processed_intensity = 100 * processed_intensity / max_val
        
        # Remove baseline
        if remove_baseline:
            processed_intensity = self._remove_baseline(processed_intensity)
        
        # Smooth the data
        if smooth and len(processed_intensity) > self.smooth_window:
            try:
                # Use Savitzky-Golay filter for smoothing
                processed_intensity = savgol_filter(
                    processed_intensity,
                    window_length=self.smooth_window,
                    polyorder=2,
                    mode='nearest'
                )
            except Exception as e:
                logger.warning(f"Smoothing failed, using Gaussian filter: {e}")
                processed_intensity = gaussian_filter1d(processed_intensity, sigma=1.0)
        
        return two_theta, processed_intensity
    
    def _remove_baseline(self, intensity: np.ndarray, percentile: float = 5) -> np.ndarray:
        """
        Remove baseline from intensity data.
        
        Args:
            intensity: Intensity values
            percentile: Percentile for baseline estimation
            
        Returns:
            Baseline-corrected intensity
        """
        # Simple baseline removal using rolling minimum
        from scipy.ndimage import minimum_filter1d
        
        window_size = max(len(intensity) // 20, 10)
        baseline = minimum_filter1d(intensity, size=window_size, mode='nearest')
        
        # Subtract baseline and ensure non-negative
        corrected = intensity - baseline
        corrected[corrected < 0] = 0
        
        return corrected
    
    def detect_peaks(self,
                    two_theta: np.ndarray,
                    intensity: np.ndarray,
                    preprocess: bool = True) -> List[Dict]:
        """
        Detect peaks in XRD spectrum.
        
        Args:
            two_theta: 2θ angles
            intensity: Intensity values
            preprocess: Whether to preprocess the data
            
        Returns:
            List of peak dictionaries with position, intensity, and properties
        """
        # Preprocess if requested
        if preprocess:
            two_theta, intensity = self.preprocess_spectrum(two_theta, intensity)
        
        # Calculate minimum distance in indices
        if len(two_theta) > 1:
            avg_step = np.mean(np.diff(two_theta))
            min_distance_indices = max(1, int(self.min_peak_distance / avg_step))
        else:
            min_distance_indices = 1
        
        # Calculate absolute thresholds
        max_intensity = np.max(intensity)
        height_threshold = self.min_peak_height * max_intensity / 100
        prominence_threshold = self.min_peak_prominence * max_intensity / 100
        
        # Find peaks using scipy
        peak_indices, properties = find_peaks(
            intensity,
            height=height_threshold,
            prominence=prominence_threshold,
            distance=min_distance_indices
        )
        
        if len(peak_indices) == 0:
            logger.warning("No peaks detected with current parameters")
            return []
        
        # Calculate peak widths
        try:
            widths, width_heights, left_ips, right_ips = peak_widths(
                intensity,
                peak_indices,
                rel_height=0.5
            )
        except Exception as e:
            logger.warning(f"Could not calculate peak widths: {e}")
            widths = np.ones(len(peak_indices))
        
        # Create peak list
        peaks = []
        for i, idx in enumerate(peak_indices):
            peak = {
                'index': int(idx),
                'two_theta': float(two_theta[idx]),
                'intensity': float(intensity[idx]),
                'height': float(properties['peak_heights'][i]),
                'prominence': float(properties['prominences'][i]),
                'width': float(widths[i]) * avg_step if len(two_theta) > 1 else 0.0
            }
            peaks.append(peak)
        
        # Sort by intensity (descending)
        peaks.sort(key=lambda x: x['intensity'], reverse=True)
        
        logger.info(f"Detected {len(peaks)} peaks")
        
        return peaks
    
    def get_top_peaks(self,
                     two_theta: np.ndarray,
                     intensity: np.ndarray,
                     n_peaks: int = 5,
                     preprocess: bool = True,
                     normalize: bool = True) -> List[Dict]:
        """
        Detect and return the top N most intense peaks.

        Args:
            two_theta: 2θ angles
            intensity: Intensity values
            n_peaks: Number of top peaks to return (default: 5)
            preprocess: Whether to preprocess the data
            normalize: Whether to normalize strongest peak to 100

        Returns:
            List of top N peak dictionaries, sorted by intensity (descending)
        """
        all_peaks = self.detect_peaks(two_theta, intensity, preprocess=preprocess)

        if len(all_peaks) == 0:
            logger.warning("No peaks detected")
            return []

        # Sort by intensity (descending) to ensure we get the strongest peaks
        all_peaks.sort(key=lambda x: x['intensity'], reverse=True)

        # Return top N peaks
        top_peaks = all_peaks[:n_peaks]

        # Normalize intensities if requested
        if normalize and len(top_peaks) > 0:
            max_intensity = max(p['intensity'] for p in top_peaks)
            if max_intensity > 0:
                for peak in top_peaks:
                    peak['intensity'] = (peak['intensity'] / max_intensity) * 100.0
                logger.info(f"Normalized peak intensities to max=100")

        logger.info(f"Returning top {len(top_peaks)} peaks out of {len(all_peaks)} detected")

        # Warn if fewer peaks than requested
        if len(top_peaks) < n_peaks:
            logger.warning(f"Only {len(top_peaks)} peaks detected, requested {n_peaks}")

        return top_peaks
    
    def extract_peak_positions_and_intensities(self, 
                                              peaks: List[Dict]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract peak positions and intensities as arrays.
        
        Args:
            peaks: List of peak dictionaries
            
        Returns:
            Tuple of (positions, intensities)
        """
        if not peaks:
            return np.array([]), np.array([])
        
        positions = np.array([p['two_theta'] for p in peaks])
        intensities = np.array([p['intensity'] for p in peaks])
        
        return positions, intensities
    
    def refine_peak_position(self,
                           two_theta: np.ndarray,
                           intensity: np.ndarray,
                           peak_index: int,
                           window: int = 5) -> float:
        """
        Refine peak position using parabolic interpolation.
        
        Args:
            two_theta: 2θ angles
            intensity: Intensity values
            peak_index: Index of the peak
            window: Window size around peak for fitting
            
        Returns:
            Refined peak position
        """
        # Get window around peak
        start = max(0, peak_index - window)
        end = min(len(intensity), peak_index + window + 1)
        
        if end - start < 3:
            return two_theta[peak_index]
        
        # Extract window data
        x = two_theta[start:end]
        y = intensity[start:end]
        
        # Find maximum in window
        max_idx = np.argmax(y)
        
        # Parabolic interpolation
        if 0 < max_idx < len(y) - 1:
            # Use three points around maximum
            y0, y1, y2 = y[max_idx - 1], y[max_idx], y[max_idx + 1]
            x0, x1, x2 = x[max_idx - 1], x[max_idx], x[max_idx + 1]
            
            # Parabolic fit
            denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
            if abs(denom) > 1e-10:
                A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
                B = (x2**2 * (y0 - y1) + x1**2 * (y2 - y0) + x0**2 * (y1 - y2)) / denom
                
                if abs(A) > 1e-10:
                    refined_pos = -B / (2 * A)
                    
                    # Check if refined position is reasonable
                    if x0 <= refined_pos <= x2:
                        return refined_pos
        
        return two_theta[peak_index]

