"""
XRD File Reader Module
======================

This module handles reading experimental XRD data from CSV files.
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

# Set logger to WARNING level by default (suppress INFO messages)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class XRDReader:
    """
    Read and parse experimental XRD data from CSV files.
    
    This class handles various CSV formats and provides methods to extract
    2θ angles and intensity values from experimental XRD measurements.
    """
    
    def __init__(self):
        """Initialize the XRD reader."""
        logger.info("XRDReader initialized")
    
    def read_csv(self, 
                 file_path: Union[str, Path],
                 two_theta_col: int = 0,
                 intensity_col: int = 1,
                 delimiter: str = ',',
                 skip_rows: int = 0,
                 header: Optional[int] = None) -> Dict[str, np.ndarray]:
        """
        Read XRD data from a CSV file.
        
        Args:
            file_path: Path to the CSV file
            two_theta_col: Column index for 2θ angles (0-based)
            intensity_col: Column index for intensity values (0-based)
            delimiter: Column delimiter (default: ',')
            skip_rows: Number of rows to skip at the beginning
            header: Row number to use as column names (None for no header)
            
        Returns:
            Dictionary with 'two_theta' and 'intensity' arrays
        """
        try:
            file_path = Path(file_path)
            
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")
            
            # Try to read the file
            try:
                # First attempt: use pandas
                df = pd.read_csv(
                    file_path,
                    delimiter=delimiter,
                    skiprows=skip_rows,
                    header=header
                )
                
                # Extract columns
                columns = df.columns.tolist()
                
                if isinstance(columns[0], str):
                    # Has header, use column indices
                    two_theta = df.iloc[:, two_theta_col].values
                    intensity = df.iloc[:, intensity_col].values
                else:
                    # No header
                    two_theta = df.iloc[:, two_theta_col].values
                    intensity = df.iloc[:, intensity_col].values
                    
            except Exception as e:
                logger.warning(f"Pandas read failed, trying numpy: {e}")
                # Fallback: use numpy
                data = np.loadtxt(
                    file_path,
                    delimiter=delimiter,
                    skiprows=skip_rows
                )
                
                if data.ndim == 1:
                    raise ValueError("Data must have at least 2 columns")
                
                two_theta = data[:, two_theta_col]
                intensity = data[:, intensity_col]
            
            # Convert to float and remove NaN values
            two_theta = np.array(two_theta, dtype=float)
            intensity = np.array(intensity, dtype=float)
            
            # Remove NaN and infinite values
            valid_mask = np.isfinite(two_theta) & np.isfinite(intensity)
            two_theta = two_theta[valid_mask]
            intensity = intensity[valid_mask]
            
            # Remove negative intensities
            positive_mask = intensity >= 0
            two_theta = two_theta[positive_mask]
            intensity = intensity[positive_mask]
            
            # Sort by 2θ
            sort_indices = np.argsort(two_theta)
            two_theta = two_theta[sort_indices]
            intensity = intensity[sort_indices]
            
            logger.info(f"Successfully read {len(two_theta)} data points from {file_path}")
            logger.info(f"2θ range: {two_theta.min():.2f}° - {two_theta.max():.2f}°")
            logger.info(f"Intensity range: {intensity.min():.2f} - {intensity.max():.2f}")
            
            return {
                'two_theta': two_theta,
                'intensity': intensity,
                'file_path': str(file_path)
            }
            
        except Exception as e:
            logger.error(f"Failed to read XRD file {file_path}: {e}")
            raise
    
    def read_auto(self, file_path: Union[str, Path]) -> Dict[str, np.ndarray]:
        """
        Automatically detect format and read XRD data.
        
        This method tries different common formats to read the file.
        
        Args:
            file_path: Path to the XRD data file
            
        Returns:
            Dictionary with 'two_theta' and 'intensity' arrays
        """
        file_path = Path(file_path)
        
        # Try different delimiters and formats
        delimiters = [',', '\t', ' ', ';']
        
        for delimiter in delimiters:
            try:
                data = self.read_csv(
                    file_path,
                    delimiter=delimiter,
                    two_theta_col=0,
                    intensity_col=1
                )
                
                # Check if data looks reasonable
                if len(data['two_theta']) > 10:
                    logger.info(f"Successfully auto-detected format with delimiter '{delimiter}'")
                    return data
                    
            except Exception:
                continue
        
        # If all attempts failed, raise error
        raise ValueError(f"Could not automatically read file {file_path}. "
                        "Please specify format parameters manually.")
    
    def normalize_intensity(self, 
                           intensity: np.ndarray,
                           method: str = 'max') -> np.ndarray:
        """
        Normalize intensity values.
        
        Args:
            intensity: Array of intensity values
            method: Normalization method ('max', 'sum', 'minmax')
            
        Returns:
            Normalized intensity array
        """
        if method == 'max':
            # Normalize to maximum = 100
            max_val = np.max(intensity)
            if max_val > 0:
                return 100 * intensity / max_val
            return intensity
            
        elif method == 'sum':
            # Normalize to sum = 1
            sum_val = np.sum(intensity)
            if sum_val > 0:
                return intensity / sum_val
            return intensity
            
        elif method == 'minmax':
            # Normalize to range [0, 1]
            min_val = np.min(intensity)
            max_val = np.max(intensity)
            if max_val > min_val:
                return (intensity - min_val) / (max_val - min_val)
            return intensity
            
        else:
            raise ValueError(f"Unknown normalization method: {method}")
    
    def smooth_data(self,
                   two_theta: np.ndarray,
                   intensity: np.ndarray,
                   window_size: int = 5) -> Tuple[np.ndarray, np.ndarray]:
        """
        Apply smoothing to XRD data.
        
        Args:
            two_theta: 2θ angles
            intensity: Intensity values
            window_size: Size of smoothing window
            
        Returns:
            Tuple of (smoothed_two_theta, smoothed_intensity)
        """
        from scipy.ndimage import uniform_filter1d
        
        smoothed_intensity = uniform_filter1d(intensity, size=window_size, mode='nearest')
        
        return two_theta, smoothed_intensity
    
    def resample_data(self,
                     two_theta: np.ndarray,
                     intensity: np.ndarray,
                     step: float = 0.02) -> Tuple[np.ndarray, np.ndarray]:
        """
        Resample XRD data to uniform 2θ spacing.
        
        Args:
            two_theta: Original 2θ angles
            intensity: Original intensity values
            step: Desired step size in 2θ (degrees)
            
        Returns:
            Tuple of (resampled_two_theta, resampled_intensity)
        """
        from scipy.interpolate import interp1d
        
        # Create uniform grid
        two_theta_min = np.min(two_theta)
        two_theta_max = np.max(two_theta)
        new_two_theta = np.arange(two_theta_min, two_theta_max, step)
        
        # Interpolate
        interpolator = interp1d(two_theta, intensity, kind='linear', 
                               bounds_error=False, fill_value=0)
        new_intensity = interpolator(new_two_theta)
        
        return new_two_theta, new_intensity

