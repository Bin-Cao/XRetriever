"""
XRD Pattern Matching Module
============================

This module implements state-of-the-art algorithms for matching experimental
XRD patterns to database entries.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional, Set, Union
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

# Set logger to WARNING level by default (suppress INFO messages)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class XRDMatcher:
    """
    Match experimental XRD patterns to database entries.
    
    This class implements robust matching algorithms that account for:
    - Peak position shifts (instrument calibration, sample displacement)
    - Intensity variations (preferred orientation, sample preparation)
    - Missing or extra peaks (impurities, detection limits)
    """
    
    def __init__(self,
                 position_tolerance: float = 0.2,
                 intensity_weight: float = 0.3,
                 position_weight: float = 0.7,
                 min_matched_peaks: int = 2,
                 scoring_method: str = 'weighted'):
        """
        Initialize the XRD matcher.

        Args:
            position_tolerance: Maximum allowed difference in 2θ (degrees)
            intensity_weight: Weight for intensity similarity (0-1)
            position_weight: Weight for position similarity (0-1)
            min_matched_peaks: Minimum number of peaks that must match
            scoring_method: Scoring method to use:
                - 'weighted': 70% position + 30% intensity similarity (default)
                - 'fom': Figure of Merit calculation
                - 'combined': Return both weighted and FOM scores
        """
        self.position_tolerance = position_tolerance
        self.intensity_weight = intensity_weight
        self.position_weight = position_weight
        self.min_matched_peaks = min_matched_peaks
        self.scoring_method = scoring_method

        # Validate scoring method
        valid_methods = ['weighted', 'fom', 'combined']
        if scoring_method not in valid_methods:
            raise ValueError(f"scoring_method must be one of {valid_methods}, got '{scoring_method}'")

        # Normalize weights
        total_weight = intensity_weight + position_weight
        self.intensity_weight = intensity_weight / total_weight
        self.position_weight = position_weight / total_weight

        logger.info(f"XRDMatcher initialized with tolerance={position_tolerance}°, "
                   f"weights=(pos:{self.position_weight:.2f}, int:{self.intensity_weight:.2f}), "
                   f"scoring_method='{scoring_method}'")
    
    def filter_by_elements(self,
                          database: Dict,
                          query_elements: List[str]) -> List[int]:
        """
        Filter database entries by chemical elements.
        
        The database entry must contain AT LEAST the query elements
        (but can contain additional elements).
        
        Args:
            database: XRD database dictionary
            query_elements: List of required elements
            
        Returns:
            List of matching entry IDs
        """
        query_set = set(query_elements)
        matching_ids = []
        
        xrd_db = database['xrd_database']
        
        for entry_id, entry in xrd_db.items():
            entry_elements = set(entry['elements'])
            
            # Check if all query elements are present
            if query_set.issubset(entry_elements):
                matching_ids.append(entry_id)
        
        logger.info(f"Found {len(matching_ids)} entries containing elements {query_elements}")
        
        return matching_ids
    
    def calculate_peak_match_score(self,
                                   exp_positions: np.ndarray,
                                   exp_intensities: np.ndarray,
                                   db_positions: np.ndarray,
                                   db_intensities: np.ndarray) -> float:
        """
        Calculate similarity score between experimental and database peaks.
        
        Uses Hungarian algorithm for optimal peak assignment.
        
        Args:
            exp_positions: Experimental peak positions (2θ)
            exp_intensities: Experimental peak intensities
            db_positions: Database peak positions (2θ)
            db_intensities: Database peak intensities
            
        Returns:
            Similarity score (0-100, higher is better)
        """
        if len(exp_positions) == 0 or len(db_positions) == 0:
            return 0.0
        
        # Normalize intensities
        exp_int_norm = exp_intensities / np.max(exp_intensities) if np.max(exp_intensities) > 0 else exp_intensities
        db_int_norm = db_intensities / np.max(db_intensities) if np.max(db_intensities) > 0 else db_intensities
        
        # Calculate position distance matrix
        pos_dist = cdist(
            exp_positions.reshape(-1, 1),
            db_positions.reshape(-1, 1),
            metric='euclidean'
        )
        
        # Calculate intensity difference matrix
        int_diff = np.abs(
            exp_int_norm.reshape(-1, 1) - db_int_norm.reshape(1, -1)
        )
        
        # Combined cost matrix (lower is better)
        # Position cost: normalized by tolerance
        pos_cost = pos_dist / self.position_tolerance
        # Intensity cost: already in [0, 1] range
        int_cost = int_diff
        
        # Weighted combination
        cost_matrix = (self.position_weight * pos_cost + 
                      self.intensity_weight * int_cost)
        
        # Apply position tolerance constraint
        # Set cost to infinity for peaks outside tolerance
        cost_matrix[pos_dist > self.position_tolerance] = np.inf

        # Check if there are any valid matches at all
        if np.all(np.isinf(cost_matrix)):
            # No peaks within tolerance
            return 0.0

        # Solve assignment problem using Hungarian algorithm
        try:
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
        except ValueError:
            # Assignment problem is infeasible (shouldn't happen after the check above)
            return 0.0

        # Calculate matched peaks
        valid_matches = cost_matrix[row_ind, col_ind] < np.inf
        n_matched = np.sum(valid_matches)

        if n_matched < self.min_matched_peaks:
            return 0.0
        
        # Calculate match quality
        matched_costs = cost_matrix[row_ind, col_ind][valid_matches]
        
        if len(matched_costs) == 0:
            return 0.0
        
        # Convert cost to similarity score
        # Lower cost = higher similarity
        avg_cost = np.mean(matched_costs)
        
        # Score based on:
        # 1. Number of matched peaks (coverage)
        # 2. Quality of matches (average cost)
        coverage_score = n_matched / max(len(exp_positions), len(db_positions))
        quality_score = np.exp(-avg_cost)  # Exponential decay of cost
        
        # Combined score
        score = 100 * coverage_score * quality_score

        return float(score)

    def calculate_fom_score(self,
                           exp_positions: np.ndarray,
                           exp_intensities: np.ndarray,
                           db_positions: np.ndarray,
                           db_intensities: np.ndarray) -> Tuple[float, int]:
        """
        Calculate Figure of Merit (FOM) score.

        FOM = (Σ I_matched / Σ I_total) × 100

        Args:
            exp_positions: Experimental peak positions (2θ)
            exp_intensities: Experimental peak intensities
            db_positions: Database peak positions (2θ)
            db_intensities: Database peak intensities

        Returns:
            Tuple of (FOM score, number of matched peaks)
        """
        if len(exp_positions) == 0 or len(db_positions) == 0:
            return 0.0, 0

        # Total intensity in database
        total_intensity = np.sum(db_intensities)

        if total_intensity == 0:
            return 0.0, 0

        # Find matched peaks
        matched_intensity = 0.0
        n_matched = 0

        for db_pos, db_int in zip(db_positions, db_intensities):
            # Check if any experimental peak is within tolerance
            pos_diff = np.abs(exp_positions - db_pos)
            min_diff = np.min(pos_diff)

            if min_diff <= self.position_tolerance:
                # Peak is matched
                matched_intensity += db_int
                n_matched += 1

        # Calculate FOM
        fom = (matched_intensity / total_intensity) * 100.0

        return float(fom), n_matched

    def match_single_entry(self,
                          exp_positions: np.ndarray,
                          exp_intensities: np.ndarray,
                          entry: Dict) -> Union[float, Dict[str, float]]:
        """
        Match experimental peaks to a single database entry.

        Args:
            exp_positions: Experimental peak positions
            exp_intensities: Experimental peak intensities
            entry: Database entry dictionary

        Returns:
            Match score (0-100) or dict with multiple scores if scoring_method='combined'
        """
        db_positions = np.array(entry['peaks']['positions'])
        db_intensities = np.array(entry['peaks']['intensities'])

        if self.scoring_method == 'weighted':
            score = self.calculate_peak_match_score(
                exp_positions,
                exp_intensities,
                db_positions,
                db_intensities
            )
            return score

        elif self.scoring_method == 'fom':
            fom, n_matched = self.calculate_fom_score(
                exp_positions,
                exp_intensities,
                db_positions,
                db_intensities
            )
            return fom

        else:  # 'combined'
            weighted_score = self.calculate_peak_match_score(
                exp_positions,
                exp_intensities,
                db_positions,
                db_intensities
            )
            fom, n_matched = self.calculate_fom_score(
                exp_positions,
                exp_intensities,
                db_positions,
                db_intensities
            )
            return {
                'weighted_score': weighted_score,
                'fom': fom,
                'n_matched_peaks': n_matched
            }
    
    def match_pattern(self,
                     exp_positions: np.ndarray,
                     exp_intensities: np.ndarray,
                     database: Dict,
                     elements: Optional[List[str]] = None,
                     top_n: int = 10) -> List[Dict]:
        """
        Match experimental XRD pattern to database.
        
        Args:
            exp_positions: Experimental peak positions (2θ)
            exp_intensities: Experimental peak intensities
            database: XRD database dictionary
            elements: Optional list of required elements for filtering
            top_n: Number of top matches to return
            
        Returns:
            List of match results, sorted by score (descending)
        """
        xrd_db = database['xrd_database']
        
        # Filter by elements if provided
        if elements is not None:
            candidate_ids = self.filter_by_elements(database, elements)
        else:
            candidate_ids = list(xrd_db.keys())
        
        if len(candidate_ids) == 0:
            logger.warning("No candidate entries found")
            return []
        
        logger.info(f"Matching against {len(candidate_ids)} candidate entries...")
        
        # Calculate scores for all candidates
        results = []

        for entry_id in candidate_ids:
            entry = xrd_db[entry_id]

            score_result = self.match_single_entry(
                exp_positions,
                exp_intensities,
                entry
            )

            # Handle different scoring methods
            if self.scoring_method == 'combined':
                # score_result is a dict
                if score_result['weighted_score'] > 0 or score_result['fom'] > 0:
                    result = {
                        'entry_id': entry_id,
                        'mpid': entry['mpid'],
                        'formula': entry['formula'],
                        'elements': entry['elements'],
                        'score': score_result['weighted_score'],
                        'fom': score_result['fom'],
                        'n_matched_peaks': score_result['n_matched_peaks'],
                        'spacegroup': entry['spacegroup_number']
                    }
                    results.append(result)
            else:
                # score_result is a float
                if score_result > 0:
                    result = {
                        'entry_id': entry_id,
                        'mpid': entry['mpid'],
                        'formula': entry['formula'],
                        'elements': entry['elements'],
                        'spacegroup': entry['spacegroup_number']
                    }

                    # Add appropriate score field
                    if self.scoring_method == 'fom':
                        result['fom'] = score_result
                        result['score'] = score_result  # For backward compatibility
                    else:  # 'weighted'
                        result['score'] = score_result

                    results.append(result)
        
        # Sort by score (descending)
        results.sort(key=lambda x: x['score'], reverse=True)
        
        # Return top N
        top_results = results[:top_n]
        
        logger.info(f"Found {len(results)} matches, returning top {len(top_results)}")
        
        if len(top_results) > 0:
            logger.info(f"Best match: {top_results[0]['formula']} "
                       f"(score: {top_results[0]['score']:.2f})")
        
        return top_results
    
    def calculate_figure_of_merit(self,
                                  exp_positions: np.ndarray,
                                  exp_intensities: np.ndarray,
                                  db_positions: np.ndarray,
                                  db_intensities: np.ndarray) -> Dict[str, float]:
        """
        Calculate detailed figure of merit metrics.
        
        Args:
            exp_positions: Experimental peak positions
            exp_intensities: Experimental peak intensities
            db_positions: Database peak positions
            db_intensities: Database peak intensities
            
        Returns:
            Dictionary with various metrics
        """
        # Basic match score
        match_score = self.calculate_peak_match_score(
            exp_positions, exp_intensities,
            db_positions, db_intensities
        )
        
        # Calculate additional metrics
        metrics = {
            'match_score': match_score,
            'n_exp_peaks': len(exp_positions),
            'n_db_peaks': len(db_positions)
        }
        
        return metrics

