"""
XRD Database Builder Module
============================

This module builds a searchable database of XRD patterns from crystal structure data.
"""

import json
import pickle
import logging
from pathlib import Path
from typing import Dict, List, Optional
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial
import gc

from .db_processor import DatabaseProcessor
from .xrd_calculator import XRDCalculator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def process_entry_worker(args):
    """
    Worker function for parallel processing of database entries.

    Args:
        args: Tuple of (entry_id, info, wavelength, n_peaks, two_theta_range)

    Returns:
        Tuple of (entry_id, entry_record) or (entry_id, None) if failed
    """
    entry_id, info, wavelength, n_peaks, two_theta_range = args

    try:
        # Initialize XRD calculator for this worker
        xrd_calc = XRDCalculator(
            wavelength=wavelength,
            two_theta_range=two_theta_range
        )

        # Calculate XRD pattern and extract peaks
        atoms = info['atoms']
        peaks = xrd_calc.calculate_and_extract_peaks(
            atoms,
            n_peaks=n_peaks,
            normalize=True
        )

        if peaks is None or len(peaks) == 0:
            return (entry_id, None)

        # Format peaks for storage
        peak_data = xrd_calc.format_peaks_for_storage(peaks)

        # Create entry record
        entry_record = {
            'entry_id': entry_id,
            'mpid': info['mpid'],
            'formula': info['formula'],
            'elements': info['elements'],
            'n_atoms': info['n_atoms'],
            'spacegroup_number': info['spacegroup_number'],
            'spacegroup_symbol': info['spacegroup_symbol'],
            'peaks': peak_data
        }

        return (entry_id, entry_record)

    except Exception as e:
        logger.warning(f"Worker failed to process entry {entry_id}: {e}")
        return (entry_id, None)


class DatabaseBuilder:
    """
    Build a searchable XRD pattern database from crystal structures.
    
    This class processes a crystal structure database, calculates XRD patterns,
    extracts top peaks, and creates an indexed database for fast retrieval.
    """
    
    def __init__(self, 
                 db_path: str,
                 wavelength: str = 'CuKa',
                 n_peaks: int = 5,
                 two_theta_range: tuple = (10, 90)):
        """
        Initialize the database builder.
        
        Args:
            db_path: Path to the ASE crystal structure database
            wavelength: X-ray wavelength for XRD calculation
            n_peaks: Number of top peaks to extract per structure
            two_theta_range: Range of 2θ angles to calculate
        """
        self.db_path = db_path
        self.wavelength = wavelength
        self.n_peaks = n_peaks
        self.two_theta_range = two_theta_range
        
        # Initialize processors
        self.db_processor = DatabaseProcessor(db_path)
        self.xrd_calculator = XRDCalculator(
            wavelength=wavelength,
            two_theta_range=two_theta_range
        )
        
        # Storage for processed data
        self.xrd_database = {}
        self.element_index = {}  # Index by elements for fast lookup
        
        logger.info(f"DatabaseBuilder initialized for {db_path}")
    
    def process_single_entry(self, entry_id: int) -> Optional[Dict]:
        """
        Process a single database entry to extract XRD peaks.
        
        Args:
            entry_id: Database entry ID
            
        Returns:
            Dictionary with entry information and XRD peaks, or None if failed
        """
        # Get structure information
        info = self.db_processor.get_structure_info(entry_id)
        
        if info is None:
            return None
        
        # Calculate XRD pattern and extract peaks
        atoms = info['atoms']
        peaks = self.xrd_calculator.calculate_and_extract_peaks(
            atoms, 
            n_peaks=self.n_peaks,
            normalize=True
        )
        
        if peaks is None or len(peaks) == 0:
            logger.warning(f"No peaks found for entry {entry_id}")
            return None
        
        # Format peaks for storage
        peak_data = self.xrd_calculator.format_peaks_for_storage(peaks)
        
        # Create entry record
        entry_record = {
            'entry_id': entry_id,
            'mpid': info['mpid'],
            'formula': info['formula'],
            'elements': info['elements'],
            'n_atoms': info['n_atoms'],
            'spacegroup_number': info['spacegroup_number'],
            'spacegroup_symbol': info['spacegroup_symbol'],
            'peaks': peak_data
        }
        
        return entry_record
    
    def build_database(self,
                      output_path: str = 'xrd_database.pkl',
                      max_entries: Optional[int] = None,
                      skip_errors: bool = True,
                      batch_size: int = 5000,
                      save_interval: int = 1000,
                      save_json: bool = False) -> Dict:
        """
        Build the complete XRD database from all entries with batch processing.

        Args:
            output_path: Path to save the database
            max_entries: Maximum number of entries to process (None for all)
            skip_errors: Whether to skip entries that fail processing
            batch_size: Number of entries to process before clearing memory
            save_interval: Save database every N entries

        Returns:
            Dictionary containing the XRD database
        """
        logger.info("Starting database build process...")

        total_entries = self.db_processor.get_total_entries()
        if max_entries is not None:
            total_entries = min(total_entries, max_entries)

        logger.info(f"Processing {total_entries} entries...")
        logger.info(f"Batch size: {batch_size}, Save interval: {save_interval}")

        processed_count = 0
        failed_count = 0

        # Load existing database if it exists
        output_file = Path(output_path)
        if output_file.exists():
            logger.info(f"Loading existing database from {output_path}")
            try:
                with open(output_path, 'rb') as f:
                    saved_data = pickle.load(f)
                    self.xrd_database = saved_data.get('xrd_database', {})
                    self.element_index = saved_data.get('element_index', {})
                    processed_count = len(self.xrd_database)
                    logger.info(f"Loaded {processed_count} existing entries")
            except Exception as e:
                logger.warning(f"Failed to load existing database: {e}")

        # Process all entries with progress bar
        with tqdm(total=total_entries, initial=processed_count,
                  desc="Building XRD database", unit="entry") as pbar:

            for i, (entry_id, info) in enumerate(self.db_processor.iterate_all_entries()):
                if max_entries is not None and i >= max_entries:
                    break

                # Skip if already processed
                if entry_id in self.xrd_database:
                    pbar.update(1)
                    continue

                try:
                    # Calculate XRD pattern and extract peaks
                    atoms = info['atoms']
                    peaks = self.xrd_calculator.calculate_and_extract_peaks(
                        atoms,
                        n_peaks=self.n_peaks,
                        normalize=True
                    )

                    if peaks is None or len(peaks) == 0:
                        if not skip_errors:
                            raise ValueError(f"No peaks found for entry {entry_id}")
                        failed_count += 1
                        pbar.update(1)
                        continue

                    # Format peaks for storage
                    peak_data = self.xrd_calculator.format_peaks_for_storage(peaks)

                    # Create entry record
                    entry_record = {
                        'entry_id': entry_id,
                        'mpid': info['mpid'],
                        'formula': info['formula'],
                        'elements': info['elements'],
                        'n_atoms': info['n_atoms'],
                        'spacegroup_number': info['spacegroup_number'],
                        'spacegroup_symbol': info['spacegroup_symbol'],
                        'peaks': peak_data
                    }

                    # Store in database
                    self.xrd_database[entry_id] = entry_record

                    # Update element index
                    element_key = tuple(sorted(info['elements']))
                    if element_key not in self.element_index:
                        self.element_index[element_key] = []
                    self.element_index[element_key].append(entry_id)

                    processed_count += 1
                    pbar.update(1)

                    # Periodic save
                    if processed_count % save_interval == 0:
                        self.save_database(output_path, save_json=False)  # Don't save JSON for checkpoints
                        logger.info(f"Checkpoint saved: {processed_count} entries")

                except Exception as e:
                    if not skip_errors:
                        raise
                    logger.warning(f"Failed to process entry {entry_id}: {e}")
                    failed_count += 1
                    pbar.update(1)

        logger.info(f"Database build complete: {processed_count} successful, {failed_count} failed")

        # Final save
        self.save_database(output_path, save_json=save_json)
        logger.info(f"Final database saved to {output_path}")

        return self.xrd_database

    def build_database_parallel(self,
                               output_path: str = 'xrd_database.pkl',
                               max_entries: Optional[int] = None,
                               skip_errors: bool = True,
                               n_workers: Optional[int] = None,
                               batch_size: int = 5000,
                               save_interval: int = 1000,
                               save_json: bool = False) -> Dict:
        """
        Build XRD database using parallel processing (multiprocessing).

        Args:
            output_path: Path to save the database
            max_entries: Maximum number of entries to process (None for all)
            skip_errors: Whether to skip entries that fail processing
            n_workers: Number of parallel workers (None = use all CPU cores)
            batch_size: Number of entries to process in each batch
            save_interval: Save database every N entries

        Returns:
            Dictionary containing the XRD database
        """
        if n_workers is None:
            n_workers = cpu_count()

        logger.info(f"Starting parallel database build with {n_workers} workers...")

        total_entries = self.db_processor.get_total_entries()
        if max_entries is not None:
            total_entries = min(total_entries, max_entries)

        logger.info(f"Processing {total_entries} entries...")
        logger.info(f"Workers: {n_workers}, Batch size: {batch_size}, Save interval: {save_interval}")

        processed_count = 0
        failed_count = 0

        # Load existing database if it exists
        output_file = Path(output_path)
        if output_file.exists():
            logger.info(f"Loading existing database from {output_path}")
            try:
                with open(output_path, 'rb') as f:
                    saved_data = pickle.load(f)
                    self.xrd_database = saved_data.get('xrd_database', {})
                    self.element_index = saved_data.get('element_index', {})
                    processed_count = len(self.xrd_database)
                    logger.info(f"Loaded {processed_count} existing entries")
            except Exception as e:
                logger.warning(f"Failed to load existing database: {e}")

        # Collect entries to process
        entries_to_process = []

        logger.info("Collecting entries to process...")
        for i, (entry_id, info) in enumerate(self.db_processor.iterate_all_entries()):
            if max_entries is not None and i >= max_entries:
                break

            # Skip if already processed
            if entry_id in self.xrd_database:
                continue

            entries_to_process.append((
                entry_id,
                info,
                self.wavelength,
                self.n_peaks,
                self.two_theta_range
            ))

        logger.info(f"Found {len(entries_to_process)} entries to process")

        if len(entries_to_process) == 0:
            logger.info("No new entries to process")
            return self.xrd_database

        # Process in batches with multiprocessing
        total_to_process = len(entries_to_process)

        with tqdm(total=total_to_process, desc="Building XRD database (parallel)", unit="entry") as pbar:

            for batch_start in range(0, total_to_process, batch_size):
                batch_end = min(batch_start + batch_size, total_to_process)
                batch = entries_to_process[batch_start:batch_end]

                logger.info(f"Processing batch {batch_start}-{batch_end} with {n_workers} workers...")

                # Process batch in parallel
                with Pool(processes=n_workers) as pool:
                    results = pool.map(process_entry_worker, batch)

                # Store results
                for entry_id, entry_record in results:
                    if entry_record is not None:
                        # Store in database
                        self.xrd_database[entry_id] = entry_record

                        # Update element index
                        element_key = tuple(sorted(entry_record['elements']))
                        if element_key not in self.element_index:
                            self.element_index[element_key] = []
                        self.element_index[element_key].append(entry_id)

                        processed_count += 1
                    else:
                        failed_count += 1

                    pbar.update(1)

                # Save checkpoint after each batch
                self.save_database(output_path, save_json=False)  # Don't save JSON for checkpoints
                logger.info(f"Checkpoint saved: {processed_count} entries processed, {failed_count} failed")

                # Clear memory
                gc.collect()

        logger.info(f"Parallel database build complete: {processed_count} successful, {failed_count} failed")

        # Final save
        self.save_database(output_path, save_json=save_json)
        logger.info(f"Final database saved to {output_path}")

        return self.xrd_database
    
    def save_database(self, output_path: str, save_json: bool = False):
        """
        Save the XRD database to disk.

        Args:
            output_path: Path to save the database
            save_json: Whether to also save as JSON (default: False, not recommended for large databases)
        """
        database_package = {
            'xrd_database': self.xrd_database,
            'element_index': self.element_index,
            'metadata': {
                'source_db': self.db_path,
                'wavelength': self.wavelength,
                'n_peaks': self.n_peaks,
                'two_theta_range': self.two_theta_range,
                'total_entries': len(self.xrd_database)
            }
        }

        # Save as pickle
        with open(output_path, 'wb') as f:
            pickle.dump(database_package, f, protocol=pickle.HIGHEST_PROTOCOL)

        logger.info(f"Database saved to {output_path}")

        # Optionally save as JSON for human readability (not recommended for large databases)
        if save_json:
            json_path = Path(output_path).with_suffix('.json')
            try:
                logger.info(f"Saving JSON version (this may take a while for large databases)...")
                with open(json_path, 'w') as f:
                    json.dump(database_package, f, indent=2)
                logger.info(f"Database also saved as JSON to {json_path}")
            except Exception as e:
                logger.warning(f"Could not save JSON version: {e}")
    
    @staticmethod
    def load_database(database_path: str) -> Dict:
        """
        Load a previously built XRD database.
        
        Args:
            database_path: Path to the database file
            
        Returns:
            Dictionary containing the database
        """
        with open(database_path, 'rb') as f:
            database_package = pickle.load(f)
        
        logger.info(f"Loaded database from {database_path}")
        logger.info(f"Total entries: {database_package['metadata']['total_entries']}")
        
        return database_package

