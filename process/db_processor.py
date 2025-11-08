"""
Database Processor Module
=========================

This module handles reading and processing crystal structure data from ASE database files.
"""

import logging
from typing import Dict, List, Optional, Tuple
from ase.db import connect
from ase.spacegroup import get_spacegroup
from ase import Atoms

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DatabaseProcessor:
    """
    Process crystal structure database files.
    
    This class provides methods to read crystal structures from ASE database files
    and extract relevant information including chemical composition, space group,
    and atomic positions.
    """
    
    def __init__(self, db_path: str):
        """
        Initialize the database processor.
        
        Args:
            db_path: Path to the ASE database file (.db)
        """
        self.db_path = db_path
        self.db = None
        self._connect()
    
    def _connect(self):
        """Establish connection to the database."""
        try:
            self.db = connect(self.db_path)
            logger.info(f"Successfully connected to database: {self.db_path}")
            logger.info(f"Total entries in database: {len(self.db)}")
        except Exception as e:
            logger.error(f"Failed to connect to database {self.db_path}: {e}")
            raise
    
    def get_total_entries(self) -> int:
        """
        Get the total number of entries in the database.
        
        Returns:
            Number of entries
        """
        return len(self.db)
    
    def get_atoms(self, entry_id: int) -> Optional[Atoms]:
        """
        Get atoms object for a specific entry.
        
        Args:
            entry_id: Database entry ID
            
        Returns:
            ASE Atoms object or None if not found
        """
        try:
            atoms = self.db.get_atoms(id=entry_id)
            return atoms
        except Exception as e:
            logger.warning(f"Failed to get atoms for entry {entry_id}: {e}")
            return None
    
    def get_structure_info(self, entry_id: int) -> Optional[Dict]:
        """
        Extract comprehensive structure information for a database entry.
        
        Args:
            entry_id: Database entry ID
            
        Returns:
            Dictionary containing structure information or None if failed
        """
        try:
            row = self.db.get(id=entry_id)
            atoms = row.toatoms()
            
            # Extract basic information
            chemical_symbols = atoms.get_chemical_symbols()
            unique_elements = sorted(list(set(chemical_symbols)))
            
            # Get space group
            try:
                spacegroup = get_spacegroup(atoms)
                spacegroup_number = spacegroup.no
                spacegroup_symbol = spacegroup.symbol
            except Exception as e:
                logger.warning(f"Could not determine space group for entry {entry_id}: {e}")
                spacegroup_number = None
                spacegroup_symbol = None
            
            # Get structural data
            cell = atoms.get_cell()
            scaled_positions = atoms.get_scaled_positions()
            cartesian_positions = atoms.get_positions()
            
            # Get material project ID if available
            mpid = row.get('mpid', None)
            if mpid is None:
                # Try alternative keys
                mpid = row.get('material_id', None)
                if mpid is None:
                    mpid = row.get('id', None)
            
            # Get chemical formula
            formula = atoms.get_chemical_formula()
            
            info = {
                'entry_id': entry_id,
                'mpid': mpid,
                'formula': formula,
                'elements': unique_elements,
                'n_atoms': len(atoms),
                'chemical_symbols': chemical_symbols,
                'spacegroup_number': spacegroup_number,
                'spacegroup_symbol': spacegroup_symbol,
                'cell': cell,
                'scaled_positions': scaled_positions,
                'cartesian_positions': cartesian_positions,
                'atoms': atoms
            }
            
            return info
            
        except Exception as e:
            logger.error(f"Failed to extract structure info for entry {entry_id}: {e}")
            return None
    
    def iterate_all_entries(self):
        """
        Generator to iterate through all database entries.
        
        Yields:
            Tuple of (entry_id, structure_info)
        """
        total = self.get_total_entries()
        logger.info(f"Starting iteration through {total} entries")
        
        for i, row in enumerate(self.db.select(), start=1):
            entry_id = row.id
            info = self.get_structure_info(entry_id)
            
            if info is not None:
                yield entry_id, info
            
            if i % 100 == 0:
                logger.info(f"Processed {i}/{total} entries ({100*i/total:.1f}%)")
    
    def get_entries_by_elements(self, elements: List[str]) -> List[int]:
        """
        Find all entries containing specific elements.
        
        Args:
            elements: List of element symbols (e.g., ['Al', 'O'])
            
        Returns:
            List of entry IDs
        """
        matching_ids = []
        
        for entry_id, info in self.iterate_all_entries():
            entry_elements = set(info['elements'])
            query_elements = set(elements)
            
            # Check if all query elements are present
            if query_elements.issubset(entry_elements):
                matching_ids.append(entry_id)
        
        logger.info(f"Found {len(matching_ids)} entries containing elements {elements}")
        return matching_ids
    
    def close(self):
        """Close the database connection."""
        if self.db is not None:
            # ASE database doesn't require explicit closing
            self.db = None
            logger.info("Database connection closed")

