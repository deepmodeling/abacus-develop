"""
Module for reading and writing XYZ format files.

This module provides functions to read from and write to XYZ format files,
which are simple coordinate files used to specify molecular geometries.
"""

from typing import List, Tuple, Union
import numpy as np

def read_xyz(filepath: str) -> Tuple[List[str], np.ndarray]:
    """
    Read an XYZ format file.
    
    Args:
        filepath: Path to the XYZ file
        
    Returns:
        Tuple containing:
            - List of atomic symbols
            - Numpy array of atomic coordinates
            
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    species = []
    coords = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        num_atoms = int(lines[0])
        # Skip comment line
        for line in lines[2:2+num_atoms]:
            parts = line.split()
            species.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
    
    return species, np.array(coords)

def write_xyz(filepath: str, 
              species: List[str], 
              coords: Union[List[List[float]], np.ndarray], 
              comment: str = '') -> None:
    """
    Write structure information to an XYZ file.
    
    Args:
        filepath: Path where the file should be written
        species: List of atomic symbols
        coords: Array of atomic coordinates
        comment: Optional comment for the second line
        
    Raises:
        ValueError: If the input data is invalid
        OSError: If there are problems writing the file
    """
    coords = np.asarray(coords)
    if len(species) != len(coords):
        raise ValueError("Number of species must match number of coordinate sets")
        
    with open(filepath, 'w') as f:
        f.write(f"{len(species)}\n")
        f.write(f"{comment}\n")
        for sym, pos in zip(species, coords):
            f.write(f"{sym:2s} {pos[0]:15.8f} {pos[1]:15.8f} {pos[2]:15.8f}\n")