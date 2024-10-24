"""
Utility functions for input/output operations.

This module provides common utility functions used across different IO modules.
"""

from typing import List
import numpy as np

def cart_to_direct(cart_coords: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    Convert Cartesian coordinates to direct (fractional) coordinates.
    
    Args:
        cart_coords: Cartesian coordinates
        lattice: Lattice vectors (3x3 matrix)
        
    Returns:
        Direct coordinates
    """
    return np.dot(cart_coords, np.linalg.inv(lattice))

def direct_to_cart(direct_coords: np.ndarray, lattice: np.ndarray) -> np.ndarray:
    """
    Convert direct (fractional) coordinates to Cartesian coordinates.
    
    Args:
        direct_coords: Direct coordinates
        lattice: Lattice vectors (3x3 matrix)
        
    Returns:
        Cartesian coordinates
    """
    return np.dot(direct_coords, lattice)

def parse_bool_vector(vector_str: str) -> List[bool]:
    """
    Parse a string of 1s and 0s into a list of booleans.
    
    Args:
        vector_str: String containing 1s and 0s
        
    Returns:
        List of boolean values
        
    Raises:
        ValueError: If the string contains invalid characters
    """
    return [bool(int(x)) for x in vector_str.split()]