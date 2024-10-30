"""
Module for reading and writing ABACUS STRU files.

This module provides functions to read from and write to STRU files used by ABACUS.
It supports various formats and features including different coordinate types,
magnetic moments, movement constraints, and initial velocities.
"""

from typing import Dict, Any
import numpy as np
import os

# Move the content from struio.py here
# Add type hints and improve documentation

def _parse_coord_line(line):
    """
    Parses a coordinate line in the ATOMIC_POSITIONS block.
    
    Args:
        line: String containing atomic coordinates and optional parameters
        
    Returns:
        Dictionary containing parsed coordinate data and optional parameters
    """
    fields = line.strip().split()
    if len(fields) < 3:
        raise ValueError(f"Invalid coordinate line: '{line}'. Expected at least 3 coordinates.")
    
    result = {'coord': [float(x) for x in fields[0:3]]}
    
    idx = 3
    while idx < len(fields):
        if fields[idx].isdigit(): # no keyword, 0/1 -> frozen atom
            result['m'] = [int(x) for x in fields[idx:idx+3]]
            idx += 3
        elif fields[idx] == 'm': # frozen atom
            result['m'] = [int(x) for x in fields[idx+1:idx+4]]
            idx += 4
        elif fields[idx] in ['v', 'vel', 'velocity']: # initial velocity
            result['v'] = [float(x) for x in fields[idx+1:idx+4]]
            idx += 4
        elif fields[idx] in ['mag', 'magmom']:
            if idx + 2 < len(fields) and fields[idx+2] == 'angle1':
                result['mag'] = ('Spherical',
                                [float(fields[idx+1]), float(fields[idx+3]), float(fields[idx+5])])
                idx += 6
            elif idx + 2 < len(fields) and fields[idx+2][0].isdigit():
                result['mag'] = ('Cartesian',
                                [float(x) for x in fields[idx+1:idx+4]])
                idx += 4
            else: # collinear
                result['mag'] = float(fields[idx+1])
                idx += 2
        else:
            raise ValueError(f'Unknown keyword: {fields[idx]} in line: {line}')
    
    return result


def _atomic_positions_gen(lines):
    """
    Iteratively generates info per species from the ATOMIC_POSITIONS block.
    
    Args:
        lines: List of lines from the ATOMIC_POSITIONS block
        
    Yields:
        Dictionary containing species information and atomic positions
    """
    current_pos = 0
    while current_pos < len(lines):
        # Read species symbol
        symbol = lines[current_pos].strip()
        
        # Read magnetic moment and number of atoms
        try:
            mag_each = float(lines[current_pos + 1])
            natom = int(lines[current_pos + 2])
        except (IndexError, ValueError) as e:
            raise ValueError(f"Error reading atomic positions for {symbol}: {e}")
        
        # Read atomic positions
        pos_start = current_pos + 3
        pos_end = pos_start + natom
        if pos_end > len(lines):
            raise ValueError(f"Not enough coordinate lines for {symbol}. Expected {natom}, found {len(lines) - pos_start}")
        
        # Parse coordinates for each atom
        try:
            atoms = [_parse_coord_line(line) for line in lines[pos_start:pos_end]]
        except Exception as e:
            raise ValueError(f"Error parsing coordinates for {symbol}: {e}")
        
        # Yield the species information
        yield {
            'symbol': symbol,
            'mag_each': mag_each,
            'natom': natom,
            'atom': atoms
        }
        
        # Move to next species block
        current_pos = pos_end

def read_stru(fpath: str) -> Dict[str, Any]:
    """Read an ABACUS STRU file and return its content as a dictionary."""
    block_title = ['ATOMIC_SPECIES',
                   'NUMERICAL_ORBITAL',
                   'LATTICE_CONSTANT',
                   'LATTICE_PARAMETER',
                   'LATTICE_VECTORS',
                   'ATOMIC_POSITIONS']

    def _trim(line):
        return line.split('#')[0].split('//')[0].strip(' \t\n')

    with open(fpath, 'r') as f:
        lines = [_trim(line).replace('\t', ' ')
                 for line in f.readlines() if len(_trim(line)) > 0]

    # break the content into blocks
    delim = [i for i, line in enumerate(lines) if line in block_title] \
            + [len(lines)]
    blocks = {lines[delim[i]] : lines[delim[i]+1:delim[i+1]]
              for i in range(len(delim) - 1)}

    stru = {}
    #============ LATTICE_CONSTANT/PARAMETER/VECTORS ============
    stru['lat'] = {'const': float(blocks['LATTICE_CONSTANT'][0])}
    if 'LATTICE_VECTORS' in blocks:
        stru['lat']['vec'] = [[float(x) for x in line.split()]
                              for line in blocks['LATTICE_VECTORS']]
    elif 'LATTICE_PARAMETER' in blocks:
        stru['lat']['param'] = [float(x)
                                for x in blocks['LATTICE_PARAMETERS'].split()]

    #============ ATOMIC_SPECIES ============
    stru['species'] = [_atomic_species_from_file(line)
                      for line in blocks['ATOMIC_SPECIES']]

    #============ NUMERICAL_ORBITAL ============
    if 'NUMERICAL_ORBITAL' in blocks:
        for i, s in enumerate(stru['species']):
            s['orb_file'] = blocks['NUMERICAL_ORBITAL'][i].strip()

    #============ ATOMIC_POSITIONS ============
    stru['coord_type'] = blocks['ATOMIC_POSITIONS'][0]
    index = {s['symbol']: i for i, s in enumerate(stru['species'])}
    try:
        for ap in _atomic_positions_gen(blocks['ATOMIC_POSITIONS'][1:]):
            stru['species'][index[ap['symbol']]].update(ap)
    except Exception as e:
        raise ValueError(f"Error processing ATOMIC_POSITIONS: {e}")

    return stru

def _atomic_species_from_file(line):
    """
    Process ATOMIC_SPECIES line to handle paths correctly.
    Stores clean filename but preserves original format for writing.
    """
    fields = line.split()
    # Remove './' prefix if present and extract just the filename
    pp_file = fields[2].replace('./', '', 1) if fields[2] else ''
    return {
        'symbol': fields[0],
        'mass': float(fields[1]),
        'pp_file': pp_file,  # Store clean filename without './'
        'pp_type': fields[3] if len(fields) > 3 else None
    }

def write_stru(job_dir: str, stru: Dict[str, Any], fname: str = 'STRU') -> None:
    """
    Write structure information to an ABACUS STRU file.
    
    Args:
        job_dir: Directory where the file should be written
        stru: Dictionary containing the structure information
        fname: Name of the output file (default: 'STRU')
        
    Raises:
        ValueError: If the structure dictionary is invalid
        OSError: If there are problems writing the file
    """
    # Copy the write_stru implementation from struio.py
    with open(os.path.join(job_dir, fname), 'w') as f:

        #============ ATOMIC_SPECIES ============
        f.write('ATOMIC_SPECIES\n')
        # Calculate width based on actual format that will be written
        width = {
            'symbol_w': max(len(s['symbol']) for s in stru['species']),
            'mass_w': max(len(str(s['mass'])) for s in stru['species']),
            'pp_file_w': max(len('./' + s['pp_file']) for s in stru['species'])
        }
        for s in stru['species']:
            # Ensure clean filename and add './' prefix
            pp_file = s['pp_file']
            if pp_file.startswith('./'):
                pp_file = pp_file[2:]
            pp_path = './' + pp_file

            f.write('{symbol:<{symbol_w}}  {mass:>{mass_w}}  {pp_file:>{pp_file_w}}'.format(
                symbol=s['symbol'],
                mass=s['mass'],
                pp_file=pp_path,
                **width))
            if 'pp_type' in s and s['pp_type']:
                f.write(f"  {s['pp_type']}")
            f.write('\n')

        #============ NUMERICAL_ORBITAL ============
        if 'orb_file' in stru['species'][0]:
            f.write('\nNUMERICAL_ORBITAL\n')
            for s in stru['species']:
                f.write(f"{s['orb_file']}\n")
        
        #============ LATTICE_CONSTANT/PARAMETER/VECTORS ============
        f.write('\nLATTICE_CONSTANT\n')
        f.write(f"{stru['lat']['const']}\n")

        if 'vec' in stru['lat']:
            f.write('\nLATTICE_VECTORS\n')
            for v in stru['lat']['vec']:
                f.write(f'{v[0]} {v[1]} {v[2]}\n')

        if 'param' in stru['lat']:
            f.write('\nLATTICE_PARAMETER\n')
            for param in stru['lat']['param']:
                f.write(f'{param} ')
            f.write('\n')

        #============ ATOMIC_POSITIONS ============
        f.write('\nATOMIC_POSITIONS\n')
        f.write(f"{stru['coord_type']}\n")

        for s in stru['species']:
            f.write(f"\n{s['symbol']}\n")
            f.write(f"{s['mag_each']}\n")
            f.write(f"{s['natom']}\n")

            for atom in s['atom']:
                f.write(' '.join(f'{x}' for x in atom['coord']))

                for key in ['m', 'v']: # frozen atom / initial velocity
                    if key in atom:
                        f.write(f' {key}' +
                                ''.join(f' {x}' for x in atom[key]))

                if 'mag' in atom:
                    if not isinstance(atom['mag'], tuple): # collinear
                        f.write(f" mag {atom['mag']}")
                    else: # non-collinear
                        mag_coord_type, mag = atom['mag']
                        assert mag_coord_type in ['Cartesian', 'Spherical']
                        if mag_coord_type == 'Cartesian':
                            f.write(f' mag {mag[0]} {mag[1]} {mag[2]}')
                        else:
                            f.write(f' mag {mag[0]} angle1 {mag[1]} angle2 {mag[2]}')

                f.write('\n')