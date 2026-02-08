#!/usr/bin/env python3
"""
Generate markdown documentation from C++ source files.

This script parses read_input_item_*.cpp files to extract documentation
from Input_Item registrations and generates input-keywords-generated.md.

Usage:
    python generate_docs_from_source.py [--source-dir DIR] [--output FILE]
"""

import argparse
import re
import sys
from pathlib import Path
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple


def extract_string_value(text: str, field: str) -> Optional[str]:
    """
    Extract a string field value from C++ code.
    Handles both regular strings and raw string literals R"(...)".
    Also handles string concatenation across lines.
    """
    # Pattern for raw string literal: item.field = R"(...)"
    raw_pattern = rf'item\.{field}\s*=\s*R"\((.*?)\)"'
    raw_match = re.search(raw_pattern, text, re.DOTALL)
    if raw_match:
        return raw_match.group(1)

    # Pattern for simple string: item.field = "..."
    # Need to handle string concatenation: "..." "..."
    simple_pattern = rf'item\.{field}\s*=\s*("(?:[^"\\]|\\.)*"(?:\s*"(?:[^"\\]|\\.)*")*)'
    simple_match = re.search(simple_pattern, text, re.DOTALL)
    if simple_match:
        # Concatenate all string parts
        full_str = simple_match.group(1)
        # Extract individual strings and concatenate
        parts = re.findall(r'"((?:[^"\\]|\\.)*)"', full_str)
        result = ''.join(parts)
        # Unescape common escape sequences
        result = result.replace('\\n', '\n')
        result = result.replace('\\t', '\t')
        result = result.replace('\\"', '"')
        result = result.replace('\\\\', '\\')
        return result

    return None


def parse_input_item_block(block: str) -> Optional[Dict[str, str]]:
    """
    Parse a single Input_Item block and extract all documentation fields.
    """
    # Extract parameter name
    name_match = re.search(r'Input_Item\s+\w+\s*\(\s*"([^"]+)"\s*\)', block)
    if not name_match:
        return None

    param = {'name': name_match.group(1)}

    # Extract documentation fields
    fields = ['category', 'type', 'description', 'default_value', 'unit', 'availability']
    for field in fields:
        value = extract_string_value(block, field)
        if value is not None:
            param[field] = value
        else:
            param[field] = ''

    return param


def extract_input_items(content: str) -> List[Dict[str, str]]:
    """
    Extract all Input_Item blocks from a C++ file content.
    """
    items = []

    # Find all blocks between { Input_Item ... this->add_item(item); }
    # Use a simpler approach: find Input_Item declarations and their surrounding braces

    # Split by Input_Item declarations
    # Allow C++ comments between { and Input_Item declaration
    pattern = r'\{(?:\s|//[^\n]*\n)*Input_Item\s+(\w+)\s*\(\s*"([^"]+)"\s*\)'

    # Find start positions of each Input_Item block
    starts = []
    for match in re.finditer(pattern, content):
        # Find the opening brace before Input_Item
        start = content.rfind('{', 0, match.start() + 1)
        if start != -1:
            starts.append((start, match.group(2)))  # (position, param_name)

    # For each start, find the matching closing brace
    for start_pos, param_name in starts:
        depth = 0
        end_pos = start_pos
        for i in range(start_pos, len(content)):
            if content[i] == '{':
                depth += 1
            elif content[i] == '}':
                depth -= 1
                if depth == 0:
                    end_pos = i + 1
                    break

        block = content[start_pos:end_pos]

        # Parse the block
        param = parse_input_item_block(block)
        if param and param.get('category'):  # Only include items with documentation
            items.append(param)

    return items


def format_description(desc: str) -> str:
    """
    Format description text for markdown output.
    - Convert * list markers to - list markers
    - Ensure proper indentation for nested lists
    - Convert [NOTE] markers to blockquotes
    """
    if not desc:
        return ''

    lines = desc.split('\n')
    result_lines = []

    for line in lines:
        # Convert * list items to - list items
        line = re.sub(r'^(\s*)\*\s+', r'\1- ', line)

        # Convert [NOTE] markers to blockquotes
        if '[NOTE]' in line:
            line = line.replace('[NOTE]', '> Note:')

        # Convert [WARNING] markers to blockquotes
        if '[WARNING]' in line:
            line = line.replace('[WARNING]', '> Warning:')

        result_lines.append(line)

    # Join and clean up
    result = '\n'.join(result_lines)

    # Ensure list items have proper indentation (2 spaces for sub-items in markdown)
    result = re.sub(r'\n- ', '\n  - ', result)
    # But not for the first item after a non-list line
    result = re.sub(r'(\n[^-\s][^\n]*)\n  - ', r'\1\n\n  - ', result)

    return result.strip()


def generate_parameter_markdown(param: Dict[str, str]) -> str:
    """
    Generate markdown for a single parameter.
    """
    lines = [f"### {param['name']}", ""]

    # Type
    if param.get('type'):
        lines.append(f"- **Type**: {param['type']}")

    # Availability (before description, as in original format)
    if param.get('availability'):
        lines.append(f"- **Availability**: *{param['availability']}*")

    # Description
    if param.get('description'):
        desc = format_description(param['description'])
        # If description has multiple lines/lists, format properly
        if '\n' in desc:
            lines.append(f"- **Description**: {desc.split(chr(10))[0]}")
            for line in desc.split('\n')[1:]:
                if line.strip():
                    lines.append(f"  {line}" if not line.startswith('  ') else line)
                else:
                    lines.append("")
        else:
            lines.append(f"- **Description**: {desc}")

    # Default
    if param.get('default_value'):
        lines.append(f"- **Default**: {param['default_value']}")

    # Unit
    if param.get('unit'):
        lines.append(f"- **Unit**: {param['unit']}")

    lines.append("")
    return '\n'.join(lines)


def generate_category_markdown(category: str, params: List[Dict[str, str]]) -> str:
    """
    Generate markdown for a category section.
    """
    lines = [f"## {category}", ""]

    for param in params:
        lines.append(generate_parameter_markdown(param))

    return '\n'.join(lines)


# Define the desired category order
CATEGORY_ORDER = [
    "System variables",
    "Input files",
    "Plane wave related variables",
    "Numerical atomic orbitals related variables",
    "Electronic structure",
    "Electronic structure (SDFT)",
    "Geometry relaxation",
    "Output information",
    "Density of states",
    "NAOs",
    "DeePKS",
    "OFDFT: orbital free density functional theory",
    "ML-KEDF: machine learning based kinetic energy density functional for OFDFT",
    "TDOFDFT: time dependent orbital free density functional theory",
    "Electric field and dipole correction",
    "Gate field (compensating charge)",
    "Exact Exchange (Common)",
    "Exact Exchange (LCAO in PW)",
    "Exact Exchange (LCAO)",
    "Exact Exchange (PW)",
    "Molecular dynamics",
    "DFT+U correction",
    "Spin-Constrained DFT",
    "vdW correction",
    "Berry phase and wannier90 interface",
    "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory",
    "Variables useful for debugging",
    "Electronic conductivities",
    "Implicit solvation model",
    "Quasiatomic Orbital (QO) analysis",
    "PEXSI",
    "Linear Response TDDFT",
    "Linear Response TDDFT (Under Development Feature)",
    "Reduced Density Matrix Functional Theory",
    "Model",
    "Other",
]


def generate_anchor(text: str) -> str:
    """
    Generate a markdown anchor from text.
    Converts to lowercase, replaces spaces with hyphens, removes special chars.
    """
    anchor = text.lower()
    anchor = re.sub(r'[^a-z0-9\s_-]', '', anchor)
    anchor = re.sub(r'\s+', '-', anchor)
    anchor = re.sub(r'-+', '-', anchor)
    return anchor.strip('-')


def generate_toc(sorted_categories: OrderedDict) -> str:
    """
    Generate a markdown table of contents matching the original format.
    """
    lines = ["- [Full List of INPUT Keywords](#full-list-of-input-keywords)"]

    for category, params in sorted_categories.items():
        cat_anchor = generate_anchor(category)
        lines.append(f"  - [{category}](#{cat_anchor})")

        for param in params:
            # Escape underscores in parameter names for TOC display
            display_name = param['name'].replace('_', r'\_')
            param_anchor = generate_anchor(param['name'])
            lines.append(f"    - [{display_name}](#{param_anchor})")

    return '\n'.join(lines)


def generate(source_dir: Path, output: Path, verbose: bool = False):
    """
    Core generation logic. Can be called from conf.py or CLI.
    """
    # Find all read_input_item_*.cpp files
    source_files = sorted(source_dir.glob('read_input_item_*.cpp'))

    if not source_files:
        print(f"Error: No read_input_item_*.cpp files found in {source_dir}")
        sys.exit(1)

    if verbose:
        print(f"Found {len(source_files)} source files")

    # Extract all parameters
    all_params = []
    for source_file in source_files:
        if verbose:
            print(f"Processing {source_file.name}...")

        content = source_file.read_text()
        items = extract_input_items(content)
        all_params.extend(items)

        if verbose:
            print(f"  Found {len(items)} documented parameters")

    print(f"Total: {len(all_params)} documented parameters")

    # Group by category
    by_category: Dict[str, List[Dict[str, str]]] = OrderedDict()
    for param in all_params:
        cat = param.get('category', 'Other')
        if cat not in by_category:
            by_category[cat] = []
        by_category[cat].append(param)

    # Sort categories according to defined order
    sorted_categories = OrderedDict()
    for cat in CATEGORY_ORDER:
        if cat in by_category:
            sorted_categories[cat] = by_category[cat]

    # Add any remaining categories not in the predefined order
    for cat in by_category:
        if cat not in sorted_categories:
            sorted_categories[cat] = by_category[cat]

    # Generate markdown
    md_parts = [
        "# Full List of INPUT Keywords",
        "",
        "<!-- This file is auto-generated by docs/generate_docs_from_source.py -->",
        "<!-- Do not edit manually - changes will be overwritten -->",
        "",
        "<!-- Table of Contents -->",
        generate_toc(sorted_categories),
        ""
    ]

    for category, params in sorted_categories.items():
        if verbose:
            print(f"Category '{category}': {len(params)} parameters")
        md_parts.append(generate_category_markdown(category, params))

    # Write output
    output_content = '\n'.join(md_parts)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(output_content)

    print(f"Generated {output}")
    print(f"Categories: {len(sorted_categories)}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate markdown documentation from C++ source files'
    )
    parser.add_argument(
        '--source-dir',
        type=Path,
        default=Path('source/source_io'),
        help='Directory containing read_input_item_*.cpp files'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('docs/advanced/input_files/input-main.md'),
        help='Output markdown file'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print verbose output'
    )

    args = parser.parse_args()
    generate(args.source_dir, args.output, args.verbose)


if __name__ == '__main__':
    main()
