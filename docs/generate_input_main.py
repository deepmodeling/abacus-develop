#!/usr/bin/env python3
"""
Generate input-main.md from a YAML parameter dump.

Usage:
    # Generate YAML first:
    abacus --generate-parameters-yaml > docs/parameters.yaml

    # Then generate markdown:
    python docs/generate_input_main.py docs/parameters.yaml [--output FILE]

Can also be imported from conf.py:
    from generate_input_main import generate
    generate(yaml_path, output)
"""

import argparse
import html
import re
import sys
from pathlib import Path
from collections import OrderedDict
from typing import Dict, List

try:
    import yaml
except ImportError:
    sys.exit("Error: PyYAML is required. Install with: pip install pyyaml")


DOC_FOLDER = Path(__file__).parent
# Desired category display order
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


def normalize_type(type_text: str) -> str:
    """
    Normalize legacy type labels for cleaner generated docs.
    """
    if not type_text:
        return ''

    normalized = type_text.strip()
    aliases = {
        "Bool": "Boolean",
        "Int*2": "Integer*2",
    }
    return aliases.get(normalized, normalized)


def escape_md_text(text: str) -> str:
    """
    Escape markdown-sensitive angle brackets while avoiding double escaping.
    """
    if text is None:
        return ''
    normalized = html.unescape(str(text))
    return normalized.replace('<', '&lt;').replace('>', '&gt;')


def format_description(desc: str) -> str:
    """
    Format description text for markdown output.
    - Convert * list markers to - list markers
    - Convert [NOTE] markers to blockquotes
    - Convert [WARNING] markers to blockquotes
    """
    if not desc:
        return ''

    # Prevent placeholder tokens like <property> from being parsed as raw HTML
    # and breaking list/heading structure in rendered docs.
    desc = escape_md_text(desc)

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

        # Normalize doubled note/warning prefixes from legacy content
        line = re.sub(r'>\s*Note:\s*Note\s*:?\s*', '> Note: ', line)
        line = re.sub(r'>\s*Warning:\s*Warning\s*:?\s*', '> Warning: ', line)

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
    if param.get('type', '') != '':
        type_text = escape_md_text(normalize_type(str(param['type'])))
        lines.append(f"- **Type**: {type_text}")

    # Availability (before description, as in original format)
    if param.get('availability', '') != '':
        availability_text = escape_md_text(str(param['availability']))
        lines.append(f"- **Availability**: *{availability_text}*")

    # Description
    if param.get('description', '') != '':
        desc = format_description(str(param['description']))
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
    if param.get('default_value', '') != '':
        default_text = escape_md_text(str(param['default_value']))
        lines.append(f"- **Default**: {default_text}")

    # Unit
    if param.get('unit', '') != '':
        unit_text = escape_md_text(str(param['unit']))
        lines.append(f"- **Unit**: {unit_text}")

    lines.append("")
    return '\n'.join(lines)


def generate_category_markdown(category: str, params: List[Dict[str, str]]) -> str:
    """
    Generate markdown for a category section.
    """
    lines = [f"## {category}", ""]

    for param in params:
        lines.append(generate_parameter_markdown(param))

    # Keep legacy navigation aid used by downstream tooling/rendering.
    lines.append("[back to top](#full-list-of-input-keywords)")
    lines.append("")

    return '\n'.join(lines)


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


def generate(yaml_path: Path, output: Path, verbose: bool = False):
    """
    Core generation logic. Can be called from conf.py or CLI.
    """
    yaml_path = Path(yaml_path)
    output = Path(output)

    if not yaml_path.exists():
        raise FileNotFoundError(f"YAML file not found: {yaml_path}")

    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)

    all_params = data.get('parameters', [])
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
        "<!-- This file is auto-generated from parameters.yaml -->",
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
        description='Generate input-main.md from YAML parameter dump'
    )
    parser.add_argument(
        'yaml_file',
        type=Path,
        help='Path to parameters.yaml (generated by abacus --generate-parameters-yaml)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path(f'{DOC_FOLDER}/advanced/input_files/input-main.md'),
        help='Output markdown file'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print verbose output'
    )

    args = parser.parse_args()
    generate(args.yaml_file, args.output, args.verbose)


if __name__ == '__main__':
    main()
