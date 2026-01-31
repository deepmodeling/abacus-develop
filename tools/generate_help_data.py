#!/usr/bin/env python3
"""
Generate C++ code with static parameter data from ABACUS INPUT documentation.

This script parses docs/advanced/input_files/input-main.md and generates
a C++ header file with static constexpr data structures for the help system.

Usage:
    python3 generate_help_data.py <input_markdown> <output_header>

Expected Input Format:
    The markdown file must follow this structure:

    1. Category Headers (optional):
       ## Category Name

       Defines a category for grouping parameters. All parameters following
       this header belong to this category until the next ## header.

    2. Parameter Headers (required):
       ### parameter_name

       Or for multi-parameter groups:
       ### param1, param2, param3

       Parameter names should be valid identifiers (alphanumeric + underscore).
       Multi-parameter headers will be split into individual parameter entries.

    3. Parameter Fields (at least Type and Description required):
       - **Type**: Integer|Real|String|Boolean
       - **Description**: Detailed description of the parameter...
         Descriptions can span multiple lines and include sub-bullets.
       - **Default**: default_value (optional)
       - **Unit**: unit_string (optional)
       - **Availability**: availability_conditions (optional)
       - **Note**: notes (ignored by parser)

    4. Format Assumptions:
       - Category headers start with '## ' (two hashes + space)
       - Parameter headers start with '### ' (three hashes + space)
       - Field markers use format: '- **FieldName**: value'
       - Descriptions should not contain lines starting with ## or ###
       - Multi-line descriptions are flattened into single lines
       - Markdown formatting (bold, italic, links, code) is stripped

    5. Example Structure:
       ## System Variables

       ### suffix

       - **Type**: String
       - **Description**: A suffix to the output directory name.
       - **Default**: ABACUS

       ### nx, ny, nz

       - **Type**: Integer
       - **Description**: FFT grid dimensions in x, y, z directions.
       - **Default**: 0
       - **Unit**: None

Design Philosophy:
    This parser is intentionally lightweight and assumes a controlled,
    well-formatted documentation file. It does not use full markdown parsing
    libraries to keep dependencies minimal. If the documentation format becomes
    significantly more complex, consider migrating to an AST-based parser.
"""

import re
import sys
from collections import Counter
from typing import List, Dict, Optional


def escape_cpp_string(text: str) -> str:
    """Escape special characters for C++ string literals."""
    # Replace backslash first to avoid double-escaping
    text = text.replace('\\', '\\\\')
    # Replace double quotes
    text = text.replace('"', '\\"')
    # Replace newlines with space (we'll use single-line strings)
    text = text.replace('\n', ' ')
    # Replace tabs with space
    text = text.replace('\t', ' ')
    # Collapse multiple spaces
    text = re.sub(r'\s+', ' ', text)
    # Strip leading/trailing whitespace
    text = text.strip()
    return text


def clean_markdown_formatting(text: str) -> str:
    """Remove markdown formatting from text."""
    # Remove LaTeX math mode expressions ($...$)
    text = re.sub(r'\$[^$]+\$', '', text)
    # Unescape markdown-escaped underscores (\_  -> _)
    text = re.sub(r'\\_', '_', text)
    # Remove bold/italic markers (use non-greedy matching for robustness)
    text = re.sub(r'\*\*(.*?)\*\*', r'\1', text)
    text = re.sub(r'\*([^*]+)\*', r'\1', text)
    # Remove inline code markers
    text = re.sub(r'`([^`]+)`', r'\1', text)
    # Remove markdown links [text](url) -> text
    text = re.sub(r'\[([^\]]+)\]\([^\)]+\)', r'\1', text)
    return text


def parse_parameter_section(lines: List[str], start_idx: int, category: str) -> Optional[List[Dict[str, str]]]:
    """
    Parse a single parameter section starting at the given line index.
    Returns a list of parameter dictionaries (one for each parameter name if multi-param),
    or None if parsing fails.
    """
    # Parameter name is in the ### header
    header_line = lines[start_idx].strip()
    if not header_line.startswith('###'):
        return None

    # Extract parameter name (remove ### and strip)
    param_name = header_line[3:].strip()

    # Handle parameters with multiple names separated by comma (e.g., "nx, ny, nz")
    # Split them into individual parameters
    if ',' in param_name:
        param_names = [n.strip() for n in param_name.split(',')]
        # Validate parameter names to catch markdown formatting issues
        for name in param_names:
            if re.search(r'[\s()]', name):
                print(f"  WARNING: Suspicious parameter name '{name}' contains spaces or parentheses")
                print(f"  This might indicate a markdown formatting issue in header: {param_name}")
        print(f"  Parsing multi-parameter: {param_name} -> {param_names}")
    else:
        param_names = [param_name]

    param = {
        'name': param_names[0] if len(param_names) == 1 else param_names,  # Store list for multi-params
        'category': category,
        'type': '',
        'description': '',
        'default': '',
        'unit': '',
        'availability': ''
    }

    # Parse bullet points after the header
    i = start_idx + 1
    current_field = None
    description_lines = []

    while i < len(lines):
        line = lines[i].strip()

        # Stop at next parameter (### header) or next category (## header)
        # Note: Assumes descriptions don't contain lines starting with ## or ###
        # This is a safe assumption for the controlled INPUT documentation format
        if line.startswith('###') or line.startswith('##'):
            break

        # Check for field markers
        if line.startswith('- **Type**:'):
            current_field = 'type'
            param['type'] = line.split(':', 1)[1].strip() if ':' in line else ''
        elif line.startswith('- **Description**:'):
            current_field = 'description'
            desc_text = line.split(':', 1)[1].strip() if ':' in line else ''
            if desc_text:
                description_lines.append(desc_text)
        elif line.startswith('- **Default**:'):
            current_field = 'default'
            param['default'] = line.split(':', 1)[1].strip() if ':' in line else ''
        elif line.startswith('- **Unit**:'):
            current_field = 'unit'
            param['unit'] = line.split(':', 1)[1].strip() if ':' in line else ''
        elif line.startswith('- **Availability**:'):
            current_field = 'availability'
            param['availability'] = line.split(':', 1)[1].strip() if ':' in line else ''
        elif line.startswith('- **Note**:'):
            # Ignore notes
            current_field = None
        elif line.startswith('-') and current_field == 'description':
            # Continuation of description (sub-bullets)
            # Remove leading "  - " or "- "
            desc_line = re.sub(r'^\s*-\s*', '', line)
            if desc_line:
                description_lines.append(desc_line)
        elif line and current_field == 'description' and not line.startswith('-'):
            # Multi-line description continuation
            description_lines.append(line)

        i += 1

    # Join description lines with spaces (intentional flattening)
    # The C++ side (input_help.cpp) has word-wrapping logic that works on single-line text
    # Preserving paragraph breaks would require more complex C++ layout handling
    if description_lines:
        param['description'] = ' '.join(description_lines)

    # Clean markdown formatting from all fields
    for key in param:
        if isinstance(param[key], str):
            param[key] = clean_markdown_formatting(param[key])

    # Validate required fields
    if not param['type']:
        name_str = ', '.join(param['name']) if isinstance(param['name'], list) else param['name']
        print(f"  WARNING: Parameter '{name_str}' missing Type field")
        return None

    if not param['description']:
        name_str = ', '.join(param['name']) if isinstance(param['name'], list) else param['name']
        print(f"  WARNING: Parameter '{name_str}' has empty description")

    # If multi-parameter, create separate entries for each name
    if isinstance(param['name'], list):
        params = []
        for name in param['name']:
            individual_param = param.copy()
            individual_param['name'] = name
            params.append(individual_param)
        return params
    else:
        return [param]


def parse_markdown(md_file: str) -> List[Dict[str, str]]:
    """Parse the markdown file and extract all parameter definitions."""
    try:
        with open(md_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"ERROR: File not found: {md_file}")
        sys.exit(1)
    except IOError as e:
        print(f"ERROR: Cannot read file {md_file}: {e}")
        sys.exit(1)

    parameters = []
    current_category = "General"

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Check for category header (## level)
        if line.startswith('## ') and not line.startswith('###'):
            # Extract category name
            category = line[3:].strip()
            # Remove markdown formatting
            category = clean_markdown_formatting(category)
            current_category = category
            print(f"Found category: {current_category}")

        # Check for parameter header (### level)
        elif line.startswith('###'):
            params = parse_parameter_section(lines, i, current_category)
            if params:
                for param in params:
                    parameters.append(param)
                    print(f"  Parsed parameter: {param['name']}")

        i += 1

    return parameters


def generate_cpp_header(parameters: List[Dict[str, str]], output_file: str):
    """Generate C++ header file with static const parameter data."""

    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            # Write header
            f.write("// Auto-generated by tools/generate_help_data.py - DO NOT EDIT\n")
            f.write("// Generated from: docs/advanced/input_files/input-main.md\n\n")
            f.write("#ifndef INPUT_HELP_DATA_H\n")
            f.write("#define INPUT_HELP_DATA_H\n\n")
            f.write("#include <cstddef>\n\n")
            f.write("namespace ModuleIO {\n\n")

            # Write struct definition
            f.write("struct ParameterInfo {\n")
            f.write("    const char* name;\n")
            f.write("    const char* category;\n")
            f.write("    const char* type;\n")
            f.write("    const char* description;\n")
            f.write("    const char* default_value;\n")
            f.write("    const char* unit;          // nullptr if no unit\n")
            f.write("    const char* availability;  // nullptr if always available\n")
            f.write("};\n\n")

            # Write parameter data array
            f.write("// Parameter data extracted from documentation\n")
            f.write("static constexpr ParameterInfo PARAMETER_DATA[] = {\n")

            for param in parameters:
                f.write("    {\n")
                f.write(f'        "{escape_cpp_string(param["name"])}",\n')
                f.write(f'        "{escape_cpp_string(param["category"])}",\n')
                f.write(f'        "{escape_cpp_string(param["type"])}",\n')
                f.write(f'        "{escape_cpp_string(param["description"])}",\n')
                f.write(f'        "{escape_cpp_string(param["default"])}",\n')
                # Use nullptr for empty optional fields
                if param['unit']:
                    f.write(f'        "{escape_cpp_string(param["unit"])}",\n')
                else:
                    f.write('        nullptr,\n')
                if param['availability']:
                    f.write(f'        "{escape_cpp_string(param["availability"])}"\n')
                else:
                    f.write('        nullptr\n')
                f.write("    },\n")

            f.write("};\n\n")

            # Write count constant
            f.write("static constexpr size_t PARAMETER_COUNT = "
                    f"sizeof(PARAMETER_DATA) / sizeof(ParameterInfo);\n\n")

            # Close namespace and header guard
            f.write("}  // namespace ModuleIO\n\n")
            f.write("#endif  // INPUT_HELP_DATA_H\n")
    except IOError as e:
        print(f"ERROR: Cannot write to file {output_file}: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 generate_help_data.py <input_markdown> <output_header>")
        sys.exit(1)

    input_md = sys.argv[1]
    output_h = sys.argv[2]

    print(f"Parsing markdown: {input_md}")
    parameters = parse_markdown(input_md)

    print(f"\nParsed {len(parameters)} parameters")

    # Check for duplicates using Counter (O(n) instead of O(n²))
    name_counts = Counter(p['name'] for p in parameters)
    duplicates = [name for name, count in name_counts.items() if count > 1]
    if duplicates:
        print(f"ERROR: Duplicate parameter names: {set(duplicates)}")
        sys.exit(1)

    # Count by category
    categories = {}
    for p in parameters:
        cat = p['category']
        categories[cat] = categories.get(cat, 0) + 1

    print("\nParameters by category:")
    for cat, count in sorted(categories.items()):
        print(f"  {cat}: {count}")

    print(f"\nGenerating C++ header: {output_h}")
    generate_cpp_header(parameters, output_h)

    print("Done!")


if __name__ == '__main__':
    main()
