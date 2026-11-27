#!/usr/bin/env python3
import os
import re
import sys
from collections import defaultdict

INCLUDE_PATTERN = re.compile(r'^\s*#\s*include\s*["<]([^">]+)[">]', re.MULTILINE)

def find_source_files(root_dir):
    cpp_files = []
    h_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith('.cpp'):
                cpp_files.append(os.path.join(dirpath, filename))
            elif filename.endswith('.h'):
                h_files.append(os.path.join(dirpath, filename))
    return cpp_files + h_files

def parse_includes(filepath, include_dirs):
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except:
        return []
    
    includes = []
    for match in INCLUDE_PATTERN.finditer(content):
        include_path = match.group(1)
        includes.append(include_path)
    return includes

def resolve_include(include_name, source_file_dir, include_dirs):
    for include_dir in [source_file_dir] + include_dirs:
        full_path = os.path.join(include_dir, include_name)
        if os.path.exists(full_path):
            return os.path.normpath(full_path)
    return None

def compute_depth(filepath, include_map, visited, depth_cache, include_dirs):
    if filepath in depth_cache:
        return depth_cache[filepath]
    if filepath in visited:
        return 0
    
    visited.add(filepath)
    max_depth = 0
    source_dir = os.path.dirname(filepath)
    
    if filepath not in include_map:
        depth_cache[filepath] = 0
        visited.remove(filepath)
        return 0
    
    for include_name in include_map[filepath]:
        resolved_path = resolve_include(include_name, source_dir, include_dirs)
        if resolved_path:
            child_depth = compute_depth(resolved_path, include_map, visited, depth_cache, include_dirs)
            max_depth = max(max_depth, child_depth + 1)
    
    depth_cache[filepath] = max_depth
    visited.remove(filepath)
    return max_depth

def get_include_chain(filepath, include_map, include_dirs, max_depth=30):
    chain = [filepath]
    visited = set()
    visited.add(filepath)
    current = filepath
    
    for _ in range(max_depth):
        source_dir = os.path.dirname(current)
        if current not in include_map or not include_map[current]:
            break
        
        next_include = None
        max_child_depth = -1
        
        for include_name in include_map[current]:
            resolved_path = resolve_include(include_name, source_dir, include_dirs)
            if resolved_path and resolved_path not in visited:
                child_depth = 0
                temp_visited = set()
                child_depth = compute_depth(resolved_path, include_map, temp_visited, {}, include_dirs)
                if child_depth > max_child_depth:
                    max_child_depth = child_depth
                    next_include = resolved_path
        
        if next_include:
            visited.add(next_include)
            chain.append(next_include)
            current = next_include
        else:
            break
    
    return chain

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Script is at tools/03_code_analysis/, repo root is two levels up
    repo_root = os.path.dirname(os.path.dirname(script_dir))
    root_dir = os.path.join(repo_root, 'source')
    
    include_dirs = [
        root_dir,
        os.path.join(root_dir, 'source_base'),
        os.path.join(root_dir, 'source_cell'),
        os.path.join(root_dir, 'source_pw'),
        os.path.join(root_dir, 'source_lcao'),
        os.path.join(root_dir, 'source_basis'),
        os.path.join(root_dir, 'source_esolver'),
        os.path.join(root_dir, 'source_hsolver'),
        os.path.join(root_dir, 'source_io'),
        os.path.join(root_dir, 'source_psi'),
        os.path.join(root_dir, 'source_relax'),
        os.path.join(root_dir, 'source_estate'),
        os.path.join(root_dir, 'source_hamilt'),
        os.path.join(root_dir, 'source_main'),
        os.path.join(root_dir, 'source_md'),
    ]
    
    all_files = find_source_files(root_dir)
    print(f"Found {len(all_files)} source files")
    
    include_map = {}
    for filepath in all_files:
        includes = parse_includes(filepath, include_dirs)
        if includes:
            include_map[filepath] = includes
    
    print(f"Parsed includes for {len(include_map)} files")
    
    depth_cache = {}
    file_depths = []
    
    for filepath in all_files:
        visited = set()
        depth = compute_depth(filepath, include_map, visited, depth_cache, include_dirs)
        file_depths.append((filepath, depth))
    
    file_depths.sort(key=lambda x: -x[1])
    
    min_depth = 18
    deep_files = [(f, d) for f, d in file_depths if d >= min_depth]
    
    output = []
    output.append(f"=== Files with #include dependency depth >= {min_depth} ===")
    output.append(f"Total found: {len(deep_files)} files")
    output.append("")
    
    top_10 = deep_files[:10]
    if top_10:
        output.append("=" * 80)
        output.append("TOP 10 deepest dependency chains (showing full include path):")
        output.append("=" * 80)
        output.append("")
        
        for filepath, depth in top_10:
            chain = get_include_chain(filepath, include_map, include_dirs)
            rel_path = os.path.relpath(filepath, root_dir)
            output.append(f"  Depth {depth}: {rel_path}")
            output.append("  Include chain:")
            for i, step in enumerate(chain):
                prefix = "         " if i == 0 else "        -> "
                output.append(f"{prefix}{os.path.relpath(step, root_dir)}")
            output.append("")
    
    output.append("=" * 80)
    output.append(f"Full list of all {len(deep_files)} files with depth >= {min_depth}:")
    output.append("=" * 80)
    output.append("")
    
    depth_groups = defaultdict(list)
    for filepath, depth in deep_files:
        depth_groups[depth].append(filepath)
    
    for depth in sorted(depth_groups.keys(), reverse=True):
        files = sorted(depth_groups[depth])
        output.append(f"  Depth {depth} ({len(files)} files):")
        for filepath in files:
            rel_path = os.path.relpath(filepath, root_dir)
            output.append(f"    {rel_path}")
        output.append("")
    
    output_path = os.path.join(script_dir, 'deep_include_analysis.txt')
    with open(output_path, 'w') as f:
        f.write('\n'.join(output))
    
    print(f"Analysis complete.")
    print(f"Results have been generated at:")
    print(f"  {output_path}")
    print(f"Total files with depth >= {min_depth}: {len(deep_files)}")

if __name__ == '__main__':
    main()
