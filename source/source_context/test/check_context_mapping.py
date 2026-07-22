#!/usr/bin/env python3
"""Verify that the legacy-state manifests cover every declared data member."""

from __future__ import print_function

import argparse
import collections
import hashlib
import re
import sys
from pathlib import Path


def strip_comments(text):
    result = []
    index = 0
    state = "code"
    quote = ""
    while index < len(text):
        char = text[index]
        following = text[index + 1] if index + 1 < len(text) else ""
        if state == "line_comment":
            if char == "\n":
                result.append(char)
                state = "code"
            index += 1
            continue
        if state == "block_comment":
            if char == "*" and following == "/":
                state = "code"
                index += 2
            else:
                if char == "\n":
                    result.append(char)
                index += 1
            continue
        if state == "string":
            result.append(char)
            if char == "\\":
                if following:
                    result.append(following)
                    index += 2
                    continue
            elif char == quote:
                state = "code"
            index += 1
            continue
        if char == "/" and following == "/":
            state = "line_comment"
            index += 2
        elif char == "/" and following == "*":
            state = "block_comment"
            index += 2
        else:
            result.append(char)
            if char in ('"', "'"):
                state = "string"
                quote = char
            index += 1
    return "".join(result)


def block_body(text, keyword, name):
    clean = strip_comments(text)
    match = re.search(r"\b" + re.escape(keyword) + r"\s+" + re.escape(name) + r"\b[^\{]*\{", clean)
    if not match:
        raise ValueError("cannot find {} {}".format(keyword, name))
    start = match.end()
    depth = 1
    index = start
    quote = ""
    while index < len(clean):
        char = clean[index]
        if quote:
            if char == "\\":
                index += 2
                continue
            if char == quote:
                quote = ""
        elif char in ('"', "'"):
            quote = char
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return clean[start:index]
        index += 1
    raise ValueError("unterminated {} {}".format(keyword, name))


def top_level_parts(text, separator):
    parts = []
    start = 0
    braces = brackets = parentheses = angles = 0
    quote = ""
    index = 0
    while index < len(text):
        char = text[index]
        if quote:
            if char == "\\":
                index += 2
                continue
            if char == quote:
                quote = ""
        elif char in ('"', "'"):
            quote = char
        elif char == "{":
            braces += 1
        elif char == "}":
            braces -= 1
        elif char == "[":
            brackets += 1
        elif char == "]":
            brackets -= 1
        elif char == "(":
            parentheses += 1
        elif char == ")":
            parentheses -= 1
        elif char == "<":
            angles += 1
        elif char == ">" and angles:
            angles -= 1
        elif char == separator and not any((braces, brackets, parentheses, angles)):
            parts.append(text[start:index])
            start = index + 1
        index += 1
    parts.append(text[start:])
    return parts


def normalize_type(type_name):
    normalized = re.sub(r"\s+", " ", type_name.strip())
    normalized = re.sub(r"\s*([<>,*&:\[\]])\s*", r"\1", normalized)
    return re.sub(r"^(?:extern|mutable)\s+", "", normalized)


def data_member_schema(text, kind, name):
    body = block_body(text, kind, name)
    body = "\n".join(line for line in body.splitlines() if not line.lstrip().startswith("#"))
    members = []
    for statement in top_level_parts(body, ";"):
        statement = statement.strip()
        if not statement or statement.endswith(":"):
            continue
        if statement.startswith(("using ", "typedef ", "static_assert", "struct ", "class ")):
            continue
        if "(" in top_level_parts(statement, "=")[0]:
            continue
        declarators = top_level_parts(statement, ",")
        base_type = ""
        for index, declarator in enumerate(declarators):
            declaration = top_level_parts(declarator, "=")[0].strip()
            match = re.search(r"([A-Za-z_]\w*)\s*(?:\[[^\]]*\]\s*)*$", declaration)
            if match:
                prefix = declaration[: match.start()].strip()
                if index == 0:
                    base_type = prefix
                    member_type = base_type
                else:
                    member_type = base_type + (" " + prefix if prefix else "")
                members.append((match.group(1), normalize_type(member_type)))
    return members


def data_members(text, kind, name):
    return [field for field, _ in data_member_schema(text, kind, name)]


def manifest_fields(text, macro, path_field=False):
    if path_field:
        pattern = re.compile(r"\b" + re.escape(macro) + r"\(\s*([^\)]+?)\s*\)")
        return [match.group(1).strip() for match in pattern.finditer(text)]
    pattern = re.compile(
        r"\b" + re.escape(macro) + r"\(\s*([A-Za-z_]\w*)\s*,\s*([A-Za-z_]\w*)\s*\)"
    )
    return [match.group(2) for match in pattern.finditer(text)]


def compare(label, declared, mapped):
    errors = []
    declared_counts = collections.Counter(declared)
    mapped_counts = collections.Counter(mapped)
    duplicate_declarations = sorted(name for name, count in declared_counts.items() if count != 1)
    duplicate_mappings = sorted(name for name, count in mapped_counts.items() if count != 1)
    missing = sorted(set(declared_counts) - set(mapped_counts))
    extra = sorted(set(mapped_counts) - set(declared_counts))
    if duplicate_declarations:
        errors.append("{} duplicate declarations: {}".format(label, ", ".join(duplicate_declarations)))
    if duplicate_mappings:
        errors.append("{} duplicate mappings: {}".format(label, ", ".join(duplicate_mappings)))
    if missing:
        errors.append("{} fields missing from manifest: {}".format(label, ", ".join(missing)))
    if extra:
        errors.append("{} manifest fields not declared by legacy type: {}".format(label, ", ".join(extra)))
    return errors


def schema_digest(schema):
    canonical = "".join("{}:{}\n".format(field, type_name) for field, type_name in sorted(schema))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def manifest_schema_digest(text):
    match = re.search(r"MODULE_CONTEXT_LEGACY_SCHEMA_SHA256\(([0-9a-f]{64})\)", text)
    return match.group(1) if match else None


def compare_schema_digest(label, schema, manifest):
    expected = manifest_schema_digest(manifest)
    actual = schema_digest(schema)
    if expected == actual:
        return []
    if expected is None:
        return ["{} manifest has no legacy schema type fingerprint (actual {})".format(label, actual)]
    return [
        "{} legacy field types changed: manifest fingerprint {}, actual {}".format(label, expected, actual)
    ]


def read(root, relative):
    return (root / relative).read_text(encoding="utf-8")


def collect_schemas(root):
    context = root / "source/source_context"
    input_header = read(root, "source/source_io/module_parameter/input_parameter.h")
    input_manifest = read(context, "input_field_mapping.inc")
    input_schema = data_member_schema(input_header, "struct", "Input_para")

    system_header = read(root, "source/source_io/module_parameter/system_parameter.h")
    system_manifest = read(context, "system_field_mapping.inc")
    system_schema = data_member_schema(system_header, "struct", "System_para")

    global_header = read(root, "source/source_base/global_variable.h")
    global_manifest = read(context, "globalv_field_mapping.inc")
    global_schema = data_member_schema(global_header, "namespace", "GlobalV")

    exx_sections = {
        "info_global": ("source/source_hamilt/module_xc/exx_info_global.h", "Exx_Info_Global"),
        "info_lip": ("source/source_hamilt/module_xc/exx_info_lip.h", "Exx_Info_Lip"),
        "info_ri": ("source/source_hamilt/module_xc/exx_info_ri.h", "Exx_Info_RI"),
        "info_opt_abfs": ("source/source_hamilt/module_xc/exx_info_opt_abfs.h", "Exx_Info_Opt_ABFs"),
    }
    exx_schema = []
    for section, (header_path, type_name) in exx_sections.items():
        exx_schema.extend(
            (section + "." + field, type_name_value)
            for field, type_name_value in data_member_schema(read(root, header_path), "struct", type_name)
        )
    exx_manifest = read(context, "exx_state_field_mapping.inc")

    restart_header = read(root, "source/source_io/module_restart/restart.h")
    restart_body = block_body(restart_header, "class", "Restart")
    save_member = re.search(r"\bInfo_Save\s+([A-Za-z_]\w*)\s*;", strip_comments(restart_body))
    load_member = re.search(r"\bInfo_Load\s+([A-Za-z_]\w*)\s*;", strip_comments(restart_body))
    restart_schema = []
    if save_member and load_member:
        restart_schema.extend(
            (save_member.group(1) + "." + field, type_name_value)
            for field, type_name_value in data_member_schema(restart_header, "struct", "Info_Save")
        )
        restart_schema.extend(
            (load_member.group(1) + "." + field, type_name_value)
            for field, type_name_value in data_member_schema(restart_header, "struct", "Info_Load")
        )
        restart_schema.append(("folder", "std::string"))
    restart_manifest = read(context, "restart_field_mapping.inc")

    return {
        "Input_para": (input_schema, input_manifest),
        "System_para": (system_schema, system_manifest),
        "GlobalV": (global_schema, global_manifest),
        "GlobalC::exx_info": (exx_schema, exx_manifest),
        "GlobalC::restart": (restart_schema, restart_manifest),
    }, save_member, load_member


def check(root):
    errors = []
    schemas, save_member, load_member = collect_schemas(root)

    input_schema, input_manifest = schemas["Input_para"]
    errors.extend(
        compare(
            "Input_para",
            [field for field, _ in input_schema],
            manifest_fields(input_manifest, "MODULE_CONTEXT_INPUT_FIELD"),
        )
    )
    errors.extend(compare_schema_digest("Input_para", input_schema, input_manifest))

    system_schema, system_manifest = schemas["System_para"]
    errors.extend(
        compare(
            "System_para",
            [field for field, _ in system_schema],
            manifest_fields(system_manifest, "MODULE_CONTEXT_SYSTEM_FIELD"),
        )
    )
    errors.extend(compare_schema_digest("System_para", system_schema, system_manifest))

    global_schema, global_manifest = schemas["GlobalV"]
    errors.extend(
        compare(
            "GlobalV",
            [field for field, _ in global_schema],
            manifest_fields(global_manifest, "MODULE_CONTEXT_GLOBALV_FIELD"),
        )
    )
    errors.extend(compare_schema_digest("GlobalV", global_schema, global_manifest))

    exx_schema, exx_manifest = schemas["GlobalC::exx_info"]
    declared_exx = [field for field, _ in exx_schema]
    mapped_exx = []
    for match in re.finditer(
        r"\bMODULE_CONTEXT_EXX_STATE_FIELD\(\s*([A-Za-z_]\w*)\s*,\s*([A-Za-z_]\w*)\s*\)",
        exx_manifest,
    ):
        mapped_exx.append(match.group(1) + "." + match.group(2))
    errors.extend(compare("GlobalC::exx_info", declared_exx, mapped_exx))
    errors.extend(compare_schema_digest("GlobalC::exx_info", exx_schema, exx_manifest))

    restart_schema, restart_manifest = schemas["GlobalC::restart"]
    if not save_member or not load_member:
        errors.append("Restart state container declarations could not be identified")
    else:
        declared_restart = [field for field, _ in restart_schema]
        errors.extend(
            compare(
                "GlobalC::restart",
                declared_restart,
                manifest_fields(restart_manifest, "MODULE_CONTEXT_RESTART_FIELD", path_field=True),
            )
        )
        errors.extend(compare_schema_digest("GlobalC::restart", restart_schema, restart_manifest))
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--print-digests", action="store_true")
    args = parser.parse_args()
    if args.print_digests:
        schemas, _, _ = collect_schemas(args.repo_root.resolve())
        for label, (schema, _) in sorted(schemas.items()):
            print("{} {}".format(label, schema_digest(schema)))
        return 0
    errors = check(args.repo_root.resolve())
    if errors:
        for error in errors:
            print("context mapping error: " + error, file=sys.stderr)
        return 1
    print("Context mapping manifests cover all legacy fields exactly once.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
