"""Per-file analysis orchestration for the ABACUS code-quality scorer."""

import re
from pathlib import Path
from typing import List, Sequence, Set, Tuple

from cqs_config import (
    AUTO_RE, CAPS, CHINESE_RE, COMMENT_LINE_RE, CYCLO_THRESHOLD,
    DEFAULT_PARAM_RE, DELETE_EXPR_RE, DOC_FILE_RE, FILE_LENGTH_STEP,
    FILE_LENGTH_THRESHOLD, FILENAME_LENGTH_LIMIT, FRIEND_RE,
    FUNCTION_PARAM_THRESHOLD, GLOBAL_DEPENDENCY_RE, GOTO_RE, GUARD_RE,
    HEADER_FUNC_BODY_STEP, HEADER_FUNC_BODY_THRESHOLD,
    HEADER_LENGTH_STEP, HEADER_LENGTH_THRESHOLD, HPP_INCLUDE_RE,
    INDENTED_PREPROC_RE, LINE_LENGTH_LIMIT, MACRO_DEFINE_RE, NEW_EXPR_RE,
    OWNED_NEW_RE, SKIP_DIRS, SOURCE_EXTENSIONS,
    UPPERCASE_CONST_RE,
    USING_NS_STD_RE, WEIGHTS,
)
from cqs_report import FileReport, Finding
from cqs_rules import (
    analyze_class_blocks, find_complex_call_arguments,
    find_function_bodies, find_header_function_bodies,
    find_high_complexity_functions, find_ifdef_in_call_args,
    find_long_function_signatures, find_param_member_inits,
    find_post_cpp11_features,
)
from cqs_text import extract_comment_blocks, strip_comments, strip_strings


def discover_files(paths: Sequence[str]) -> List[Path]:
    """Walk input paths and return source files, skipping generated dirs."""
    result: List[Path] = []
    for p in paths:
        path = Path(p)
        if path.is_file():
            if path.suffix in SOURCE_EXTENSIONS:
                result.append(path)
        elif path.is_dir():
            for sub in path.rglob("*"):
                if not sub.is_file():
                    continue
                if any(part in SKIP_DIRS for part in sub.parts):
                    continue
                if sub.suffix in SOURCE_EXTENSIONS:
                    result.append(sub)
    return result


def analyze_file(path: Path) -> FileReport:
    """Analyze one file and return its FileReport."""
    findings: List[Finding] = []
    name = path.name
    stem = path.stem
    suffix = path.suffix

    # filename rules
    if len(stem) > FILENAME_LENGTH_LIMIT:
        findings.append(Finding(
            rule="filename_too_long",
            line=None,
            reason=(
                f"filename '{name}' stem has {len(stem)} chars "
                f"(limit {FILENAME_LENGTH_LIMIT})"
            ),
            deduction=WEIGHTS["filename_too_long"],
        ))
    if any(c.isupper() for c in stem):
        findings.append(Finding(
            rule="filename_uppercase",
            line=None,
            reason=f"filename '{name}' contains uppercase letters",
            deduction=WEIGHTS["filename_uppercase"],
        ))
    if suffix == ".hpp":
        findings.append(Finding(
            rule="hpp_implementation",
            line=None,
            reason=".hpp implementation header prohibited (use .cpp + .h split)",
            deduction=WEIGHTS["hpp_implementation"],
        ))

    # read content
    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        findings.append(Finding(
            rule="read_error",
            line=None,
            reason=f"failed to read file: {e}",
            deduction=0,
        ))
        return FileReport(path=str(path), score=100, findings=findings)

    lines = content.split("\n")

    # file length rule (no cap: large files legitimately scale deduction)
    total_lines = len(lines)
    if total_lines > FILE_LENGTH_THRESHOLD:
        excess = total_lines - FILE_LENGTH_THRESHOLD
        blocks = (excess + FILE_LENGTH_STEP - 1) // FILE_LENGTH_STEP
        findings.append(Finding(
            rule="file_too_long",
            line=None,
            reason=(
                f"file has {total_lines} lines "
                f"(exceeds {FILE_LENGTH_THRESHOLD} by {excess})"
            ),
            deduction=blocks * WEIGHTS["file_too_long"],
        ))

    # .h header length rule (cap 5: headers should stay compact)
    if path.suffix == ".h" and total_lines > HEADER_LENGTH_THRESHOLD:
        header_excess = total_lines - HEADER_LENGTH_THRESHOLD
        header_blocks = (
            (header_excess + HEADER_LENGTH_STEP - 1) // HEADER_LENGTH_STEP
        )
        header_deduction = min(
            header_blocks * WEIGHTS["header_too_long"],
            CAPS["header_too_long"],
        )
        findings.append(Finding(
            rule="header_too_long",
            line=None,
            reason=(
                f".h header has {total_lines} lines "
                f"(exceeds {HEADER_LENGTH_THRESHOLD} by {header_excess})"
            ),
            deduction=header_deduction,
        ))

    # non-inline, non-template function bodies in .h headers
    if path.suffix == ".h":
        for fn_line, fn_name, fn_lines in find_header_function_bodies(
            content
        ):
            excess = fn_lines - HEADER_FUNC_BODY_THRESHOLD
            blocks = (
                (excess + HEADER_FUNC_BODY_STEP - 1)
                // HEADER_FUNC_BODY_STEP
            )
            deduction = min(
                blocks * WEIGHTS["header_function_body"],
                CAPS["header_function_body"],
            )
            findings.append(Finding(
                rule="header_function_body",
                line=fn_line,
                reason=(
                    f"non-inline function `{fn_name}` defined in .h header "
                    f"spans {fn_lines} lines "
                    f"(exceeds {HEADER_FUNC_BODY_THRESHOLD}; move to .cpp)"
                ),
                deduction=deduction,
            ))

    # line-based rules with caps
    tab_count = sum(1 for l in lines if "\t" in l)
    long_line_count = sum(1 for l in lines if len(l) > LINE_LENGTH_LIMIT)
    using_ns_count = sum(
        1 for l in lines if USING_NS_STD_RE.search(strip_comments(l))
    )
    chinese_count = sum(1 for l in lines if CHINESE_RE.search(l))

    stripped_content = strip_comments(content)
    # uppercase constants: strip strings too (so default-value strings like
    # "ABACUS" do not trigger), and the regex excludes member accesses
    # (e.g. GlobalV::MY_RANK, x.UPPER) via a negative lookbehind.
    stripped_for_upper = strip_strings(stripped_content)
    upper_const_count = sum(
        len(UPPERCASE_CONST_RE.findall(l))
        for l in stripped_for_upper.split("\n")
    )

    # interface & dependency rules
    global_dep_count = len(GLOBAL_DEPENDENCY_RE.findall(stripped_content))
    hpp_include_count = sum(1 for l in lines if HPP_INCLUDE_RE.match(l))
    friend_count = len(FRIEND_RE.findall(stripped_content))
    goto_count = len(GOTO_RE.findall(stripped_content))
    # `auto` keyword: counted on comment-and-string-stripped content so
    # occurrences inside comments or string literals do not deduct points.
    auto_count = len(AUTO_RE.findall(stripped_for_upper))

    # default parameter: scan the full stripped content so multi-line
    # declarations (parameter list spanning multiple lines) are detected;
    # `[^();]*` in the regex already matches newlines.
    default_param_count = len(DEFAULT_PARAM_RE.findall(stripped_content))

    # unpaired new/delete (file-level): rough heuristic, cap applied below.
    # All `new` expressions are counted first. `delete` expressions are the
    # obvious symmetric ownership cancel, but there is a second class of
    # cancellations that never produces a `delete` keyword at all: a
    # C++11 (or later) smart pointer that owns the object from construction.
    # Recognised ownership-transfer patterns are aggregated in OWNED_NEW_RE
    # and subtracted before comparing with `delete` so the leak heuristic
    # stays meaningful on modern code that prefers std::unique_ptr /
    # std::shared_ptr to manual delete. Examples that now cancel:
    #     ptr.reset(new Foo());
    #     std::unique_ptr<Foo> p(new Foo(arg));
    # These cancellations deliberately do NOT touch `raw_new_count` below
    # because the "raw_new_keyword" rule charges per direct `new` regardless
    # of ownership — it is the stylistic counterpart to std::make_unique /
    # std::make_shared, so `reset(new Foo)` and `unique_ptr<T>(new T)`
    # must still contribute to that count.
    new_count = len(NEW_EXPR_RE.findall(stripped_content))
    delete_count = len(DELETE_EXPR_RE.findall(stripped_content))
    owned_new_count = len(OWNED_NEW_RE.findall(stripped_content))
    cancelled_new = delete_count + owned_new_count
    unpaired_new = max(0, new_count - cancelled_new)

    # raw `new` keyword usage: each occurrence costs 1 (no cap)
    raw_new_count = new_count

    def append_capped(rule: str, count: int) -> None:
        if count == 0:
            return
        cap = CAPS.get(rule)
        actual = min(count, cap) if cap else count
        cap_text = f" (capped at {cap})" if cap and count > cap else ""
        findings.append(Finding(
            rule=rule,
            line=None,
            reason=f"{count} occurrence(s){cap_text}",
            deduction=actual * WEIGHTS[rule],
        ))

    append_capped("tab_indentation", tab_count)
    append_capped("line_too_long", long_line_count)
    append_capped("using_namespace_std", using_ns_count)
    append_capped("chinese_comment", chinese_count)
    append_capped("uppercase_constant", upper_const_count)
    append_capped("default_parameter", default_param_count)
    append_capped("global_dependency", global_dep_count)
    append_capped("hpp_include", hpp_include_count)
    append_capped("friend_keyword", friend_count)
    append_capped("unpaired_new_delete", unpaired_new)
    append_capped("raw_new_keyword", raw_new_count)
    append_capped("goto_keyword", goto_count)
    append_capped("auto_keyword", auto_count)

    # --- documentation style rules -------------------------------------
    # indented preprocessor directive (raw lines, `#` not at column 0)
    indented_preproc_count = sum(1 for l in lines if INDENTED_PREPROC_RE.match(l))
    append_capped("indented_preprocessor", indented_preproc_count)

    # complex expressions (?:, &&, ||) written inline in a call's
    # argument list: one finding per offending call, cap per file
    complex_arg_hits = find_complex_call_arguments(content)
    cap_complex = CAPS.get("complex_argument_expression", len(complex_arg_hits))
    for hit_idx, (line_no, callee) in enumerate(complex_arg_hits):
        if hit_idx >= cap_complex:
            break
        findings.append(Finding(
            rule="complex_argument_expression",
            line=line_no,
            reason=(
                f"call to '{callee}' embeds a ternary/&&/|| expression "
                f"in its argument list (assign to a named local first)"
            ),
            deduction=WEIGHTS["complex_argument_expression"],
        ))

    # #if/#ifdef/#ifndef block wrapping individual call arguments
    ifdef_arg_lines = find_ifdef_in_call_args(content)
    cap_ifdef = CAPS.get("ifdef_in_call_args", len(ifdef_arg_lines))
    for hit_idx, line_no in enumerate(ifdef_arg_lines):
        if hit_idx >= cap_ifdef:
            break
        findings.append(Finding(
            rule="ifdef_in_call_args",
            line=line_no,
            reason=(
                "#if/#ifdef block inside a function call's argument list "
                "(select the argument before the call instead)"
            ),
            deduction=WEIGHTS["ifdef_in_call_args"],
        ))

    # new `#define` macros in .h/.cpp (guard-style macros exempt)
    if path.suffix in (".h", ".cpp"):
        macro_hits = 0
        cap_macro = CAPS.get("new_macro_definition", 0)
        for idx, l in enumerate(lines):
            if macro_hits >= cap_macro:
                break
            m = MACRO_DEFINE_RE.match(l)
            if not m:
                continue
            name = m.group(1)
            rest = l[m.end():].strip()
            # include guard: NAME ends in _H/_HPP and carries no
            # parameter list or replacement text
            if (name.endswith("_H") or name.endswith("_HPP")) and not rest:
                continue
            macro_hits += 1
            findings.append(Finding(
                rule="new_macro_definition",
                line=idx + 1,
                reason=(
                    f"new `#define {name}` macro requires approval "
                    "(AGENTS.md rule 20)"
                ),
                deduction=WEIGHTS["new_macro_definition"],
            ))

    # member variables initialized from PARAM in .h/.hpp (rule 21)
    if path.suffix in (".h", ".hpp"):
        param_init_lines = find_param_member_inits(content)
        cap_param_init = CAPS.get("param_member_init", 0)
        for hit_idx, line_no in enumerate(param_init_lines):
            if hit_idx >= cap_param_init:
                break
            findings.append(Finding(
                rule="param_member_init",
                line=line_no,
                reason=(
                    "member variable initialized from PARAM in a header "
                    "(pass the value through the constructor instead; "
                    "AGENTS.md rule 21)"
                ),
                deduction=WEIGHTS["param_member_init"],
            ))

    # comment block(s) above the first include guard (headers only). One
    # contiguous comment block counts once regardless of its length.
    if path.suffix == ".h":
        comment_above_guard = 0
        in_block = False
        for l in lines:
            if GUARD_RE.match(l):
                break
            if COMMENT_LINE_RE.match(l):
                if not in_block:
                    comment_above_guard += 1
                    in_block = True
            elif l.strip():
                in_block = False
        append_capped("comment_above_guard", comment_above_guard)

    # `@file xxx` name vs the actual filename
    doc_file_match = DOC_FILE_RE.search(content)
    if doc_file_match and doc_file_match.group(1) != path.name:
        findings.append(Finding(
            rule="doc_file_mismatch",
            line=None,
            reason=(
                f"@file names '{doc_file_match.group(1)}' but the file is "
                f"'{path.name}'"
            ),
            deduction=WEIGHTS["doc_file_mismatch"],
        ))

    # duplicated comment text across comment blocks (e.g. the same formula
    # repeated in the file-header and in a function doxygen block). Only
    # lines of >= 40 chars are compared. Excluded from comparison:
    #   - doxygen command lines (@param/@return/...): repeated `@param pv`
    #     across overloads is normal, not redundant prose;
    #   - decorative separator lines (=== Part N ===): >50% non-alnum chars;
    #   - ALL-CAPS banner lines: visual section headers, not content.
    def _is_dup_candidate(l: str) -> bool:
        if len(l) < 40 or re.match(r"@\w+", l):
            return False
        alnum = sum(1 for ch in l if ch.isalnum())
        if alnum * 2 < len(l):
            return False
        if l.isupper():
            return False
        return True

    comment_blocks = extract_comment_blocks(content)
    if len(comment_blocks) > 1:
        seen: dict = {}
        duplicate_doc_lines = 0
        for block in comment_blocks:
            block_lines = set(l for l in block if _is_dup_candidate(l))
            for l in block_lines:
                if l in seen:
                    duplicate_doc_lines += 1
                else:
                    seen[l] = True
        append_capped("duplicate_doc_block", duplicate_doc_lines)

    # too-many-parameters rule: each function's deduction is capped
    # individually (per-function cap, not per-file). A file with many
    # over-parameterised functions accumulates deductions across all of
    # them, instead of one bad function hiding the rest.
    long_param_funcs = find_long_function_signatures(content, FUNCTION_PARAM_THRESHOLD)
    cap_params = CAPS.get("too_many_parameters")
    for line_no, fname, pcount in long_param_funcs:
        excess = pcount - FUNCTION_PARAM_THRESHOLD
        raw_deduction = excess * WEIGHTS["too_many_parameters"]
        if cap_params is not None and raw_deduction > cap_params:
            per_deduction = cap_params
            cap_text = f" (capped at {cap_params})"
        else:
            per_deduction = raw_deduction
            cap_text = ""
        findings.append(Finding(
            rule="too_many_parameters",
            line=line_no,
            reason=(
                f"function '{fname}' has {pcount} parameters "
                f"(exceeds {FUNCTION_PARAM_THRESHOLD} by {excess}){cap_text}"
            ),
            deduction=per_deduction,
        ))

    # high cyclomatic complexity rule: each function's deduction is capped
    # individually (per-function cap, not per-file).
    high_cyclo_funcs = find_high_complexity_functions(content, CYCLO_THRESHOLD)
    cap_cyclo = CAPS.get("high_cyclomatic_complexity")
    for line_no, fname, cyclo in high_cyclo_funcs:
        excess = cyclo - CYCLO_THRESHOLD
        raw_deduction = excess * WEIGHTS["high_cyclomatic_complexity"]
        if cap_cyclo is not None and raw_deduction > cap_cyclo:
            per_deduction = cap_cyclo
            cap_text = f" (capped at {cap_cyclo})"
        else:
            per_deduction = raw_deduction
            cap_text = ""
        findings.append(Finding(
            rule="high_cyclomatic_complexity",
            line=line_no,
            reason=(
                f"function '{fname}' has cyclomatic complexity {cyclo} "
                f"(exceeds {CYCLO_THRESHOLD} by {excess}){cap_text}"
            ),
            deduction=per_deduction,
        ))

    # post-C++11 feature rule: base deduction for any usage + per-feature
    # incremental deduction (capped). Total cap matches the previous
    # one-shot -80, but a single misuse no longer swamps the whole file.
    post_cpp11_hits = find_post_cpp11_features(content)
    if post_cpp11_hits:
        feature_list = ", ".join(
            f"'{label}' (line {ln})" for ln, label in post_cpp11_hits
        )
        findings.append(Finding(
            rule="post_cpp11_feature",
            line=post_cpp11_hits[0][0],
            reason=(
                f"uses post-C++11 feature(s): {feature_list} "
                f"(repo baseline is C++11)"
            ),
            deduction=WEIGHTS["post_cpp11_feature"],
        ))
        num_features = len(post_cpp11_hits)
        cap_features = CAPS.get("post_cpp11_per_feature", num_features)
        capped = min(num_features, cap_features)
        if capped > 0:
            cap_text = (
                f" (capped at {cap_features})" if num_features > cap_features else ""
            )
            findings.append(Finding(
                rule="post_cpp11_per_feature",
                line=post_cpp11_hits[0][0],
                reason=(
                    f"{num_features} distinct post-C++11 feature(s) detected"
                    f"{cap_text}"
                ),
                deduction=capped * WEIGHTS["post_cpp11_per_feature"],
            ))

    # class-based rules
    (
        pub_findings,
        static_findings,
        long_func_findings,
        conflict_findings,
    ) = analyze_class_blocks(content)
    findings.extend(pub_findings)
    # static member findings: cap per file
    if len(static_findings) > CAPS.get("static_member_variable", len(static_findings)):
        cap_static = CAPS["static_member_variable"]
        static_findings = static_findings[:cap_static]
    findings.extend(static_findings)
    findings.extend(long_func_findings)
    # name conflict findings already respect cap via per-file cap on the rule
    if len(conflict_findings) > CAPS.get("member_local_name_conflict", len(conflict_findings)):
        cap = CAPS["member_local_name_conflict"]
        conflict_findings = conflict_findings[:cap]
    findings.extend(conflict_findings)

    # compute score: real_score is the unclamped value (may be negative);
    # score is clamped to a minimum of 0 for display.
    real_score = 100
    for f in findings:
        real_score -= f.deduction
    score = max(0, real_score)

    return FileReport(
        path=str(path),
        score=score,
        findings=findings,
        real_score=real_score,
    )
