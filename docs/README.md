# ABACUS Documentation

These are the sources of the [ABACUS documentation](https://abacus.deepmodeling.com/). They are
built by Read the Docs according to [this configuration file](../.readthedocs.yaml).

To build a local version of the documentation, perform the following steps from the `docs/`
directory:

1. Create and activate a [virtual Python environment](https://docs.python.org/3/tutorial/venv.html):

   ```bash
   python3 -m venv ../docs_venv
   source ../docs_venv/bin/activate
   ```

1. Install the required Python packages:

   ```bash
   pip3 install -r ./requirements.txt
   ```

1. (optional but recommended) Build an ABACUS binary and use it to generate the `parameters.yaml`
   file:

   ```bash
   ../bin/abacus --generate-parameters-yaml > ./parameters.yaml
   ```

1. (optional but recommended) Generate Markdown pages from the `parameters.yaml` file:

   ```bash
   python3 ./generate_input_main.py ./parameters.yaml
   ```

1. Run Sphinx:

   ```bash
   make html
   ```

1. Browse the HTML output in the `build/html` directory.

______________________________________________________________________

# Syntax Cheat Sheet

The ABACUS documentation uses Sphinx with the [MyST parser](https://myst-parser.readthedocs.io) for
Markdown support. The following gives a quick overview of the syntax:

## Headings

```
# A first-level heading

## A second-level heading

### A third-level heading
```

## Basic Text Formatting

```
**bold text**

_italic text_ (alternatively, *italic text*)

~~strikethrough text~~

`inline code`
```

For all typography options, see the
[MyST documentation](https://myst-parser.readthedocs.io/en/latest/syntax/typography.html).

## Links

- Another page: `[](advanced/install.md)`
  - Use a relative path. The `.md` file extension can be specified or omitted when linking to the
    page itself.
- Subsection in the current page: `[](#cross-references)`
- Subsection in another page: `[](advanced/install.md#build-with-cuda-support)`
  - Use a relative path and include the `.md` suffix before the `#` sign. Headings up to level four
    have anchors generated automatically; see [](#cross-references) below for a more stable
    alternative.
- External URL: `<https://abacus.deepmodeling.com/>`
- External URL with label: `[ABACUS documentation](https://abacus.deepmodeling.com/)`

## Cross References

For references that should remain stable when headings or files are renamed, prefer an explicit MyST
target over an automatically generated heading anchor. Define a target immediately before the
heading:

```text
(input-structure)=
### Structure of the INPUT file
```

Reference it from the same or another page with the Sphinx `ref` role:

```text
{ref}`input-structure`
```

A custom link label can be specified as well:

```text
{ref}`INPUT structure <input-structure>`
```

Explicit target names are global within the documentation and should therefore be unique. Use
ordinary relative Markdown links for simple page links and explicit targets for cross-references
that are expected to be long-lived.

For more details, see the
[MyST cross-reference documentation](https://myst-parser.readthedocs.io/en/latest/syntax/cross-referencing.html).

## Lists

For a numbered list:

```
1. First enumerated item
1. Second enumerated item
1. And the third item
```

> [!NOTE]
>
> Every item uses `1.` intentionally in the markdown source file, and the formatting tool in the
> precommit check will apply this style if detected. The list will still be rendered with the
> intended numbers as indices on github and the final HTML documentation page, but there is no
> longer the need to track and edit the numbers manually. For more information, see `mdformat` docs
> on [ordered lists](https://mdformat.readthedocs.io/en/stable/users/style.html#ordered-lists).

For an unordered list:

```
- A bullet point
- Another bullet point
  - Indented bullet point
- Yet another bullet point
* An asterisk is also okay
```

When nesting levels of lists, watch out for indentations and newlines.

```
1. 3D periodicity
  - XYZ
1. 2D periodicity
  - XY
  - YZ
  - XZ
1. 1D periodicity
  - X
  - Y
  - Z
1. non-periodic
```

## Tables

```
| foo | bar |
| --- | --- |
| baz | bim |
```

For more table formatting options, see the
[MyST documentation](https://myst-parser.readthedocs.io/en/latest/syntax/tables.html).

## Math

```
Inline math: $A_{ia,jb}$.

Math block:
$$ \begin{align}
    A_{ia,jb} &= (\varepsilon_a^{GW}-\varepsilon_i^{GW})\delta_{ij}\delta_{ab}
    B_{ia,jb} &= 2 v_{ia,bj} - W_{ib,aj} \quad .
\end{align} $$
```

See also the
[MyST](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#math-shortcuts) and
[MathJax](https://docs.mathjax.org/en/latest/input/tex/index.html) documentation.

## Notes and Warnings

````
```{note}
A note box.
```

```{warning}
A warning box.
```
````

For all available admonitions see the
[MyST documentation](https://myst-parser.readthedocs.io/en/latest/syntax/admonitions.html).

## Code Blocks

````
```python
for i in range(10):
  print("Hello World")
```
````

````
```text
calculation             cell-relax
symmetry                1
basis_type              lcao
ecutwfc                 100
```
````

The language identifiers like `python` determine syntax highlighting; `text` is the choice for a
plain display. Details can be found at
[MyST documentation](https://myst-parser.readthedocs.io/en/latest/syntax/code_and_apis.html).
