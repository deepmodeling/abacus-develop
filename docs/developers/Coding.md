# Coding

- [Getting started](#getting-started)
- [Structure of the package](#structure-of-the-package)
- [Code conventions](#code-conventions)
- [Doxygen](#doxygen)
- [Debugging](#debugging)
- [Github respository](#github-respository)
- [Issue](#issue)
- [Commit](#commit)
- [Log](#log)
- [Profiling](#profiling)
- [Release Checklist](#release-checklist)

## Getting started

First of all, thank you for taking time making contributions to ABACUS!
This file provides a guideline for it.
Please referring to our GitHub [issue tracker](https://github.com/deepmodeling/abacus-develop/issues), and our developers are willing to help.
If you find a bug, you can help us by submitting an issue to our GitHub Repository. Even better, you can submit a Pull Request with a patch. You can request a new feature by submitting an issue to our GitHub Repository.
If you would like to implement a new feature, please submit an issue with a proposal for your work first, and that ensures your work collaborates with our development road map well. For a major feature, first open an issue and outline your proposal so that it can be discussed. This will also allow us to better coordinate our efforts, prevent duplication of work, and help you to craft the change so that it is successfully accepted into the project.

The ABACUS [Github respository](https://github.com/deepmodeling/abacus-develop). 

The ABACUS [releases](https://github.com/deepmodeling/abacus-develop/releases).

### Submitting a Pull Request

1. [Fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the ABACUS repo.

2. Pull your forked repository, and create a new git branchmake to your changes in it:

     ```shell
     git checkout -b my-fix-branch
     ```

3. Coding your patch, including appropriate test cases and docs.
To run a subset of unit test, use `ctest -R <test-match-pattern>` to perform tests with name matched by given pattern.

4. After tests passed, commit your changes [with a proper message](#commit-message-guidelines).

5. Push your branch to GitHub:

    ```shell
    git push origin my-fix-branch
    ```

6. In GitHub, send a pull request with `deepmodeling/abacus-develop:develop` as the base repository.

7. After your pull request is merged, you can safely delete your branch and sync the changes from the main (upstream) repository:

- Delete the remote branch on GitHub either through the GitHub web UI or your local shell as follows:

    ```shell
    git push origin --delete my-fix-branch
    ```

- Check out the master branch:

    ```shell
    git checkout develop -f
    ```

- Delete the local branch:

    ```shell
    git branch -D my-fix-branch
    ```

- Update your master with the latest upstream version:

    ```shell
    git pull --ff upstream develop
    ```

## Structure of the package

Please refer to [our instructions](./quick_start/easy_install.md) on how to installing ABACUS.
The source code of ABACUS is based on several modules. Under the ABACUS root directory, there are the following folders:

- `cmake`: relevant files for finding required packages when compiling the code with cmake;
- `docs`: documents and supplementary info about ABACUS;
- `examples`: some examples showing the usage of ABACUS;
- `source`: the source code in separated modules, under which a `test` folder for its unit tests;
- `tests`: End-to-end test cases;
- `tools`: the script for generating the numerical atomic orbitals.

### Core module guides

This part will introduce core modules of ABACUS in detail. It should include basic introduction, API manual and examples.

## Code conventions

### Formatting style

We use `clang-format` as our code formatter. The `.clang-format` file in root directory describes the rules to conform with.
For Visual Studio Code developers, the [official extension of C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) provided by Microsoft can help you format your codes following the rules. With this extension installed, format your code with `shift+command/alt+f`.
Configure your VS Code settings as `"C_Cpp.clang_format_style": "file"` (you can look up this option by pasting it into the search box of VS Code settings page), and all this stuff will take into effect. You may also set `"editor.formatOnSave": true` to avoid formatting files everytime manually.

### Comment Style for documentation

ABACUS uses Doxygen to generate docs directly from `.h` and `.cpp` code files.

For comments that need to be shown in documents, these formats should be used -- **Javadoc style** (as follow) is recommended, though Qt style is also ok. See it in [official manual](https://www.doxygen.nl/manual/docblocks.html).

A helpful VS Code extension -- [Doxygen Documentation Generator](https://marketplace.visualstudio.com/items?itemName=cschlosser.doxdocgen), can help you formating comments.

An practical example is class [LCAO_Deepks](https://github.com/deepmodeling/abacus-develop/blob/deepks/source/module_deepks/LCAO_deepks.h), the effects can be seen on [readthedocs page](https://abacus-deepks.readthedocs.io/en/latest/DeePKS_API/classLCAO__Descriptor.html#exhale-class-classLCAO-Descriptor)

- Tips
  - Only comments in .h file will be visible in generated  by Doxygen + Sphinx;
  - Private class members will not be documented;
  - Use [Markdown features](https://www.doxygen.nl/manual/markdown.html), such as using a empty new line for a new paragraph.

- Detailed Comment Block

    ```cpp
    /**
    * ... text ...
    */
    ```

- Brief + Detailed Comment Block

    ```cpp
    /// Brief description which ends at this dot. Details follow
    /// here.

    /// Brief description.
    /** Detailed description. */
    ```

- Comments After the Item: Add a "<"

    ```cpp
    int var; /**<Detailed description after the member */
    int var; ///<Brief description after the member
    ```

- Parameters
    usage: `[in],[out],[in,out] description`
    *e.g.*

    ```cpp
    void foo(int v/**< [in] docs for input parameter v.*/);
    ```

    or use `@param` command.

- Formula
  - inline: `\f$myformula\f$`
  - separate line: `\f[myformula\f]`
  - environment: `\f{environment}{myformula}`
  - *e.g.*

    ```latex
    \f{eqnarray*}{
            g &=& \frac{Gm_2}{r^2} \\
            &=& \frac{(6.673 \times 10^{-11}\,\mbox{m}^3\,\mbox{kg}^{-1}\,
                \mbox{s}^{-2})(5.9736 \times 10^{24}\,\mbox{kg})}{(6371.01\,\mbox{km})^2} \\
            &=& 9.82066032\,\mbox{m/s}^2
    \f}
    ```

## Doxygen

This will be the entrance of doxygen documentation.

## Debugging

This part will include several important guidlines regarding how to debug efficiently.

## Issue

Before you submit an issue, please search the [issue tracker](https://github.com/deepmodeling/abacus-develop/issues), and maybe your problem has been discussed and fixed. You can [submit new issues](https://github.com/deepmodeling/abacus-develop/issues/new/choose) by filling our issue forms.
To help us reproduce and confirm a bug, please provide a test case and building environment in your issue.


## Commit

### [Commit tracker](https://github.com/deepmodeling/abacus-develop/commits)

### Commit Message Guidelines

A well-formatted commit message leads a more readable history when we look through some changes, and helps us generate change log.
We follow up [The Conventional Commits specification](https://www.conventionalcommits.org) for commit message format.
The commit message should be structured as follows:

```text
<type>[optional scope]: <description>

[optional body]

[optional footer(s)]
```

- Header
  - type: The general intention of this commit
    - `feat`: A new feature
    - `fix`: A bug fix
    - `docs`: Only documentation changes
    - `style`: Changes that do not affect the meaning of the code
    - `refactor`: A code change that neither fixes a bug nor adds a feature
    - `perf`: A code change that improves performance
    - `test`: Adding missing tests or correcting existing tests
    - `build`: Changes that affect the build system or external dependencies
    - `ci`: Changes to our CI configuration files and scripts
    - `revert`: Reverting commits
  - scope: The scope could be the module which this commit changes; for example, `orbital`
  - description: A short summary of the code changes: tell others what you did in one sentence.
- Body: optional, providing additional contextual information about the code changes, e.g. the motivation of this commit, referenced materials, and so on.
- Footer: optional, reference GitHub issues, PRs that this commit closes or is related to

Here is an example:

```text
fix(lcao): use correct scalapack interface.

`pzgemv_` and `pzgemm_` used `double*` for alpha and beta parameters but not `complex*` , this would cause error in GNU compiler.

Fix #753.
```

## Log

This part will introduce how to trigger different ranks of log info in ABACUS, and how to use them to analyse the performance of ABACUS and to debug, et al.

## Profiling

This part is about how to analyse the time table reported in log file of ABACUS.

## Release Checklist

This part introduces the release checlist before releasing a new version or subversion of ABACUS.
