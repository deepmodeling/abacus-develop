# For developers

- [Raising issues on GitHub](#raising-issues-on-github)
- [Modularization and module tests](#modularization-and-module-tests)
- [Contributing to ABACUS](#contributing-to-abacus)
  - [Making pull requests](#making-pull-requests)
  - [Providing uniy tests](#providing-unit-tests)
  - [Upating documentation](#updating-documentation)
  - [Macros](#macros)
  - [Comment style for documentation](#comment-style-for-documentation)
  - [Code Formatting style](#code-formatting-style)
    [back to main page](../README.md)

## Creating issues on GitHub

Creating issues on GitHub is a convernient way to notify the develper team about bugs and feature requests of the ABACUS code. We provide a few templates for issues.

[back to top](#for-developers)

## Modularization and module tests

The ABACUS code is refactored to several self-contained modules. A description of modules can be found in the [installation guide](install.md#structure-of-source-code). We also provide module unit tests.

### Add a unit test

We use GoogleTest as our test framework. Write your test under the corresponding module folder at `abacus-develop/tests`, then append the test to `tests/CMakeLists.txt`. If there are currently no unit tests provided for the module, do as follows. `module_base` provides a simple demonstration.

- Add a folder named `test` under the module.
- Append the content below to `CMakeLists.txt` of the module:

```cmake
IF (BUILD_TESTING)
  add_subdirectory(test)
endif()
```

- Add a blank `CMakeLists.txt` under `module*/test`.

To add a unit test:

- Write your test under `GoogleTest` framework.
- Add your testing source code with suffix `*_test.cpp` in `test` directory.
- Append the content below to `CMakeLists.txt` of the module:

```cmake
AddTest(
  TARGET <module_name>_<test_name> # this is the executable file name of the test
  SOURCES <test_name>.cpp

  # OPTIONAL: if this test requires external libraries, add them with "LIBS" statement.
  LIBS math_libs # `math_libs` includes all math libraries in ABACUS.
)
```

- Build with `-D BUILD_TESTING=1` flag. You can find built testing programs under `build/source/<module_name>/test`.
- Follow the installing procedure of CMake. The tests will move to `build/test`.

[back to top](#for-developers)

## Contributing to ABACUS

### Making pull requests

1. [Fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the abacus repo.

2. In your forked repository, make your changes in a new git branch:

     ```shell
     git checkout -b my-fix-branch
     ```

3. Create your patch, including appropriate test cases.

4. Commit Your changes.

5. Push your branch to GitHub:

    ```shell
    git push origin my-fix-branch
    ```

6. In GitHub, send a pull request to `abacus-develop:develop`.

7. After your pull request is merged

After your pull request is merged, you can safely delete your branch and pull the changes from the main (upstream) repository:

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


### Updating documentation

### Macros

### Comment Style for Documentation

ABACUS uses Doxygen to generate docs directly from `.h` and `.cpp` code files.

For comments that need to be shown in documents, these formats should be used -- **Javadoc style** (as follow) is recommended, though Qt style is also ok. See it in [official manual](https://www.doxygen.nl/manual/docblocks.html).

A helpful VS Code extension -- Doxygen Documentation Generator, can help you formating comments.

An practical example is class [LCAO_Deepks](https://github.com/deepmodeling/abacus-develop/blob/deepks/source/module_deepks/LCAO_deepks.h), the effects can be seen on [readthedocs page](https://abacus-deepks.readthedocs.io/en/latest/DeePKS_API/classLCAO__Descriptor.html#exhale-class-classLCAO-Descriptor)

- Detailed Comment Block

    ```cpp
    /**
    * ... text ...
    */
    ```

    ```cpp
    ///
    /// ... text ...
    ///
    ```

    or ( set JAVADOC_BANNER = YES )

    ```cpp
    /********************************************//**
    *  ... text
    ***********************************************/

    /////////////////////////////////////////////////
    /// ... text ...
    /////////////////////////////////////////////////

    /************************************************
    *  ... text
    ***********************************************/
    ```

- Brief + Detailed Comment Block

    ```cpp
    /// Brief description which ends at this dot. Details follow
    /// here.

    /// Brief description.
    /** Detailed description. */
    ```

​

- Comments After the Item: Add a "<"

    ```cpp
    int var; /**<Detailed description after the member */
    int var; ///<Brief description after the member
    ```

- Parameters
    usage: `[in],[out],[in,out] description`
    eg.

    ```cpp
    void foo(int v/**< [in] docs for input parameter v.*/);
    ```

    or use `@param` commond
​
- List
    eg.1

    ```cpp
    /**
    * Text before the list
    * - list item 1
    *   - sub item 1
    *     - sub sub item 1
    *     - sub sub item 2
    *     .
    *     The dot above ends the sub sub item list.
    *
    *     More text for the first sub item
    *   .
    *   The dot above ends the first sub item.
    *
    *   More text for the first list item
    *   - sub item 2
    *   - sub item 3
    * - list item 2
    * .
    * More text in the same paragraph.
    *
    * More text in a new paragraph.
    */
    ```

    eg.2

    ```cpp
    /*!
    *  A list of events:
    *    - mouse events
    *         -# mouse move event
    *         -# mouse click event\n
    *            More info about the click event.
    *         -# mouse double click event
    *    - keyboard events
    *         1. key down event
    *         2. key up event
    *
    *  More text here.
    */
    ```

- Formula

  - inline:  `\f$myformula\f$`
  - separate line:   `\f[myformula\f]`
  - environment: `\f{environment}{myformula}`
  - eg.

    ```latex
    \f{eqnarray*}{
            g &=& \frac{Gm_2}{r^2} \\
            &=& \frac{(6.673 \times 10^{-11}\,\mbox{m}^3\,\mbox{kg}^{-1}\,
                \mbox{s}^{-2})(5.9736 \times 10^{24}\,\mbox{kg})}{(6371.01\,\mbox{km})^2} \\
            &=& 9.82066032\,\mbox{m/s}^2
    \f}
    ```

- Tips
  - Only comments in .h file will be visible in generated  by Doxygen + Sphinx;
  - Private class members will not be documented;
  - Use [Markdown features](https://www.doxygen.nl/manual/markdown.html), such as using a empty new line for a new paragraph.
​
[back to top](#for-developers)

### Code formatting style

We use `clang-format` as our code formatter. The `.clang-format` file in root directory describes the rules to conform.
For Visual Studio Code developers, the [official extension of C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) provided by Microsoft can help you format your codes following the rules. Configure your VS Code settings as `"C_Cpp.clang_format_style": "file"` (you can look up this option by pasting it into the search box of VS Code settings page), and the clang-format will take into effect. You may also set `"editor.formatOnSave": true` to avoid formatting files everytime manually.
