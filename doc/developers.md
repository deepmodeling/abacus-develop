# For developers

- [Raising issues on GitHub](#raising-issues-on-github)
- [Modularization and module tests](#modularization-and-module-tests)
- [Contributing to ABACUS](#contributing-to-abacus)
    - [Making pull requests](#making-pull-requests)
    - [Providing uniy tests](#providing-unit-tests)
    - [Upating documentation](#updating-documentation)
    - [Macros](#macros)
    - [Comment style for documentation](#comment-style-for-documentation)

[back to main page](../README.md)

## Raising issues on GitHub
Raise issues on GitHub is a convernient way to notify the develper team about bugs and feature requests of the ABACUS code. We provide a few templates for issues.

[back to top](#for-developers)

## Modularization and module tests
The ABACUS code is reconstructed to form several self-contained modules. A description of modules can be found in the [installation guide](install.md#structure-of-source-code). We also provide module tests for the modules.

[back to top](#for-developers)

## Contributing to ABACUS

### Making pull requests

### Providing unit tests

### Updating documentation

### Macros

### Comment Style for Documentation
ABACUS uses Doxygen to generate docs directly from `.h ` and `.cpp` code files. 

For comments that need to be shown in documents, these formats should be used -- **Javadoc style** (as follow) is recommended, though Qt style is also ok. See it in [official manual](https://www.doxygen.nl/manual/docblocks.html).

A helpful VS Code extension -- Doxygen Documentation Generator, can help you formating comments.

- Detailed Comment Block
    ```
    /**
    * ... text ...
    */
    ```
    ```
    ///
    /// ... text ...
    ///
    ```
    or ( set JAVADOC_BANNER = YES )
    ```
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
    ```
    /// Brief description which ends at this dot. Details follow
    /// here.

    /// Brief description.
    /** Detailed description. */
    ```
​

- Comments After the Item: Add a "<"
    ```
    int var; /**<Detailed description after the member */
    int var; ///<Brief description after the member
    ```

- Parameters
    usage: `[in],[out],[in,out] description`
    eg. 
    ```
    void foo(int v/**< [in] docs for input parameter v.*/);
    ```
    or use `@param` commond
​
- List
    eg.1
    ```
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
    ```
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

    - inline: 	`\f$myformula\f$`
    - separate line: 	 `\f[myformula\f]`
    - environment:	`\f{environment}{myformula}`
    - eg.
    ```
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