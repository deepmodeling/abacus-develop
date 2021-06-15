
## Style of Comments for Doxygen
for comments that need to be shown in documents, use this format: 
Javadoc style recommended, though Qt style is ok. see it in [official manual](https://www.doxygen.nl/manual/docblocks.html)
VS Code Extension: Doxygen Documentation Generator

### detailed comment block
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
### brief + detailed
```
/// Brief description which ends at this dot. Details follow
/// here.

/// Brief description.
/** Detailed description. */
```
​

### Comments after the item：add a "<"

```
int var; /**<Detailed description after the member */
int var; ///<Brief description after the member
```

### Parameter
usage: `[in],[out],[in,out] description`
eg. 
```
void foo(int v/**< [in] docs for input parameter v.*/);
```
or use `@param` commond
​
### list
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
### Formula

inline: 	`\f$myformula\f$`
separate line: 	 `\f[myformula\f]`
environment:	`\f{environment}{myformula}`
eg.
```
   \f{eqnarray*}{
        g &=& \frac{Gm_2}{r^2} \\ 
          &=& \frac{(6.673 \times 10^{-11}\,\mbox{m}^3\,\mbox{kg}^{-1}\,
              \mbox{s}^{-2})(5.9736 \times 10^{24}\,\mbox{kg})}{(6371.01\,\mbox{km})^2} \\ 
          &=& 9.82066032\,\mbox{m/s}^2
   \f}
```

## Tips
- only comments in .h file will be visible in generated  by doxygen + Sphinx; 
- private class members will not be documented; 
- use [Markdown features](https://www.doxygen.nl/manual/markdown.html), such as using a empty new line for a new paragraph. 
​