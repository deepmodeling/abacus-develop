# Basic Tool Classes Design Guide for ABACUS

## Overview

This document provides guidelines for designing and implementing basic tool classes in the ABACUS codebase, focusing on best practices for memory management, code style, and testing. These guidelines apply to all basic mathematical and utility classes, including but not limited to:

- vector3.h
- matrix.h
- timer.h
- ndarray.h
- realarray.h
- complexarray.h
- complexmatrix.h
- matrix3.h
- intarray.h
- formatter.h
- math_chebyshev.h

While this guide uses `IntArray` as an example for illustration purposes, the principles and practices described here are applicable to all basic tool classes in ABACUS.

## Memory Management

### 1. Exception Handling for Memory Allocation

Always use try-catch blocks when allocating memory to handle `std::bad_alloc` exceptions gracefully:

### 2. Two-Stage Memory Allocation

When reallocating memory (e.g., in `create` methods), use a two-stage approach to ensure that the original object remains valid if memory allocation fails.

### 3. Null Pointer Checks

Always check for null pointers before accessing memory, especially in methods that might be called on objects with failed memory allocation.

## Class Design

### 1. Copy Constructor

Implement a copy constructor to avoid shallow copy issues.

### 2. Move Semantics

Implement move constructor and move assignment operator to improve performance.

### 3. Boundary Checks

Add boundary checks to prevent out-of-bounds access.

## Code Style

### 1. Brace Style

Use separate lines for braces, and always use braces for "if" and "for" statements, even if they contain one line of code

### 2. Indentation

Use spaces instead of tabs for indentation (4 spaces per indent level).

### 3. Comments

Use English for comments and document important functionality. Follow Doxygen-style documentation for classes and methods.

## Code Quality

### 1. Named Constants

Avoid using magic numbers. Instead, define named constants for numerical values:

### 2. Header Includes

Ensure all necessary header files are included, especially for functions like `assert`:

```cpp
#include <cassert>
```

## Testing

### 1. Unit Tests

Write comprehensive unit tests for all classes, including:
- Constructor tests
- Method tests
- Exception handling tests
- Edge case tests

### 2. Test Class Initialization

Use constructor initialization lists for test classes to improve compatibility:

```cpp
class IntArrayTest : public testing::Test
{
protected:
    ModuleBase::IntArray a2, a3, a4, a5, a6;
    int aa;
    int bb;
    int count0;
    int count1;
    const int zero;

    IntArrayTest() : aa(11), bb(1), zero(0)
    {
    }
};
```

## Best Practices

1. **Single Responsibility Principle**: Each class should have a single, well-defined responsibility.
2. **Encapsulation**: Hide implementation details and expose only necessary interfaces.
3. **Error Handling**: Handle errors gracefully, especially memory allocation failures.
4. **Performance**: Use move semantics and other performance optimizations where appropriate.
5. **Testing**: Write comprehensive tests for all functionality.
6. **Code Style**: Follow consistent code style guidelines, including:
   - Always use braces for if and for statements
   - Use separate lines for braces
   - Use spaces instead of tabs for indentation
   - Use English for comments
7. **Code Quality**: Maintain high code quality by:
   - Using named constants instead of magic numbers
   - Ensuring all necessary header files are included
   - Adding boundary checks to prevent out-of-bounds access
8. **Documentation**: Document classes and methods to improve maintainability.
9. **Compatibility**: Ensure code is compatible with C++11 standard.
10. **Portability**: Write code that works across different platforms.
11. **Reusability**: Design classes to be reusable in different contexts.

## Application to Other Basic Tool Classes

While this guide uses `IntArray` as an example, these principles apply to all basic tool classes in ABACUS. For example:

- **vector3.h**: Apply the same memory management and error handling principles, with additional focus on vector operations and operator overloading.
- **matrix.h**: Extend the memory management practices to 2D arrays, with additional considerations for matrix operations.
- **timer.h**: Focus on static member management and time measurement accuracy.
- **ndarray.h**: Apply the same principles to multi-dimensional arrays, with additional considerations for shape manipulation.
- **formatter.h**: Focus on string manipulation and formatting, with attention to performance and usability.
- **math_chebyshev.h**: Apply the principles to template classes, with additional focus on mathematical algorithm implementation.

By following these guidelines, you can ensure that all basic tool classes in ABACUS are well-designed, robust, and maintainable.
