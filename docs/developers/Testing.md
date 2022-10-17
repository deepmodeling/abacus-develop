# Testing

- [Adding a unittest](#adding-a-unittest)
- [Regression tests](#regression-tests)
- [Test coverage](#test-coverage)
- [Dashboard](#dashboard)
- [Compiling](#compiling)

## Adding a unittest

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

## Regression tests

How to run E2E tests, how to add a E2E test, how to understand results.

## Test coverage

Test coverage report table.

## Dashboard

Clusters, platforms supported table, respective test coverage and running time.

## Compiling

Compilation environment supportting table.
