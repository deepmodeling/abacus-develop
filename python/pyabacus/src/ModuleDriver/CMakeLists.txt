# CMakeLists.txt for ModuleDriver Python bindings
# This module provides Python bindings for the ABACUS Driver
#
# IMPORTANT: This module requires the full ABACUS library to be built with -fPIC.
# Currently, the standard ABACUS build does not support this.
#
# To use pyabacus.abacus(), you need to:
# 1. Build ABACUS as a shared library (not yet supported upstream)
# 2. Or use the subprocess-based fallback implementation
#
# For now, this module is disabled and the Python-side implementation
# uses subprocess to call the abacus executable.

message(STATUS "PyABACUS driver module: C++ bindings require ABACUS to be built as a shared library")
message(STATUS "PyABACUS driver module: Using subprocess-based fallback implementation instead")
message(STATUS "PyABACUS driver module: Skipping C++ driver build")

# The driver module is not built - Python fallback is used instead
return()


