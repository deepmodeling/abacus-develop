cmake_minimum_required(VERSION 3.14)

project(
  ABACUS
  DESCRIPTION "ABACUS build system experiment"
  LANGUAGES CXX
  VERSION 0.0.1
)

set(CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED True)

add_subdirectory(ABACUS.develop/source)
