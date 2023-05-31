# PACE 2023: TinyWidth Solver

## Description

This repository contains the source code for the TinyWidth Solver, a program that computes the twin-width of graphs, for the [2023 PACE Challenge](https://pacechallenge.org/2023/).
TinyWidth Solver focuses on graphs of small twin-width.  
**Authors**: Gabriel Bathie, Jérôme Boillot, Nicolas Bousquet, Théo Pierron. 

## Requirements

Your machine needs to have a `C++` compiler with support for `C++20`/`C++2a` installed.
Moreover, our solver makes use of AVX256 SIMD instructions, and can only run on platforms supporting them (which include almost all recent intel processors).
If you need to run it on another platform, please adapt the `long_bitset.hpp` file or contact us.

## Building & running

This repository contains file to build our solver with the `cmake` toolchain.

To build the executable, run the following commands at the root of the repository:
```bash
mkdir -p build
cd build
cmake ..
cmake --build .
```
The `build` repository will then contain the executable, named `main`.

You can the run our solver using `./main < path/to/input_file.gr`.
