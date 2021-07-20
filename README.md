# MuCHSALSA

**Mu**lti-**C**ore **H**ybrid **S**hort- **A**nd **L**ong-read **S**equence **A**ssembler

Based on:
Gatter, Thomas, et al. "Economic genome assembly from low coverage Illumina and Nanopore data." 20th International
Workshop on Algorithms in Bioinformatics (WABI 2020). Schloss Dagstuhl-Leibniz-Zentrum für Informatik, 2020.

# Binary

A statically linked binary for x86-64 can be downloaded from the releases page on GitHub.

# Usage

In contrast to [LazyB](https://github.com/TGatter/LazyB) this tool has an additional parameter to control the level of
parallelization. Simply plug it into the LazyB-pipeline and set the additional parameter.
**ATTENTION**: If no value is given all available cores will be used.

# Building

Building the project requires **cmake**, **clang-12** and **libc++-12** and can be done by running the following
commands:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If you want to link _libc++_ statically:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_BUILD_FAT_EXE=1 ..
make
```

Tests can be executed by running the following command afterwards:

```bash
make test
```

The library documentation can be built using the command:

```bash
make doc
```

# Formatting

All the code within this project is uniformly formatted using clang-format:

```bash
find . -regex '.*\.\(cpp\|h\)' -exec clang-format-12 -style=file -i {} \;
```

# Standard Library

This project uses _libc++_ as standard library. For information on building _libc++_
see [here](https://libcxx.llvm.org/docs/BuildingLibcxx.html).