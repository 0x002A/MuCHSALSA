# MuCHSALSA

**Mu**lti-**C**ore **H**ybrid **S**hort- **A**nd **L**ong-read **S**equence **A**ssembler

Based on:
Gatter, T., von Löhneysen, S., Fallmann, J. et al. LazyB: fast and cheap genome assembly. Algorithms Mol Biol 16, 8 (2021).

# Prerequisites

`MuCHSALSA` requires the following tools to be available in the `PATH`:

- jellyfish 2
- bbduk
- Abyss 2 (using 'abyss-pe')
- minimap2
- Python 3.6.9 (at least)

# Download

A zip file containing a statically linked binary for x86-64 alongside with the assembly pipeline and associated scripts
can be downloaded from the releases page on GitHub.

# Usage

The pipeline can be run using the following command:

```bash
sh pipeline.sh [k-mer-size-filter] [k-mer-size-assembly] [name] [illumina-inputfile-1] [illumina-inputfile-2] [nanopore-inputfile] [output-folder]
```

**Note**: The level of parallelization used for parts (default value: 8) of the pipeline is set via a variable
in `pipeline.sh`.

[k-mer-size-filter] specifies the k-mer size for k-mer counting in raw illumina data. Reads with highly abundant k-mers
are removed from the data. Starting at k=50 is recommended.

[k-mer-size-assembly] specifies the k-mer size during illumina assembly (here using Abyss). Starting at k=90 is
recommended.

# Building

Building the project requires **git**, **cmake**, **clang-15**, **libc++-15**, **Doxygen**, and **Graphviz** and can be done by running the following
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
find . -regex '.*\.\(cpp\|h\)' -exec clang-format-15 -style=file -i {} \;
```

# Standard Library

This project uses `libc++-15` as standard library. For information on installing the `LLVM` ecosystem
see [here]([https://libcxx.llvm.org/docs/BuildingLibcxx.html](https://llvm.org)).
