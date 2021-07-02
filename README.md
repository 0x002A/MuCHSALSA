# LazyBastardOnMate

LazyBastard which is running faster after consuming a lot of mate iced tea.

# Building

Building the project requires cmake and can be done by running the following commands:

```bash
mkdir build
cd build
cmake ..
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