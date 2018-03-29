# Setup

## Autotools
```
cd  $(BIGMPI_SOURCE_DIR) && ./autogen.sh
cd $(BIGMPI_BUILD_DIR) && $(BIGMPI_SOURCE_DIR)/configure CC=mpicc
```
(add `--prefix` if you want to `make install` somewhere)

## CMake
```
cmake $(BIGMPI_SOURCE_DIR) -DCMAKE_C_COMPILER=mpicc
```

# Build

```
make
make check
```
