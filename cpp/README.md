## MdCS CPP Progetto 1
Cholesky direct method with Eigen library on .mtx matrix.

### Dependencies Installation
Run the following command inside the `libs` subfolder:

```sh
git clone https://gitlab.com/libeigen/eigen.git eigen
```

### Matrix conversion
Since the given .mtx matrix files are symmetric matrixes and Eigen doesn't support parsing MTX files defined as symmetric, matrixes needs to be regenerated through Matlab with a custom script.

For further details check the `matlab` directory of this repository.

### Compilation
The program has been compiled and tested under g++ v9.2.0 on Windows and Linux and requires at least c++17 standards.

The compilation is managed by a Makefile that can be run from the `cpp` folder with the following command:
```sh
make clean && make debug
```

You can also compile an optimized version with the following command:
```sh
make clean && make release
```

In both cases the binaries will be placed in the corresponding subfolder inside `build`.

### Usage
Before the actual execution, make sure to convert the mat files in the matlab folder to mtx with the given script.
If you want to skip a matrix file for a particular execution, move it in any subfolder.

To run the program:
```sh
./build/debug/main.exe
```

or:
```sh
./build/release/main.exe
```

### Output
The program generates a CSV output that can be further analyzed and compared.

The output is composed with the following informations gathered at runtime:
- **filename**: the name of the matrix file (includes the folder name prefix).
- **size**: the number of rows.
- **memory_delta**: the process memory used to calculate the cholesky solution (in bytes).
- **chol_time**: time taken to solve the system (cholesky decomposition + solution calculation) (in seconds).
- **relative_error**: the relative error between the calculated solution and the exact ones solution.
