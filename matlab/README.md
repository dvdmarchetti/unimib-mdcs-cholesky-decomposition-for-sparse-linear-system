## MdCS Matlab Progetto 1
Cholesky direct method with Matlab and .mat files.

### Matrix download
The needed matrix .mat files can be downloaded from the following links:
- [Flan_1565](https://sparse.tamu.edu/Janna/Flan_1565)
- [StocF-1465](https://sparse.tamu.edu/Janna/StocF-1465)
- [cfd2](https://sparse.tamu.edu/Rothberg/cfd2)
- [cfd1](https://sparse.tamu.edu/Rothberg/cfd1)
- [G3_circuit](https://sparse.tamu.edu/AMD/G3_circuit)
- [parabolic_fem](https://sparse.tamu.edu/Wissgott/parabolic_fem)
- [apache2](https://sparse.tamu.edu/GHS_psdef/apache2)
- [shallow_water1](https://sparse.tamu.edu/MaxPlanck/shallow_water1)
- [ex15](https://sparse.tamu.edu/FIDAP/ex15)

### Matrix conversion
Given the issue with Eigen explained in `cpp/README.md`, the .mtx files needs to be regenerated once.

To do so, after all the files are placed in the `matrix` folder, open matlab in the current directory and execute the file `convert_mat_mtx.m`.

The process will take quite a while because all the .mat files needs to be read and fully saved to file from scratch. Once the conversion is complete, the files can be copied in the `cpp/matrix` folder to let the c++ program solve them.

### Usage
TODO.

### Output
TODO.
