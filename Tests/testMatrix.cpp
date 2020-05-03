#include <iostream>
#include <vector>
#include <armadillo>
#include <armadillo_bits/fn_chol.hpp>

using namespace std;
using namespace arma;

int main() {
    sp_mat A;
    string filename = "./matrix/ex15.mtx";
    //string filename = "./matrix/cfd1.mtx";
    //string filename = "./matrix/nos1.mtx";

    // load matrix only value != 0
    cout << "Load matrix: " << filename <<endl;
    cout << "format: MatLab " << endl;
    A.load(filename, coord_ascii);  

    // In Armadillo matrix indices start at 0 (due to C++ conventions),
    // while in the matrix market file they start at 1
    A = A.tail_rows(A.n_rows - 1);
    A = A.tail_cols(A.n_cols - 1);

    cout << "## Creating xe vector ##" << endl;
    dvec x_es = ones(A.n_rows);
    
    cout << "## Calculating b vector terms ##" << endl;
    dvec b(A.n_rows);
    b = A*x_es;

    // cout << A << endl;
    // cout << xe << endl;
    // cout << b << endl; 
    try {
        cout << "## Solve matrix ##" << endl;
        sp_mat R;
        R = chol(A);
    }
    catch (std::runtime_error &e) {

        cout << "Error Chol: " << e.what() << endl;

    }
    //dvec x_ap = spsolve(_chol, b);

    // cout << "## Results ##" << endl;
    // cout << "Errore relativo: " << normalise(x_es - x_ap) / normalise(x_es) << endl; 
    
    return 0;
}