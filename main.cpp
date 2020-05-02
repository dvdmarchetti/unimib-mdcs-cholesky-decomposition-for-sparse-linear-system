#include <chrono>
#include <iostream>
#include <vector>

#include <Eigen/Sparse>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>

#define DEBUG(x)  if (debugging_enabled) { std::cout << x << std::endl; }

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

// void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);
// void saveAsBitmap(const Eigen::VectorXd& x, int n, const char* filename);

int main()
{
    SpMat A;
    std::string filename = "matrix/nos1.mtx";		//nomefile (path)

    // A.from
    // A.1,0,0, 0,1,0, 0,0,1;

    std::cout << "## Loading matrix file: " << filename << " ##";
    Eigen::loadMarket(A, filename);

    std::cout << "## Creating x_es vector ##" << std::endl;
    Eigen::VectorXd x_es = Eigen::VectorXd::Ones(A.rows());	//soluzione esatta: vettore di 1 con N righe

    std::cout << "## Calculating b vector terms ##" << std::endl;
    Eigen::VectorXd b(A.rows());	//calcolo vettore termini noti
    b = A*x_es;
	
//stampa soluzioni
    std::cout << A << std::endl;
    std::cout << x_es << std::endl;
    std::cout << b << std::endl;

    std::cout << "## Solve matrix ##" << std::endl;
    Eigen::SimplicialCholesky<SpMat> chol(A);	//creo risolutore per A
    Eigen::VectorXd x_ap = chol.solve(b);	//lo uso su B
//x_ap vettore X = soluzione chol
    std::cout << "## Results ##" << std::endl;
    std::cout << "Errore relativo: " << (x_es - x_ap).norm() / x_es.norm() << std::endl;

    return 0;
}
