#include <chrono>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <vector>

#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <memory_xplatform.h>

#ifdef DEBUG
    #define D(x) (std::cerr << x << std::endl)
#else
    #define D(x) do{}while(0)
#endif

#define CSV_EOL "\n"

#if defined(_WIN32) || defined(__CYGWIN__)
    #define OUTPUT_FILE "windows-output.csv"
#elif defined(__linux__) || defined(unix) || defined(__unix__) || defined(__unix)
    #define OUTPUT_FILE "unix-output.csv"
#else
    #error Unknown environment!
#endif

typedef unsigned long long ull;
typedef Eigen::SparseMatrix<double> SpMat;

struct result {
    unsigned int size;
    ull memory_delta;
    std::chrono::duration<double> solve_time;
    double relative_error;
};

result analyze_matrix(std::string filename);
std::ostream& operator<<(std::ostream& stream, const result& r);

int main() {
    result _result;
    std::ofstream output(OUTPUT_FILE, std::ofstream::out);

    // Write headers
    output << "filename" << ","
        << "size" << ","
        << "memory_delta" << ","
        << "solve_time" << ","
        << "relative_error" << CSV_EOL;

    // Look for matrix in ../matlab/matrix_mtx folder
    std::string path = "../matlab/matrix_mtx";
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".mtx") {
            if (! output.is_open()) {
                output.open(OUTPUT_FILE, std::ofstream::app);
            }

            std::cout << entry.path() << std::endl;
            _result = analyze_matrix(entry.path().string());

            D("Writing output...");
            output << entry.path().stem() << "," << _result << CSV_EOL;
            D("");

            output.close();
        }
    }

    return 0;
}

/**
 * Analyze a single matrix given it's filename
 * The matrix will be imported assuming it's in .mtx format.
 *
 * @param std::string filename
 * @return result
 */
result analyze_matrix(std::string filename) {
    result r;
    ull start_tot_virtual, start_proc_virtual, start_proc_physical, start_tot_physical,
        end_tot_virtual, end_proc_virtual, end_proc_physical, end_tot_physical;

    SpMat A; // Eigen::SparseMatrix<double>
    D("Loading matrix file: " << filename);
    Eigen::loadMarket(A, filename);

    // Debug memory usage to cout
    D("Memory Usage (proc/total):");
    start_proc_virtual = memory::process_current_virtual();
    start_tot_virtual = memory::total_virtual();
    D("> Virtual: " << start_proc_virtual << " / " << start_tot_virtual);
    start_proc_physical = memory::process_current_physical();
    start_tot_physical = memory::total_physical();
    D("> Physical: " << start_proc_physical << " / " << start_tot_physical);
    D("");

    D("Solve:");
    D("> Calculating b vector...");
    Eigen::VectorXd x_es = Eigen::VectorXd::Ones(A.rows());
    Eigen::VectorXd b(A.rows());
    b = A*x_es;

    D("> Applying CholeskySimplicial solver...");
    auto chol_start = std::chrono::high_resolution_clock::now();
    Eigen::SimplicialCholesky<SpMat> chol(A);
    Eigen::VectorXd x_ap = chol.solve(b);
    auto chol_finish = std::chrono::high_resolution_clock::now();
    D("");

    D("Memory Usage (proc/total):");
    end_proc_virtual = memory::process_current_virtual();
    end_tot_virtual = memory::total_virtual();
    D("> Virtual: " << end_proc_virtual << " / " << end_tot_virtual);
    end_proc_physical = memory::process_current_physical();
    end_tot_physical = memory::total_physical();
    D("> Physical: " << end_proc_physical << " / " << end_tot_physical);
    D("");

    r.size = A.rows();
    r.memory_delta = end_proc_physical - start_proc_physical;
    r.solve_time = chol_finish - chol_start;
    r.relative_error = (x_ap - x_es).norm() / x_es.norm();

    D("Results");
    D("> Solve time in seconds: " << r.solve_time.count());
    D("> Relative error: " << r.relative_error);
    D("");

    return r;
}

/**
 * Override the << operator to allow an easy print of
 * the results struct to stdout or output file.
 *
 * @params std::ostream& stream
 * @params const result& r
 * @return std::ostream&
 */
std::ostream& operator<<(std::ostream& stream, const result& r) {
    stream << r.size << ","
        << r.memory_delta << ","
        << r.solve_time.count() << ","
        << r.relative_error;

    return stream;
}
