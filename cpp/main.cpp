#include <chrono>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <vector>

#include <Eigen/Sparse>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>

#include <csv_utils.h>
#include <memory_xplatform.h>

#ifdef DEBUG
#define D(x) (std::cerr << x << std::endl)
#else
#define D(x) do{}while(0)
#endif

#define CSV_EOL "\n"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor, long int> SpMat;

std::ofstream output;

void analyze_matrix(std::string filename);
void print_csv_row(std::ofstream &stream, std::vector<std::string> columns) {
    stream << join(columns, ",") << CSV_EOL;
}

int main() {
    output.open("output.csv", std::ofstream::out);
    print_csv_row(output, {"filename", "size", "proc_memory_start", "proc_memory_end", "b_time", "chol_time", "relative_error"});

    std::string path = "matrix";
    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".mtx") {
            if (! output.is_open()) {
                output.open("output.csv", std::ofstream::app);
            }

            std::cout << entry.path() << std::endl;
            analyze_matrix(entry.path().string());

            output.close();
        }
    }

    return 0;
}

void analyze_matrix(std::string filename) {
    unsigned long long start_tot_virtual, start_proc_virtual, start_proc_physical, start_tot_physical,
                       end_tot_virtual, end_proc_virtual, end_proc_physical, end_tot_physical;

    D("Memory Usage (proc/total):");
    start_proc_virtual = memory::process_current_virtual();
    start_tot_virtual = memory::total_virtual();
    D("> Virtual: " << start_proc_virtual << " / " << start_tot_virtual);
    start_proc_physical = memory::process_current_physical();
    start_tot_physical = memory::total_physical();
    D("> Physical: " << start_proc_physical << " / " << start_tot_physical);
    D("");

    SpMat A;
    D("Loading matrix file: " << filename);
    Eigen::loadMarket(A, filename);

    D("Solve:");
    D("> Calculating b vector...");
    auto b_start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd x_es = Eigen::VectorXd::Ones(A.rows());
    Eigen::VectorXd b(A.rows());
    b = A*x_es;
    auto b_finish = std::chrono::high_resolution_clock::now();

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

    std::chrono::duration<double> b_time = b_finish - b_start;
    std::chrono::duration<double> chol_time = chol_finish - chol_start;

    double relative_error = (x_ap - x_es).norm() / x_es.norm();
    D("Results");
    D("> Timing (b/chol/total) in seconds: " << b_time.count() << " / " << chol_time.count() << " / " << (b_time + chol_time).count());
    D("> Relative error: " << relative_error);
    D("");

    D("Writing output...");
    output << filename << ","
        << A.rows() << ","
        << start_proc_virtual << ","
        << end_proc_virtual << ","
        << b_time.count() << ","
        << chol_time.count() << ","
        << relative_error
        << CSV_EOL;
    D("");
}
