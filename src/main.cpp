#include "func.h"

int main() {
    vector<double> b;
    sparse_matrix matrix;

    init_system(matrix, b);
    gauss_jordan(matrix, b, eq_count, 1e-8);
    plot_heatmap(b, n);

    return 0;
}
