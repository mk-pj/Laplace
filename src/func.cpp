#include <map>
#include "func.h"
#include "matplotlibcpp.h"

inline double boundary_condition(const int i, const int j) {
    if (i == 0 && (j >= 1 && j <= inner)) return 150;
    if (j == 0 && (i >= 1 && i <= inner)) return 100;
    if (i == upper_bound && (j >= 1 && j <= inner)) return 200;
    if (j == upper_bound && (i >= 1 && i <= inner)) return 50;
    return 0.0;
}

void init_system(sparse_matrix &A, vector<double> &b) {
    A.resize(eq_count);
    b.resize(eq_count, 0.0);
    for (int i = 1; i < upper_bound; ++i) {
        for (int j = 1; j < upper_bound; ++j) {
            int ix = (i - 1) * inner + (j - 1);
            A[ix][ix] = -4;
            pair<int, int> index_pairs[] = {
                {i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1 }
            };
            for (auto& [ni, nj] : index_pairs) {
                if (ni > 0 && ni < upper_bound && nj > 0 && nj < upper_bound) {
                    int n_ix = (ni - 1) * inner + (nj - 1);
                    A[ix][n_ix] = 1.0;
                } else {
                    b[ix] -= boundary_condition(ni, nj);
                }
            }
        }
    }
}

void gauss_jordan(sparse_matrix &A, vector<double> &b, const unsigned N, const double eps) {
    for (int i = 0; i < N; ++i) {
        bool pivot_found = abs(A[i][i]) > eps;
        if (!pivot_found) {
            for (int j = i + 1; j < N; ++j) {
                if (abs(A[j][i]) > eps) {
                    std::swap(A[i], A[j]);
                    std::swap(b[i], b[j]);
                    pivot_found = true;
                    break;
                }
            }
        }
        if (!pivot_found)
            continue;

        const double pivot = A[i][i];
        for (auto& [ind, val] : A[i])
            val /= pivot;
        b[i] /= pivot;

       for (int j = 0; j < N; ++j) {
           if (i == j) continue;
           if (A[j].count(i) == 0)continue;
           const double factor = A[j][i];
           for (auto& [k, val] : A[i]) {
               A[j][k] -= factor * val;
               if (abs(A[j][k]) < eps)
                   A[j].erase(k);
           }
           b[j] -= factor * b[i];
       }
    }
}


void plot_heatmap(const std::vector<double>& b, const unsigned N) {
    namespace plt = matplotlibcpp;
    std::vector<std::vector<double>> grid(N, std::vector<double>(N, 0.0));

    for (int i = 1; i < upper_bound; ++i)
        for (int j = 1; j < upper_bound; ++j)
            grid[i][j] = b[(i - 1) * inner + (j - 1)];


    for (int k = 1; k < upper_bound; ++k) {
        grid[0][k] = 150;
        grid[k][0] = 100;
        grid[upper_bound][k] = 200;
        grid[k][upper_bound] =  50;
    }

    grid[0][0] = 125;
    grid[0][upper_bound] = 150;
    grid[upper_bound][0] = 100;
    grid[upper_bound][upper_bound] = 125;

    plt::imshow(grid, {{"cmap","jet"},
                {"origin","lower"} });
    plt::colorbar();
    plt::title("Laplace heatmap");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::show();
}
