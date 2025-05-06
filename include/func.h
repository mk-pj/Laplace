#ifndef FUNC_H
#define FUNC_H


#include <cmath>
#include <vector>
#include <unordered_map>

using std::abs, std::pair, std::vector, std::unordered_map;

using sparse_matrix = vector<unordered_map<int, double>>;

constexpr int n = 41;
constexpr int upper_bound = n - 1;
constexpr int inner = n - 2;
constexpr int eq_count = inner * inner;

double boundary_condition(int, int);
void init_system(sparse_matrix&, vector<double>&);
void gauss_jordan(sparse_matrix&, vector<double>&, unsigned, double);
void plot_heatmap(const vector<double>&, unsigned);
void heatmap(const vector<double>&, unsigned);


#endif //FUNC_H
