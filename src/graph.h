#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>
#include <random>

#include "dsdp/dsdp5.h"

#define N 3000

class Graph {
  public:
    Graph(std::string filepath, int type);
    void print();
    void color();
    void greedy_color();
    void find_max_clique();

  private:
    unsigned int n, m, maxr_ind;
    double maxr = 0.0;
    double k;
    std::vector<int> edges[N];
    std::vector<int> not_edges[N];

    std::random_device rd;
    std::mt19937 gen;
    std::normal_distribution<double> ndist;
    double values[2][3];
    int indices[2*N][3];
    double* vecColoring;
    int colors[N];

    void load_adjecency(std::ifstream& fs);
    void vector_color();
    void check_colors_valid();
    void check_clique_valid(int* clique, int len);
};