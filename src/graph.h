#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>
#include <random>

#include "dsdp/dsdp5.h"

#define N 1000

class Graph {
  public:
    Graph(std::string filepath, int type);
    void print();
    void color();
    void greedy_color();
    void find_max_clique();

  private:
    unsigned int n, m;
    double maxr = 0.0;
    double k;
    std::vector<int> edges[N];
    std::vector<int> not_edges[N];

    double values[2][3];
    int indices[2*N][3];
    double* vecColoring;
    int colors[N];
    std::random_device rd;
    std::mt19937 gen;
    std::normal_distribution<double> ndist;

    void load_adjecency(std::ifstream& fs);
    void vector_color();
};