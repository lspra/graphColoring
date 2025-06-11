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
    ~Graph();
    void print();
    int color(int* colors, double delta);
    int greedy_color(int* colors);
    int find_max_clique(int* clique);
    int greedy_clique(int* clique);
    double lovasz();
    void color_init(int* colors);
    double get_delta();
    unsigned next_color(int* colors, int nc, unsigned colored, double delta);
    unsigned int n, m, m_;
    std::vector<int> edges[N];

  private:
    unsigned int maxr_ind;
    double maxr = 0.0;
    double k;
    std::vector<int> not_edges[N];

    std::random_device rd;
    std::mt19937 gen;
    std::normal_distribution<double> ndist;
    double values[2][3];
    int indices[2*N][3];
    double* vecColoring;
    bool vecColor = 0;

    void load_adjecency(std::ifstream& fs);
    void vector_color();
    bool check_colors_valid(int* color, int len);
    bool check_clique_valid(int* clique, int len);
};