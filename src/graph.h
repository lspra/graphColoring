#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

#include "sdpa/sdpa_call.h"

#define N 1000

class Graph {
  private:
    unsigned int n, m;
    std::vector<int> edges[N];

    SDPA repr;

    void load_adjecency(std::ifstream& fs);
    
  public:
    void make_SDPA();
    Graph(std::string filepath, int type);
    void print();
};