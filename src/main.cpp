#include "graph.h"
#include <random>

int main() {
  Graph a("examples/random_prob/graph35_0.10", 0);
  // a.print();
  // a.color();
  // a.greedy_color();
  a.find_max_clique();
}