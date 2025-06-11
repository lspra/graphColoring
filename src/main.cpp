#include "graph.h"
#include <random>

int main() {
  Graph a("examples/random_prob/graph5_0.40", 0);
  a.print();
  // a.color();
  // a.greedy_color();
  int clique[5];
  a.find_max_clique(clique);
  // a.greedy_clique(clique);
  // printf("%lf\n", a.lovasz());
}