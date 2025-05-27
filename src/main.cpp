#include "graph.h"
#include <random>

int main() {
  Graph a("examples/graph1", 0);
  a.print();
  a.color();
}