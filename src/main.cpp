#include "graph.h"

int main() {
  Graph a("examples/graph1", 0);
  a.print();
  a.make_SDPA();
}