#include "graph.h"

unsigned uint_from_str(std::string val, int index, int* next_index) {
  unsigned ret = 0;
  while (*next_index < val.length() && !iswspace(val[*next_index])) {
    ret *= 10;
    ret += (val[*next_index] - '0');
    (*next_index)++;
  }
  return ret;
}

void Graph::load_adjecency(std::ifstream& fs) {
  std::string line;
  if(std::getline(fs, line)) {
    int index = 0;
    n = uint_from_str(line, index, &index);
  }
  int i = 0;
  while (std::getline(fs, line) && i < n) {
    if(line.length() < 2*n - 1) {
      fprintf(stderr, "Invalid file format line %d too short %d %d\n", i, line.length(), 2*n);
      exit(1);
    }
    for(int j = i; j < n; j++) {
      if(line[j*2] == '1') {
        edges[i].push_back(j);
        m++;
      }
    }
    i++;
  }
  if (i < n) {
    fprintf(stderr, "Invalid file format\n");
    exit(1);
  }
  m /= 2;
}

/*
  @type specifies the format of input file:
    0 - adjacency matrix
    1 - parovi susjednih vrhova
*/
Graph::Graph(std::string filename, int type) {

  std::ifstream f(filename);
  if (!f.is_open()) {
    fprintf(stderr, "Failed to open a file %s\n", filename);
    return;
  }
  if (type == 0)
    load_adjecency(f);
  
  f.close();
}

void Graph::print() {
  printf("%u %u\n", n, m);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      printf("%d ", edges[i][j]);
    printf("\n");
  }
}

void Graph::make_SDPA() {
  repr.setParameterType(SDPA::PARAMETER_DEFAULT);
  repr.inputConstraintNumber(n+m);
  repr.inputBlockNumber(2);
  // set first block size and type
  repr.inputBlockSize(1, n);
  repr.inputBlockType(1, SDPA::SDP);
  // set second block size and type
  repr.inputBlockSize(2, m+1);
  repr.inputBlockType(2, SDPA::LP);

  repr.initializeUpperTriangleSpace();
  for(int i=0; i < n; i++) 
    repr.inputCVec(i+1, 1.0);

  repr.inputElement(0, 2, 1, 1, 1.0);
  
  int i;
  // first n contraints (for each vertix)
  for(i = 0; i < n; i++) {
    repr.inputElement(i+1, 1, i, i, 1.0);
  }
  // next m contraints (for each edge)
  for(int j = 0; j < n; j++) {
    for(int k = 0; k < edges[j].size(); k++) {
      repr.inputElement(i+1, 1, j, edges[j][k], 1.0);
      repr.inputElement(i+1, 1, edges[j][k], j, 1.0);
      repr.inputElement(i+1, 2, 1, 1, -2.0);
      repr.inputElement(i+1, 2, (i-n)+2, (i-n)+2, 1.0);
      i++;
    }
  }
  repr.initializeUpperTriangle();
  repr.initializeSolve();

  repr.solve();
  double* Y = repr.getResultYMat(1);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      printf("%d ", Y[i*n+j]);
    }
    printf("\n");
  }
}