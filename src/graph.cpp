#include "graph.h"
#include "openblas64/lapack.h"

#define max(x, y) ((x) > (y) ? (x) : (y))

unsigned uint_from_str(std::string val, int index, size_t* next_index) {
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
    size_t index = 0;
    n = uint_from_str(line, index, &index);
  }
  unsigned i = 0;
  while (std::getline(fs, line) && i < n) {
    if(line.length() < 2*n - 1) {
      fprintf(stderr, "Invalid file format line %d too short %ld %d\n", i, line.length(), 2*n);
      exit(1);
    }
    double r = 0.0;
    for(unsigned j = 0; j < n; j++) {
      if(line[j*2] == '1') {
        r++;
        if(j > i) {
          edges[i].push_back(j);
          m++;
        }
      }
    }
    maxr = max(maxr, r);
    i++;
  }
  if (i < n) {
    fprintf(stderr, "Invalid file format\n");
    exit(1);
  }
}

static void unpack_matrix(double * X, int n) {
  for(int i = n-1; i >= 0; i--) {
    for(int j = i; j >= 0; j--) {
      X[i*n + j] = X[i*(i+1)/2+j];
    }  
  }
  for(int i = 0; i < n; i++) {
    for(int j = i+1; j < n; j++) {
      X[i*n+j] = X[j*n+i];
    }
  }
}

static void swap_row(double* A, int i, int j, int n) {
  for(int k = 0; k < n; k++) {
    double x = A[i*n+k];
    A[i*n+k] = A[j*n+k];
    A[j*n+k] = x;
  }
}

static void swap_column(double* A, int i, int j, int n) {
  for(int k = 0; k < n; k++) {
    double x = A[k*n+i];
    A[k*n+i] = A[k*n+j];
    A[k*n+j] = x;
  }
}


// algorithm from https://eprints.maths.manchester.ac.uk/1199/1/covered/MIMS_ep2008_116.pdf
static void cholesky(double* A, int n) {
  int P[n];
  for(int i = 0; i < n; i++)
    P[i] = i;

  for(int k = 0; k < n; k++) {
    int max_i; double max_val=-100;
    for(int i = k; i < n; i++) {
      if(A[i*n+i] > max_val) {
        max_i = i; max_val = A[i*n+i];
      }
    }
    if(max_i != k) {
      swap_row(A, k, max_i, n);
      swap_column(A, k, max_i, n);
      int tmp = P[max_i];
      P[max_i] = P[k];
      P[k] = tmp;
    }

    printf("%lf\n", A[k*n+k]);
    A[k*n+k] = sqrt(A[k*n+k]);

    for(int j = k+1; j < n; j++)
      A[j*n+k] = A[j*n+k] / A[k*n+k];

    for(int j = k+1; j < n; j++) {
      for(int i = k+1; i <= j; i++) {
        A[j*n+i] -= A[i*n+k]*A[j*n+k];
      }
    }
  }

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < i; j++)
      A[j*n+i] = 0;
  }

  for(int i = 0; i < n; i++) {
    while (P[i] != i) {
      swap_column(A, i, P[i], n);
      int temp = P[i];
      P[i] = P[temp];
      P[temp] = temp;
    }
  }
}

/*
  @type specifies the format of input file:
    0 - adjacency matrix
    1 - parovi susjednih vrhova
*/
Graph::Graph(std::string filename, int type) : gen(rd()), ndist(0.0, 1.0) {

  std::ifstream f(filename);
  if (!f.is_open()) {
    fprintf(stderr, "Failed to open a file %s\n", filename.data());
    return;
  }
  if (type == 0)
    load_adjecency(f);
  
  f.close();

  vecColoring = (double*) malloc(max(n*n, (n+m+1)*(n+m+2)/2)*sizeof(double));
  
  values[0][0] = 1.0;
  values[1][0] = 1.0;
  values[1][1] = 2.0;
  values[1][2] = 1.0;
}

void Graph::print() {
  printf("%u %u\n", n, m);
  for(unsigned i = 0; i < n; i++) {
    for(unsigned j = 0; j < edges[i].size(); j++)
      printf("%d %d\n", i, edges[i][j]);
  }
}


void Graph::make_SDP() {
  DSDPCreate(n+m, &solver);
  for(unsigned i = 1; i <= n; i++)
    DSDPSetDualObjective(solver, i, 1.0);

  // creating cone with all constraints and C matrix
  DSDPCreateSDPCone(solver, 1, &cone);

  // setting the C matrix
  // n*(n+1) / 2 - 1 predstavnja n-ti elementi dijagonale
  indices[0][0] = (n+1)*(n+2)/2-1;
  SDPConeSetASparseVecMat(cone, 0, 0, n+m+1,
    -1.0, 0,
    indices[0], values[0], 1);
  SDPConeViewDataMatrix(cone, 0, 0);

  unsigned i = 1;
  for(; i <= n; i++) {
    // setting the matrix for first n contraints, relating to each vertix
    printf("%d\n", i);
    indices[i][0] = i*(i+1)/2-1;

    SDPConeSetASparseVecMat(cone, 0, i, n+m+1,
      1.0, 0,
      indices[i], values[0], 1);
    SDPConeViewDataMatrix(cone, 0, i);
  }

  for(unsigned j = 0; j < n; j++) {
    for(size_t k = 0; k < edges[j].size(); k++) {
      printf("%d\n", i);
      // setting the matrix for next m contraints, relating to each edge
      indices[i][0] = (edges[j][k])*(edges[j][k]+1)/2 + j;
      indices[i][1] = (n+1)*(n+2)/2-1;
      indices[i][2] = (i+1)*(i+2)/2-1;
      
      SDPConeSetASparseVecMat(cone, 0, i, n+m+1,
        1.0, 0,
        indices[i], values[1], 3);
      SDPConeViewDataMatrix(cone, 0, i);
      i++;
    }
  }

  DSDPSetup(solver);
  DSDPSolve(solver);
  DSDPComputeX(solver);
  int len;
  SDPConeGetXArray(cone, 0, &vecColoring, &len);
  double x;
  DSDPGetDObjective(solver, &x);
  // x = -1/(k-1)
  k = -1/x + 1;
  printf("objective: %lf, k: %lf\n", x, k);

  unpack_matrix(vecColoring, n);
  for(unsigned i = 0; i < n; i++) {
    for(unsigned j = 0; j < n; j++) {
      printf("%lf ", vecColoring[i*n+j]);
    }
    printf("\n");
  }
  cholesky(vecColoring, n);
}

void generate_rand_uvector(double* r, int n, std::normal_distribution<double>& nd, std::mt19937& gen) {
  for(int i = 0; i < n; i++)
    r[i] = nd(gen);

  double len = 0;
  for(int i = 0; i < n; i++)
    len += r[i] * r[i];
  for(int i = 0; i < n; i++)
    r[i] /= sqrt(len);
}

double scalar(double* x, double* y, int n) {
  double prod = 0.;
  for(int i = 0; i < n; i++) {
    prod += x[i] * y[i];
  }
  return prod;
}

void Graph::color() {
  unsigned colored = 0;
  memset(colors, 0, n*sizeof(int));
  
  make_SDP();

  int it = 1, nc = 1;
  double c = 0.5;//sqrt(2 * (k-2) / k * log(maxr));
  printf("%lf %lf %lf\n", c, k, maxr);
  while (colored != n) {
    bool used = 0;
    double r[n];
    generate_rand_uvector(r, n, ndist, gen);
  
    for(unsigned i = 0; i < n; i++) {
      if(colors[i]) continue;

      if(scalar(r, vecColoring + i*n, n) > c) {
        used = 1;
        colors[i] = nc;
        colored++;
      }
    }

    for(unsigned i = 0; i < n; i++) {
      if(colors[i] != it) continue;
      for(size_t j = 0; j < edges[i].size(); j++) {
        if(colors[edges[i][j]] == it) {
          colors[edges[i][j]] = 0;
          colored--;
        }
      }
    }
    it++;
    if(used) nc++;
  }
  
  printf("Graf uspjesno obojan s %d boja\n", nc-1);
  for(unsigned i = 0; i < n; i++) {
    printf("%d ", colors[i]);
  }
}