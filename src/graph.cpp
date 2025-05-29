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
      } else if(j > i) {
        not_edges[i].push_back(j);
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

static void convert_lapack_out(int n, double* eigval, double* vecColoring) {
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      vecColoring[i*n+j] *= sqrt(eigval[i]);
      if(i==j) continue;
      vecColoring[j*n+i] *= sqrt(eigval[j]);
      double temp = vecColoring[i*n+j];
      vecColoring[i*n+j] = vecColoring[j*n+i];
      vecColoring[j*n+i] = temp;
    }
  }
}

void Graph::vector_color() {
  DSDP solver;
  SDPCone cone;
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

  int lda = n;
  double eigval[n];
  double work[3*n-1];
  int lwork = 3*n-1, info;
  char v = 'V', u = 'U';
  printf("\n");
  LAPACK_dsyev(&v, &u, &lda, vecColoring, &lda, eigval, work, &lwork, &info);

  for(unsigned i = 0; i < n; i++)
    printf("%lf ", eigval[i]);
  printf("\n");

  convert_lapack_out(n, eigval, vecColoring);

  for(unsigned i = 0; i < n; i++) {
    for(unsigned j = 0; j < n; j++) {
      printf("%lf ", vecColoring[i*n+j]);
    }
    printf("\n");
  }

  DSDPTerminationReason reason;
  DSDPStopReason(solver,&reason); 
  DSDPSolutionType pdfeasible;
  DSDPGetSolutionType(solver, &pdfeasible);
  if(pdfeasible == DSDP_PDFEASIBLE)
    printf("feasible\n");
  if(reason == DSDP_CONVERGED)
    printf("converged\n");
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
  
  vector_color();

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
  printf("\n");
}

void Graph::greedy_color() {
  memset(colors, 0, n*sizeof(int));
  for(int i = n-1; i >= 0; i--) {
    bool taken[int(maxr)+1];
    memset(taken, 0, (maxr+1)*sizeof(bool));
    for(size_t j = 0; j < edges[i].size(); j++) {
      int l = edges[i][j];
      if(colors[l] != 0)
        taken[colors[l]-1] = 1;
    }

    for(int j = 0; j < maxr+1; j++) {
      if(!taken[j]) {
        colors[i] = j+1;
        break;
      }
    }
  }

  for(unsigned i = 0; i < n; i++) {
    printf("%d ", colors[i]);
  }
  printf("\n");
}

void Graph::find_max_clique() {
  DSDP solver;
  SDPCone cone;
  // inverse m
  int m_ = (n-1)*n/2 - m;
  DSDPCreate((n+1)+m_, &solver);
  for(unsigned i = 1; i <= (n+1)+m_; i++)
    DSDPSetDualObjective(solver, i, 1.0);

  // creating cone with all constraints and C matrix
  DSDPCreateSDPCone(solver, 1, &cone);

  // setting the C matrix
  double Cval[2*n];
  int Cind[2*n];
  for(unsigned i = 0; i < n; i++) {
    Cval[i] = 0.5;
    Cind[i] = (i+1)*(i+2)/2-1;
  }
  for(unsigned i = 0; i < n; i++) {
    Cval[i+n] = 0.25;
    Cind[i+n] = n*(n+1)/2+i;
  }
  SDPConeSetASparseVecMat(cone, 0, 0, n+1,
    1.0, 0,
    Cind, Cval, 2*n);
  SDPConeViewDataMatrix(cone, 0, 0);

  unsigned i = 1;
  double Aval[1] = {1.};
  int Aind[n+2][1];
  for(; i <= n+1; i++) {
    // setting the matrix for first n contraints, relating to each vertix
    printf("%d\n", i);
    Aind[i][0] = i*(i+1)/2-1;

    SDPConeSetASparseVecMat(cone, 0, i, n+1,
      1.0, 0,
      Aind[i], Aval, 1);
    SDPConeViewDataMatrix(cone, 0, i);
  }

  double Aval2[6] = {1., 1., 1., 1., 1., 1.};
  int Aind2[m_][6];
  for(unsigned j = 0; j < n; j++) {
    for(size_t k = 0; k < not_edges[j].size(); k++) {
      // setting the matrix for next m_ contraints, relating to each non-edge
      // (j, l) is a non-edge -- j is the smaller number
      int l = not_edges[j][k];
      printf("%d: (%d, %d)\n", i, j, l);
      Aind2[i][0] = (j+1)*(j+2)/2-1;
      Aind2[i][1] = l*(l+1)/2+j;
      Aind2[i][2] = (l+1)*(l+2)/2-1;
      Aind2[i][3] = n*(n+1)/2+j;
      Aind2[i][4] = n*(n+1)/2+l;
      Aind2[i][5] = (n+1)*(n+2)/2-1;
      
      SDPConeSetASparseVecMat(cone, 0, i, n+1,
        1.0, 0,
        Aind2[i], Aval2, 6);
      SDPConeViewDataMatrix(cone, 0, i);
      i++;
    }
  }

  DSDPSetup(solver);
  DSDPSolve(solver);
  DSDPComputeX(solver);
  int len = (n+1)*(n+1);
  double* X = (double*) malloc(len * sizeof(double));
  SDPConeGetXArray(cone, 0, &X, &len);
  double x;
  DSDPGetDObjective(solver, &x);
  printf("objective: %lf", x);

  unpack_matrix(X, n+1);

  int lda = n+1;
  double eigval[n+1];
  double work[3*(n+1)-1];
  int lwork = 3*(n+1)-1, info;
  char v = 'V', u = 'U';
  LAPACK_dsyev(&v, &u, &lda, X, &lda, eigval, work, &lwork, &info);

  convert_lapack_out(n+1, eigval, X);

  int vsign[n+1];
  double uvec[n+1];
  generate_rand_uvector(uvec, n+1, ndist, gen);
  for(unsigned i = 0; i < n+1; i++) {
    double res = scalar(X+i*(n+1), uvec, n+1);
    if(res >= 0) vsign[i] = 1;
    else vsign[i] = -1;
  }

  for(unsigned i = 0; i < n; i++) {
    for(size_t j = 0; j < not_edges[i].size(); j++) {
      int l = not_edges[i][j];
      if(abs(vsign[i] + vsign[l] + vsign[n])) {
        if(abs(vsign[i]-vsign[n]) > abs(vsign[l]-vsign[n]))
          vsign[i] = -vsign[i];
        else
          vsign[l] = -vsign[l];
      }
    }
  }

  printf("Maximum clique is:\n");
  for(unsigned i = 0; i < n; i++) {
    if(vsign[i] == vsign[n])
      printf("%d ", i);
  }
}