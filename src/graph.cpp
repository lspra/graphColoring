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
  m=m_=0;
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
        m_++;
      }
    }
    if(r > maxr) {
      maxr = r;
      maxr_ind = i;
    }
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
Graph::Graph(std::string filename, int type) : rd{}, gen{rd()}, ndist{0.0, 1.0} {
  std::ifstream f(filename);
  if (!f.is_open()) {
    fprintf(stderr, "Failed to open a file %s\n", filename.data());
    return;
  }
  if (type == 0)
    load_adjecency(f);
  
  f.close();

  vecColoring = (double*) malloc(max(n*n, (n+1)*(n+2))*sizeof(double));
  
  values[0][0] = 1.0;
  values[1][0] = 1.0;
  values[1][1] = 2.0;
  values[1][2] = 1.0;
}

Graph::~Graph() {
  free(vecColoring);
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
  DSDPSetGapTolerance(solver, 1e-2);
  for(unsigned i = 1; i <= n; i++)
    DSDPSetDualObjective(solver, i, 1.0);

  // creating cone with all constraints and C matrix
  DSDPCreateSDPCone(solver, 1, &cone);
  SDPConeSetXArray(cone, 0, n+1, vecColoring, (n+1)*(n+2));

  // setting the C matrix
  // n*(n+1) / 2 - 1 predstavnja n-ti elementi dijagonale
  indices[0][0] = (n+1)*(n+2)/2-1;
  SDPConeSetASparseVecMat(cone, 0, 0, n+1,
    -1.0, 0,
    indices[0], values[0], 1);

  unsigned i = 1;
  for(; i <= n; i++) {
    // setting the matrix for first n contraints, relating to each vertix
    //printf("%d\n", i);
    indices[i][0] = i*(i+1)/2-1;

    SDPConeSetASparseVecMat(cone, 0, i, n+1,
      1.0, 0,
      indices[i], values[0], 1);
  }

  for(unsigned j = 0; j < n; j++) {
    for(size_t k = 0; k < edges[j].size(); k++) {
      //// printf("%d\n", i);
      // setting the matrix for next m contraints, relating to each edge
      indices[i][0] = (edges[j][k])*(edges[j][k]+1)/2 + j;
      indices[i][1] = (n+1)*(n+2)/2-1;
      indices[i][2] = (i+1)*(i+2)/2-1;
      
      SDPConeSetASparseVecMat(cone, 0, i, n+1,
        1.0, 0,
        indices[i], values[1], 2);
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
  vecColor = 1;
}

double Graph::lovasz() {
  DSDP solver;
  SDPCone cone;
  DSDPCreate(n+m, &solver);
  DSDPSetDualObjective(solver, 1, 1.0);

  for(unsigned i = 2; i < 1+m; i++)
    DSDPSetDualObjective(solver, i, 0.0);
  
  // creating cone with all constraints and C matrix
  DSDPCreateSDPCone(solver, 1, &cone);

  // setting the C matrix
  double Cval[(n+1)*n/2];
  int Cind[(n+1)*n/2];
  for(unsigned i = 0; i < (n+1)*n/2; i++) {
    Cval[i] = 1.;
    Cind[i] = i;
  }
  SDPConeSetASparseVecMat(cone, 0, 0, n,
    -1.0, 0,
    Cind, Cval, 2*n);
  //SDPConeViewDataMatrix(cone, 0, 0);

  double Aval[n];
  int Aind[n];
  for(unsigned i = 0; i < n; i++) {
    // setting the matrix for first n contraints, relating to each vertix
    //printf("%d\n", i);
    Aval[i] = 1.;
    Aind[i] = (i+1)*(i+2)/2-1;
  }
  SDPConeSetASparseVecMat(cone, 0, 1, n,
    1.0, 0,
    Aind, Aval, n);

  double Aval2 = 0.;
  int Aind2[m];
  int i = 2;
  for(unsigned j = 0; j < n; j++) {
    for(size_t k = 0; k < edges[j].size(); k++) {
      // setting the matrix for next m contraints, relating to each edge
      // (j, l) is a edge -- j is the smaller number
      Aind2[i-2] = (edges[j][k])*(edges[j][k]+1)/2 + j;
      
      SDPConeSetASparseVecMat(cone, 0, i, n,
        1.0, 0,
        &Aind2[i-2], &Aval2, 1);
      i++;
    }
  }

  DSDPSetup(solver);
  DSDPSolve(solver);
  double x;
  DSDPGetDObjective(solver, &x);
  return -x;
}


void generate_rand_uvector(double* r, int n, std::normal_distribution<double>& nd, std::mt19937& gen) {
  for(int i = 0; i < n; i++)
    r[i] = nd(gen);
}

void normalize(double* r, int n) {
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

void Graph::color_init(int* color) {
  memset(color, 0, n*sizeof(int));
  if(!vecColor)
    vector_color();
}

unsigned Graph::next_color(int* color, int nc, unsigned colored, double delta) {
  if (colored != n) {
    bool used = 0;
    while (!used) {
      double r[n];
      generate_rand_uvector(r, n, ndist, gen);
      normalize(r,n);
      printf("\nrandom vector %d is:\n", nc);
      for(unsigned i = 0; i < n; i++)
        printf("%.3lf ", r[i]);
      printf("\n");
      for(unsigned i = 0; i < n; i++) {
        if(color[i]) continue;

        printf("%.3lf ", scalar(r, vecColoring + i*n, n) );
        if(scalar(r, vecColoring + i*n, n) >= delta) {
          used = 1;
          color[i] = nc;
          colored++;
        }
      }

      for(unsigned i = 0; i < n; i++) {
        if(color[i] != nc) continue;
        for(size_t j = 0; j < edges[i].size(); j++) {
          if(color[edges[i][j]] == nc) {
            if(scalar(r, vecColoring + i*n, n) > scalar(r, vecColoring + edges[i][j]*n, n)) {
              color[edges[i][j]] = 0;
              colored--;
            } else {
              color[i] = 0;
              colored--;
              break;
            }
          }
        }
      }
    }    
  }
  return colored;
}

int Graph::color(int* color, double delta) {
  color_init(color);

  unsigned colored = 0;

  int nc = 1;
  while (colored != n) {
    colored = next_color(color, nc++, colored, delta);
  }
  
  printf("Graf uspjesno obojan s %d boja\n", nc-1);

  for(unsigned i = 0; i < n; i++) {
    printf("%d ", color[i]);
  }
  printf("\n");
  if(!check_colors_valid(color, n))
    return -1;
  return nc-1;
}

int Graph::greedy_color(int* color) {
  int highest = 0;
  memset(color, 0, n*sizeof(int));
  for(int i = n-1; i >= 0; i--) {
    bool taken[int(maxr)+1];
    memset(taken, 0, (maxr+1)*sizeof(bool));
    for(size_t j = 0; j < edges[i].size(); j++) {
      int l = edges[i][j];
      if(color[l] != 0)
        taken[color[l]-1] = 1;
    }

    for(int j = 0; j < maxr+1; j++) {
      if(!taken[j]) {
        color[i] = j+1;
        if(j+1 > highest) highest = j+1;
        break;
      }
    }
  }

  for(unsigned i = 0; i < n; i++) {
    printf("%d ", color[i]);
  }
  printf("\n");
  if(!check_colors_valid(color, n))
    return -1;
  return highest;
}

bool Graph::check_colors_valid(int* colors, int len) {
  for(int i = 0; i < len; i++){
    for(size_t j = 0; j < edges[i].size(); j++) {
      int l = edges[i][j];
      if(colors[i] == colors[l]) {
        printf("coloring is not valid!\n");
        return 0;
      }
    }
  } 
  printf("coloring is valid!\n");
  return 1;
}

bool Graph::check_clique_valid(int* clique, int len) {
  for(int i = 0; i < len; i++) {
    for(int j = i; j < len; j++) {
      for(unsigned k = 0; k < not_edges[clique[i]].size(); k++) {
        if(not_edges[clique[i]][k] == clique[j]) {
          printf("clique not valid, no edges %d %d\n", clique[i], clique[j]);
          return 0;
        }
      }
    }
  }
  printf("clique valid\n");
  return 1;
}

//return size of clique
int Graph::greedy_clique(int* clique) {
  int size = 1;
  clique[0] = 0;
  for(unsigned i = 1; i < n; i++) {
    bool valid = true;
    for(int j = 0; j < size; j++) {
      int index = clique[j];
      for(size_t k = 0; k < not_edges[index].size(); k++) {
        if(not_edges[index][k] == i) {
          valid = 0; break;
        }
      }
      if(!valid) break;
    }
    if(valid) {
      clique[size++] = i;
    }
  }
  check_clique_valid(clique, size);
  printf("found clique of size %d:\n", size);
  for(int i = 0; i < size; i++)
    printf("%d ", clique[i]);
  printf("\n");
  return size;
}

int Graph::find_max_clique(int* clique) {
  DSDP solver;
  SDPCone cone;
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
  // TODO ovo bi vjerojatno trebalo biti pomnoženo s -1 jer tražimo minimum a v def problema je maksimum 
  SDPConeSetASparseVecMat(cone, 0, 0, n+1,
    -1.0, 0,
    Cind, Cval, 2*n);
  SDPConeViewDataMatrix(cone, 0, 0);

  unsigned i = 1;
  double Aval[1] = {1.};
  int Aind[n+2][1];
  for(; i <= n+1; i++) {
    // setting the matrix for first n contraints, relating to each vertix
    //printf("%d\n", i);
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
      Aind2[i-n-2][0] = (j+1)*(j+2)/2-1;
      Aind2[i-n-2][1] = l*(l+1)/2+j;
      Aind2[i-n-2][2] = (l+1)*(l+2)/2-1;
      Aind2[i-n-2][3] = n*(n+1)/2+j;
      Aind2[i-n-2][4] = n*(n+1)/2+l;
      Aind2[i-n-2][5] = (n+1)*(n+2)/2-1;
      
      SDPConeSetASparseVecMat(cone, 0, i, n+1,
        1.0, 0,
        Aind2[i-n-2], Aval2, 6);
      SDPConeViewDataMatrix(cone, 0, i);
      i++;
    }
  }
  
  int len = (n+1)*(n+1);
  double* X = (double*) malloc(len * sizeof(double));
  SDPConeSetXArray(cone, 0, n+1, X, (n+1)*(n+1));
  DSDPSetup(solver);
  DSDPSolve(solver);
  DSDPComputeX(solver);
  SDPConeGetXArray(cone, 0, &X, &len);
  printf("%d\n", len);
  double x;
  DSDPGetDObjective(solver, &x);
  printf("objective: %lf\n", x);

  unpack_matrix(X, n+1);

  for(unsigned i = 0; i < n+1; i++) {
    for(unsigned j = 0; j < n+1; j++) {
      printf("%.3lf ", X[i*(n+1)+j]);
    }
    printf("\n");
  }
  printf("\n");

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
  normalize(uvec, n+1);
  for(unsigned i = 0; i < n+1; i++) {
    printf("%.3lf ", uvec[i]);
  }
  printf("\n");

  for(unsigned i = 0; i < n+1; i++) {
    double res = scalar(X+i*(n+1), uvec, n+1);
    if(res >= 0) vsign[i] = 1;
    else vsign[i] = -1;
  }

  for(int i = 0; i < n+1; i++) {
    printf("%d ", vsign[i]);
  }
  printf("\n");

  for(unsigned i = 0; i < n; i++) {
    for(size_t j = 0; j < not_edges[i].size(); j++) {
      int l = not_edges[i][j];
      if(abs(vsign[i] + vsign[l] + vsign[n]) != 1) {
        if(abs(vsign[i]-vsign[n]) > abs(vsign[l]-vsign[n]))
          vsign[i] = -vsign[i];
        else
          vsign[l] = -vsign[l];
      }
    }
  }

  int count = 0;

  printf("Maximum clique is:\n");
  for(unsigned i = 0; i < n; i++) {
    if(vsign[i] == vsign[n]) {
      clique[count++] = i;
      printf("%d ", i);
    }
  }
  printf("\nSize: %d\n", count);
  free(X);
  if(!check_clique_valid(clique, count))
    return -1;
  return count;
}

double Graph::get_delta() {
  return maxr;
}