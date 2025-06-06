#include <iostream>
#include <fstream>
#include <format>
#include <random>

void generate_graph(int n, double p) {
  std::ofstream file;
  file.open(std::format("examples/random_prob/graph{}_{:.2f}", n, p));
  file << n << '\n';
  int matrix[n][n];
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      if((double)rand() / RAND_MAX <= p) {
        matrix[i][j] = 1;
        matrix[j][i] = 1;
      } else {
        matrix[i][j] = 0;
        matrix[j][i] = 0;
      }
    }
  }
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      file << matrix[i][j];
      if(j == n-1)
        file << '\n';
      else
        file << ' ';
    }
  }
  file.flush();
  file.close();
}

int main() {
  for(int n = 5; n < 1000; n+=5) {
    for(double p = 0; p < 1; p+=0.05) {
      generate_graph(n, p);
    }
  }
}