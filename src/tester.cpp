#include "graph.h"
#include <format>

#define MIN(x, y) ((x) > (y) ? (y) : (x))
#define MAX(x, y) ((x) < (y) ? (y) : (x))

int colors_greedy[19][20];
int colors_sdp1[19][20];
int colors_sdp2[19][20];
int colors_sdp3[19][20];
int clique_greedy[19][20];
int clique_sdp[19][20];
double delta[19][20];

void test_graph(int n, double p) {
  Graph graph(std::format("examples/random_prob/graph{}_{:.2f}", n, p), 0);
  int ind_n = n/5 -1;
  int ind_p = p/0.05;
  int array[n];
  colors_greedy[ind_n][ind_p] = graph.greedy_color(array);
  int min1 = n, min2 = n, min3 = n;
  for(int i = 0; i < 3; i++) {
    min1 = MIN(min1, graph.color(array, -0.2));
    min2 = MIN(min2, graph.color(array, 0.));
    min3 = MIN(min3, graph.color(array, 0.2));
  }
  colors_sdp1[ind_n][ind_p] = min1;
  colors_sdp2[ind_n][ind_p] = min2;
  colors_sdp3[ind_n][ind_p] = min3;
  clique_greedy[ind_n][ind_p] = graph.greedy_clique(array);
  min1 = 0;
  for(int i = 0; i < 3; i++) {
    min1 = MAX(min1, graph.find_max_clique(array));
  }
  clique_sdp[ind_n][ind_p] = min1;
  delta[ind_n][ind_p] = graph.get_delta();
}

int main() {
  for(int n = 5; n < 50; n+=5) {
    for(double p = 0; p < 1; p+=0.05) {
      test_graph(n, p);
    }
  }

  std::ofstream file;
  file.open("examples/results_by_n");
  
  int sum1[20] = {0}, sum2[20] = {0}, sum3[20] = {0}, sum4[20] = {0}, sum5[20] = {0}, sum6[20] = {0};
  double sum7[20] = {0};
  memset(sum1, 0, 20*sizeof(int));
  memset(sum2, 0, 20*sizeof(int));
  memset(sum3, 0, 20*sizeof(int));
  memset(sum4, 0, 20*sizeof(int));
  memset(sum5, 0, 20*sizeof(int));
  memset(sum6, 0, 20*sizeof(int));
  memset(sum7, 0, 20*sizeof(int));
  for(int n = 5; n < 50; n+=5) {
    int ind_n = n/5 -1;
    for(double p = 0; p < 1; p+=0.05) {
      int ind_p = p/0.05;
      sum1[ind_n] += colors_greedy[ind_n][ind_p];
      sum2[ind_n] += colors_sdp1[ind_n][ind_p];
      sum3[ind_n] += colors_sdp2[ind_n][ind_p];
      sum4[ind_n] += colors_sdp3[ind_n][ind_p];
      sum5[ind_n] += clique_greedy[ind_n][ind_p];
      sum6[ind_n] += clique_sdp[ind_n][ind_p];
      sum7[ind_n] += delta[ind_n][ind_p];
     }
    file << n << " " << sum1[ind_n] / 20.  << " " << sum2[ind_n] / 20. << " "
      << sum3[ind_n] / 20.  << " " << sum4[ind_n] / 20. << " "
      << sum5[ind_n] / 20.  << " " << sum6[ind_n] / 20. << " "
      << sum7[ind_n] / 20.  << std::endl;
  }
  file.close();

  {
    std::ofstream file;
    file.open("examples/results_by_p");
    
    memset(sum1, 0, 20*sizeof(int));
    memset(sum2, 0, 20*sizeof(int));
    memset(sum3, 0, 20*sizeof(int));
    memset(sum4, 0, 20*sizeof(int));
    memset(sum5, 0, 20*sizeof(int));
    memset(sum6, 0, 20*sizeof(int));
    memset(sum7, 0, 20*sizeof(double));
    for(double p = 0; p < 1; p+=0.05) {
      int ind_p = p/0.05;
      for(int n = 5; n < 50; n+=5) {
        int ind_n = n/5-1;
        sum1[ind_p] += colors_greedy[ind_n][ind_p];
        sum2[ind_p] += colors_sdp1[ind_n][ind_p];
        sum3[ind_p] += colors_sdp2[ind_n][ind_p];
        sum4[ind_p] += colors_sdp3[ind_n][ind_p];
        sum5[ind_p] += clique_greedy[ind_n][ind_p];
        sum6[ind_p] += clique_sdp[ind_n][ind_p];
        sum7[ind_p] += delta[ind_n][ind_p];
      }
      file << p << " " << sum1[ind_p] / 9.  << " " << sum2[ind_p] / 9. << " "
        << sum3[ind_p] / 9.  << " " << sum4[ind_p] / 9. << " "
        << sum5[ind_p] / 9.  << " " << sum6[ind_p] / 9. << " "
        << sum7[ind_p] / 9.  << std::endl;
    }
    file.close();
  }
  return 0;
}