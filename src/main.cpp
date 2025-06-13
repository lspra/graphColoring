#include "graph.h"

void write_help_message() {
  printf("usage: main.out [options] file\n");
  printf("options:\n");
  printf("--help         display help\n");
  printf("-print         print graph\n");
  printf("-delta_color   color the graph with delta+1 algorithm\n");
  printf("-color         color the graph with SDP algorithm\n");
  printf("-greedy_clique find clique using greedy algorithm\n");
  printf("-clique        find clique using SDP algorithm\n");
  exit(0);
}

void process_args(int argv, char** argc, int* flags) {
  int count = 0;
  for(int i = 1; i < argv-1; i++) {
    if(strcmp(argc[i], "--help") == 0)
      write_help_message();
    if(strcmp(argc[i], "-print") == 0) {
      flags[0] = 1;
      count++;
    }
    if(strcmp(argc[i], "-delta_color") == 0) {
      flags[1] = 1;
      count++;
    }
    if(strcmp(argc[i], "-color") == 0) {
      flags[2] = 1;
      count++;
    }
    if(strcmp(argc[i], "-greedy_clique") == 0) {
      flags[3] = 1;
      count++;
    }
    if(strcmp(argc[i], "-clique") == 0) {
      flags[4] = 1;
      count++;
    }
  }
  if(argv <= 1 || count != argv-2) {
    printf("Invalid input\n");
    write_help_message();
  }
}

// usage: main.out [options] file
// options:
// --help         display help
// -print         print graph
// -delta_color   color the graph with delta+1 algorithm
// -color         color the graph with SDP algorithm
// -greedy_clique find clique using greedy algorithm
// -clique        find clique using SDP algorithm
int main(int argv, char** argc) {
  int flags[5] = {0};
  process_args(argv, argc, flags);

  Graph g(argc[argv-1], 0);
  int array[g.n];
  if(flags[0])
    g.print();
  if(flags[1]) {
    int result = g.greedy_color(array);
    if(result == -1)
      fprintf(stderr, "error while (delta+1)-coloring graph\n");
    else {
      printf("(delta+1)-colored graph with %d colors\n Colors are:\n", result);
      for(int i = 0; i < g.n; i++)
        printf("%d ", array[i]);
      printf("\n");
    }
  }
  if(flags[2]) {
    int result = g.color(array, 0.);
    if(result == -1)
      fprintf(stderr, "error while SDP-coloring graph\n");
    else {
      printf("SDP-colored graph with %d colors\n Colors are:\n", result);
      for(int i = 0; i < g.n; i++)
        printf("%d ", array[i]);
      printf("\n");
    }
  }
  if(flags[3]) {
    int result = g.greedy_clique(array);
    if(result == -1)
      fprintf(stderr, "error while finding clique using greedy algorithm\n");
    else {
      printf("Greedy algorithm found clique of size %d\n Vertices in clique are:\n", result);
      for(int i = 0; i < result; i++)
        printf("%d ", array[i]);
      printf("\n");
    }
  }
  if(flags[4]) {
    int result = g.find_max_clique(array);
    if(result == -1)
      fprintf(stderr, "error while finding clique using SDP algorithm\n");
    else {
      printf("SDP algorithm found clique of size %d\n Vertices in clique are:\n", result);
      for(int i = 0; i < result; i++)
        printf("%d ", array[i]);
      printf("\n");
    }
  }
  return 0;
}