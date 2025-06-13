#include "graph.h"
#include <raylib.h>
#include <format>

Color boje[11] = {
  BLACK,
  RED,
  GREEN,
  BLUE,
  ORANGE,
  YELLOW,
  DARKBLUE,
  DARKPURPLE,
  VIOLET,
  SKYBLUE,
  BROWN
};

void set_graph_pos(Graph* graph, Vector2* pos) {
  for(int i = 0; i < graph->n; i++) {
    pos[i].x = 500 + 300 * cos(2*PI/graph->n * i); 
    pos[i].y = 400 + 300 * sin(2*PI/graph->n * i); 
  }
}

void draw_graph_colored(Graph* graph, Vector2* pos, int* colors) {
  BeginDrawing();
  for(int i = 0; i < graph->n; i++) {
    for(int j = 0; j < graph->edges[i].size(); j++) {
      int k = graph->edges[i][j];
      DrawLineEx(pos[i], pos[k], 3.0, BLACK);
    }
  }

  for(int i = 0; i < graph->n; i++) {
    DrawCircle(pos[i].x, pos[i].y, 30, boje[colors[i]]);
    char index[2];
    index[0] = i + '0'; index[1] = '\0';
  }
  
  for(int i = 0; i < graph->n; i++) {
    char index[2];
    index[0] = i + '0'; index[1] = '\0';
    DrawText(index, pos[i].x, pos[i].y, 20, WHITE);
  }
  EndDrawing();
}

void draw_graph_clique(Graph* graph, Vector2* pos, int* clique, int size) {
  BeginDrawing();
  for(int i = 0; i < graph->n; i++) {
    for(int j = 0; j < graph->edges[i].size(); j++) {
      int k = graph->edges[i][j];
      DrawLineEx(pos[i], pos[k], 3.0, BLACK);
    }
  }

  for(int i = 0; i < graph->n; i++)
    DrawCircle(pos[i].x, pos[i].y, 30, BLACK);
  
  for(int i = 0; i < size; i++) {
    for(int j = i; j < size; j++) {
      DrawLineEx(pos[clique[i]], pos[clique[j]], 3.0, RED);
    }
  }

  for(int i = 0; i < size; i++)
    DrawCircle(pos[clique[i]].x, pos[clique[i]].y, 30, RED);

  for(int i = 0; i < graph->n; i++) {
    char index[2];
    index[0] = i + '0'; index[1] = '\0';
    DrawText(index, pos[i].x, pos[i].y, 20, WHITE);
  }
  EndDrawing();
}

void write_help_message() {
  printf("usage: visualize.out option file");
  printf("options:");
  printf("--help         display help");
  printf("-delta_color   color the graph with delta+1 algorithm");
  printf("-color         color the graph with SDP algorithm");
  printf("-greedy_clique find clique using greedy algorithm");
  printf("-clique        find clique using SDP algorithm");
  exit(0);
}

void process_args(int argv, char** argc, int* flags) {
  for(int i = 1; i < argv-1; i++) {
    if(strcmp(argc[i], "--help") == 0)
      write_help_message();
    if(strcmp(argc[i], "-delta_color") == 0)
      flags[0] = 1;
    if(strcmp(argc[i], "-color") == 0)
      flags[1] = 1;
    if(strcmp(argc[i], "-greedy_clique") == 0)
      flags[2] = 1;
    if(strcmp(argc[i], "-clique") == 0)
      flags[3] = 1;
  }
  if(argv != 3) {
    printf("Invalid input\n");
    write_help_message();
  }
}

// usage: visualize.out option file
// options:
// --help         display help
// -delta_color   color the graph with delta+1 algorithm
// -color         color the graph with SDP algorithm
// -greedy_clique find clique using greedy algorithm
// -clique        find clique using SDP algorithm
int main(int argv, char** argc) {
  int flags[4];
  process_args(argv, argc, flags);
  
  Graph a(argc[2], 0);
  Vector2 pos[a.n];
  set_graph_pos(&a, pos);

  int array[a.n], n;
  
  if(flags[0])
    n = a.greedy_color(array);
  else if(flags[1])
    n = a.color(array, 0.);
  else if(flags[2])
    n = a.greedy_clique(array);
  else if(flags[3])
    n = a.find_max_clique(array);
  else {
    fprintf(stderr, "Invalid input\n");
    write_help_message();
  }

  InitWindow(1000, 800, "Graphs");
  SetTargetFPS(60);
  while(WindowShouldClose() == false) {
    ClearBackground({255, 255,255});
    if(flags[0] || flags[1])
      draw_graph_colored(&a, pos, array);
    if(flags[2] || flags[3])
      draw_graph_clique(&a, pos, array, n);
  }
  CloseWindow();
  return 0;
}