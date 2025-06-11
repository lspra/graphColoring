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
  static int nc = 1, colored = 0;
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
  if(colored < graph->n) {
    colored = graph->next_color(colors, nc++, colored);
  } else
    nc++;
  Image im = LoadImageFromScreen();
  ExportImage(im, std::format("../docs/images/colors_it{}.png", nc).data());
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
  Image im = LoadImageFromScreen();
  ExportImage(im, std::format("../docs/images/max_clique.png").data());
}

int main() {
  Graph a("examples/random_prob/graph10_0.50", 0);
  Vector2 pos[a.n];
  int clique[a.n], n;
  n = a.find_max_clique(clique);
  set_graph_pos(&a, pos);
  InitWindow(1000, 800, "Graphs");
  SetTargetFPS(60);
  while(WindowShouldClose() == false) {
    ClearBackground({255, 255,255});
    draw_graph_clique(&a, pos, clique, n);
  }
  CloseWindow();
  return 0;
}