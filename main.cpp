#include <allegro.h>
#include <fstream>
#include <math.h>
using namespace std;
#define screen_w 1300
#define screen_h 800
#define PI 3.14


void read_file();
void centre();
void roberts();
void lighting();
void draw();
void move_up(float W);
void move_down(float S);
void move_left(float A);
void move_right(float D);
void rotate_x(float K);
void rotate_y(float M);
void rotate_z(float H);
void approximation(float P);
void distance(float N);
void control();


int tops, facets;
int G[100][100];
float X[100], Y[100], Z[100];
int R[10];
float Center_X, Center_Y, Center_Z;
float Xc = 0, Yc = 0, Zc = 0;
float X1, Y1, Z1;
float A[10],B[10],C[10],D[10];
float V[4][10];
float plane;
float argument;
float E[4] = {0, 0, -1, 0};
float light[3] = {1, 1, 1};
float argument_light;
float correct_V[10];
float shine[10];


int main() {
  allegro_init();
  install_keyboard();
  set_color_depth(32);
  set_gfx_mode(GFX_AUTODETECT_WINDOWED, screen_w, screen_h, 0, 0);
  read_file();
  draw();
  control();
  allegro_exit();

  return 0;
}


void read_file() {
  ifstream file;
  float F_x, F_y, F_z, F_p;
  file.open("6polygon.txt");
  file >> tops >> facets;
  for (int i = 0; i < tops; i++) {
    file >> F_x;
    X[i] = F_x;
  }
  for (int i = 0; i < tops; i++) {
    file >> F_y;
    Y[i] = F_y;
  }
  for (int i = 0; i < tops; i++) {
    file >> F_z;
    Z[i] = F_z;
  }
  for (int i = 0; i < facets; i++) {
    for (int j = 0; j < 7; j++) {
      file >> F_p;
      G[i][j] = F_p;
      int k;
      if (j == 0) {
        k = F_p;
      }
      R[i] = j + 1;
      if (j != 0 && F_p == k) {
        break;
      }
    }
  }
  file.close();
}
void centre() {
  Center_X = 0, Center_Y = 0, Center_Z = 0;
  for (int j = 0; j < tops; j++) {
    Center_X += X[j];
    Center_Y += Y[j];
    Center_Z += Z[j];
  }

  Xc = Center_X / tops;
  Yc = Center_Y / tops;
  Zc = Center_Z / tops;
}

void roberts() {
  for (int i = 0; i < facets; i++) {
    float a1 = X[G[i][1]] - X[G[i][0]];
    float b1 = Y[G[i][1]] - Y[G[i][0]];
    float c1 = Z[G[i][1]] - Z[G[i][0]];
    float a2 = X[G[i][2]] - X[G[i][0]];
    float b2 = Y[G[i][2]] - Y[G[i][0]];
    float c2 = Z[G[i][2]] - Z[G[i][0]];
    A[i] = b2 * c1 - b1 * c2;
    B[i] = a1 * c2 - c1 * a2;
    C[i] = a2 * b1 - a1 * b2;
    D[i] = -(A[i] * X[G[i][0]] + B[i] * Y[G[i][0]] + C[i] * Z[G[i][0]]);
  }
  centre();
  for (int i = 0; i < facets; i++) {
    plane = A[i] * Xc + B[i] * Yc + C[i] * Zc + D[i];
    if (plane < 0) {
      A[i] = -A[i];
      B[i] = -B[i];
      C[i] = -C[i];
      D[i] = -D[i];
    }
  }
  for (int i = 0; i < facets; i++) {
    argument = sqrt(pow(A[i], 2) + pow(B[i], 2) + pow(D[i], 2));
    A[i] = A[i] / argument;
    B[i] = B[i] / argument;
    C[i] = C[i] / argument;
  }
  for (int i = 0; i < facets; i++) {
    V[0][i] = A[i];
    V[1][i] = B[i];
    V[2][i] = C[i];
    V[3][i] = D[i];
  }
  for (int i = 0; i < facets; i++) {
    correct_V[i] =
        V[0][i] * E[0] + V[1][i] * E[1] + V[2][i] * E[2] + V[3][i] * E[3];
  }
}
void lighting() {
  for (int i = 0; i < facets; i++) {
    float a1 = X[G[i][1]] - X[G[i][0]];
    float b1 = Y[G[i][1]] - Y[G[i][0]];
    float c1 = Z[G[i][1]] - Z[G[i][0]];
    float a2 = X[G[i][2]] - X[G[i][0]];
    float b2 = Y[G[i][2]] - Y[G[i][0]];
    float c2 = Z[G[i][2]] - Z[G[i][0]];
    A[i] = b2 * c1 - b1 * c2;
    B[i] = a1 * c2 - c1 * a2;
    C[i] = a2 * b1 - a1 * b2;
    D[i] = -(A[i] * X[G[0][i]] + B[i] * Y[G[0][i]] + C[i] * Z[G[0][i]]);
  }
  centre();
  for (int i = 0; i < facets; i++) {
    plane = A[i] * Xc + B[i] * Yc + C[i] * Zc + D[i];
    if (plane < 0) {
      A[i] = A[i] * (-1);
      B[i] = B[i] * (-1);
      C[i] = C[i] * (-1);
      D[i] = D[i] * (-1);
    }
  }
  for (int i = 0; i < facets; i++) {
    argument = sqrt(A[i] * A[i] + B[i] * B[i] + C[i] * C[i]);
    A[i] /= argument;
    B[i] /= argument;
    C[i] /= argument;
  }
  for (int i = 0; i < facets; i++) {
    V[0][i] = A[i];
    V[1][i] = B[i];
    V[2][i] = C[i];
    V[3][i] = D[i];
  }
  argument_light =
      sqrt(light[0] * light[0] + light[1] * light[1] + light[2] * light[2]);
  for (int i = 0; i < 3; i++) {
    light[i] /= argument_light;
  }
  for (int i = 0; i < 7; i++) {
    shine[i] = 0.25 + 0.75 * ((V[0][i] * light[0] + V[1][i] * light[1] +
                              V[2][i] * light[2] + 1) /
                             2);
  }
}
void draw() {
  roberts();
  lighting();
  BITMAP *buf = create_bitmap(screen_w, screen_h);
  int pol[100];
  clear_to_color(buf, makecol(255, 255, 255));
  int k;
  for (int row = 0; row < facets; row++) {
    k = 0;
    for (int col = 0; col < facets; col++) {
      pol[k++] = int(X[G[row][col]]);
      pol[k++] = int(Y[G[row][col]]);
    }
    if (correct_V[row] >= 0) {
      polygon(buf, R[row], pol,
              makecol(255 * shine[row], 0 * shine[row], 0 * shine[row]));
    }
  }
  blit(buf, screen, 0, 0, 0, 0, screen_w, screen_h);
  destroy_bitmap(buf);
}

void move_up(float W) {
  for (int i = 0; i < tops; i++) {
    Y1 = Y[i] - W;
    Y[i] = Y1;
  }
  draw();
}

void move_down(float S) {
  for (int i = 0; i < tops; i++) {
    Y1 = Y[i] + S;
    Y[i] = Y1;
  }
  draw();
}

void move_left(float A) {
  for (int k = 0; k < tops; k++) {
    X1 = X[k] - A;
    X[k] = X1;
  }
  draw();
}

void move_right(float D) {
  for (int g = 0; g < tops; g++) {
    X1 = X[g] + D;
    X[g] = X1;
  }
  draw();
}
void rotate_x(float K) {
  centre();
  for (int i = 0; i < tops; i++) {
    Y1 = Y[i] * cos(K) - Z[i] * sin(K) + Yc * (1 - cos(K)) + Zc * sin(K);
    Z1 = Y[i] * sin(K) + Z[i] * cos(K) + Zc * (1 - cos(K)) - Yc * sin(K);
    Y[i] = Y1;
    Z[i] = Z1;
  }
  draw();
}

void rotate_y(float M) {
  centre();
  for (int i = 0; i < tops; i++) {
    X1 = X[i] * cos(M) + Z[i] * sin(M) + Xc * (1 - cos(M)) - Zc * sin(M);
    Z1 = -X[i] * sin(M) + Z[i] * cos(M) + Zc * (1 - cos(M)) + Xc * sin(M);
    X[i] = X1;
    Z[i] = Z1;
  }
  draw();
}

void rotate_z(float H) {
  centre();
  for (int i = 0; i < tops; i++) {
    X1 = X[i] * cos(H) - Y[i] * sin(H) + Xc * (1 - cos(H)) + Yc * sin(H);
    Y1 = X[i] * sin(H) + Y[i] * cos(H) + Yc * (1 - cos(H)) - Xc * sin(H);
    X[i] = X1;
    Y[i] = Y1;
  }
  draw();
}
void approximation(float P) {
  centre();
  for (int i = 0; i < tops; i++) {
    X[i] = X[i] * P + Xc * (1 - P);
    Y[i] = Y[i] * P + Yc * (1 - P);
    Z[i] = Z[i] * P + Zc * (1 - P);
  }
  draw();
}
void distance(float N) {
  centre();
  for (int i = 0; i < tops; i++) {
    X[i] = X[i] * 1 / N + Xc * (1 - 1 / N);
    Y[i] = Y[i] * 1 / N + Yc * (1 - 1 / N);
    Z[i] = Z[i] * 1 / N + Zc * (1 - 1 / N);
  }
  draw();
}
void control() {
  while (!key[KEY_ESC]) {

    if (key[KEY_W]) {
      move_up(1.2);
    }
    if (key[KEY_S]) {
      move_down(1.2);
    }
    if (key[KEY_A]) {
      move_left(1.2);
    }
    if (key[KEY_D]) {
      move_right(1.2);
    }
    if (key[KEY_Y]) {
      rotate_x((1 * PI) / 180);
    }
    if (key[KEY_H]) {
      rotate_x(-((1 * PI) / 180));
    }
    if (key[KEY_G]) {
      rotate_y((1 * PI) / 180);
    }
    if (key[KEY_B]) {
      rotate_y(-(1 * PI) / 180);
    }
    if (key[KEY_J]) {
      rotate_z((1 * PI) / 180);
    }
    if (key[KEY_N]) {
      rotate_z(-((1 * PI) / 180));
    }
    if (key[KEY_T]) {
      approximation(1.01);
    }
    if (key[KEY_U]) {
      distance(1.01);
    }
  }
}

END_OF_MAIN();

