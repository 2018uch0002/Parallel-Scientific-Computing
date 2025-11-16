#include <stdio.h>

void func(int** u){
  for (int i = 0; i < 10; i++){
    printf("u[0] = %d \t u[1] = %d \n", u[0][i], u[1][i]);
  }
}

int main()
{
  int *u[2];
  u[0] = new int [10];
  u[1] = new int [10];

  for (int i = 0; i < 10; i++){
    u[0][i] = i;
    u[1][i] = 100+i;
  }

  func(u);
}
