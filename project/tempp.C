#include <iostream>

void fun(double* &vec){
  vec = new double[10];

  for (std::size_t i = 0; i < 10; i++){
    vec[i] = 10.*i;
  }
}

int main(){
  double *vec;
  fun(vec);
  for (std::size_t i = 0; i < 10; i++){
    std::cout << vec[i] << "\t";
  }
  std::cout << std::endl;
}