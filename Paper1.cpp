#include <stdlib.h>     
#include <iostream>
#include <fstream>

const int N = 10; //Tamaño del arreglo
const int T = 100; // Tiempo máximo
const int Q = 5; //Número de estados

class Material
{
  int h_old[N][N];
  int h_new[N][N];
  int areas[Q];

public:
  void fill (void);
  //void evolution (int &h_old[N][N], int &h_new[N][N]);
  //void size_count (int &areas [Q], &int h_new[N][N];
  void print(const char * Arreglo);
};
  
void Material::fill (void){
  srand(N);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (i == (N-1)){
	h_old[j][i] = h_old[j][0];
      }
      else if (j == (N-1)){
	h_old[j][i] = h_old[0][i];
      }
      else {
	h_old[j][i] = rand() % Q+1;
      }
    }
  }
}
void Material::print(const char * Arreglo){
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      MiArchivo<<h_old[j][i]<<"\t";
    }
    MiArchivo<<std::endl;
  }
  MiArchivo.close();
}

int main (void){
  Material paper1;
  paper1.fill();
  paper1.print("Arreglo.dat");
  
}
