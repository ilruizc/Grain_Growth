#include <stdlib.h>     
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>

const short int N = 20; //Tamaño del arreglo
const short int T = 200; // Tiempo máximo
const short int Q = 10; //Número de estados
const short int r = 6; //Radio en celdas de la región caliente
const float E_A = 2; // Energía de activación, como multiplo de la constante de Boltzmann 

class Material
{
private:
  short int h_old[N][N];
  short int h_new[N][N];
  float Delta_E[N][N];
  short int neighborhood [24]; //Numero de segundos vecinos de Moore
  float tempt [N][N];
  int areas[Q];

public:
  void fill (void);
  void fill_tempt (float T_min, float T_max);
  float probability (int i, int j, float M_max, float D_Emin);
  // void size_count (void);
  void print_array(const char * Arreglo);
  void print_gradient(const char * Arreglo);
  // int getArea (int a);
};

void Material::fill (void){
  std::mt19937 gen(N);
  std::mt19937 gen1(N+1);
  std::uniform_int_distribution<> rand(1, Q);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      h_old[i][j] = 1;
      Delta_E [i][j] = rand(gen)-rand(gen1);
    }
  }
  
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (i == (N-1)){
	h_old[i][j] = h_old[0][j];
      }
      else if (j == (N-1)){
	h_old[i][j] = h_old[i][0];
      }
      else {
	  h_old[i][j]= rand(gen);
	
      }
    }
  }
  }
void Material::fill_tempt(float T_min, float T_max){
  float a = r*(T_max-T_min);
  float b = T_min;
  float x = 0.0;
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      tempt [i][j] = 500;
    }
  }
  
  for(int i=1; i<(N-1); i++){
    for(int j=1; j<(N-1); j++){
      if((((i-(N/2))*(i-(N/2)))+((j-(N/2))*(j-(N/2))))<(r*r)){
	tempt [i][j] = T_max ;
      }
       else{
	 x = sqrt(((1.0*i)-(1.0*N)/2)*((1.0*i)-(1.0*N)/2)+((1.0*j)-(1.0*N)/2)*((1.0*j)-(1.0*N)/2));
	 tempt [i][j] = a/(x)+b;
       }
    }
  }
}

float Material::probability (int i, int j, float M_max, float D_Emin){
  float M_ij = exp ((-E_A)/tempt [i][j]);
  if (Delta_E[i][j] < 0){
  return ((M_ij/M_max)*(Delta_E[i][j]/D_Emin));
  }
  else {
    return 0;
  }
}
void Material::print_array(const char * Arreglo){
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      MiArchivo<<tempt [i][j]<<"  ";
    }
    MiArchivo<<std::endl;
  }
  MiArchivo.close();
}
void Material::print_gradient(const char * Arreglo)
{
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<=N/2;i++){
    MiArchivo<<i<<'\t'<<tempt [i][N/2]<<std::endl;;
  }
  MiArchivo.close();
}
int main (void){
   Material granos;
   granos.fill_tempt(800.0,1000.0);
   granos.fill();
   for(int i = 0; i < N; i++){
     for( int j = 0; j < N; j++){
       std::cout<<granos.probability(i,j, 0.8, -2*Q)<<'\t';
     }
     std::cout<<std::endl;
   }
   
  return 0;
  
}
