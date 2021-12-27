#include <stdlib.h>     
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>

const short int N = 20; //Tamaño del arreglo
const short int T = 200; // Tiempo máximo
const short int Q = 5; //Número de estados
const short int r = 6; //Radio en celdas de la región caliente
const float E_A = 2; // Energía de activación, como multiplo de la constante de Boltzmann 

class Material
{
private:
  short int h_old[N][N];
  short int h_new[N][N];
  float Delta_E[N][N];
  short int neighborhood [25]; //Numero de segundos vecinos de Moore
  float tempt [N][N];
  int areas[Q];

public:
  void fill (void);
  void fill_tempt (float T_min, float T_max);
  float probability (int i, int j, float M_max, float D_Emin);
  void neighbor_def (int i, int j);
  int neighbor_get (int ii);
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
  h_old[N/2][N/2] = 0;
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
void Material::neighbor_def (int i, int j){
  int n_index=0; //indice que reccore los vecinos
  for (n_index = 0; n_index<25; n_index++){
    neighborhood[n_index] = 0;
  }
  // Frontera interior
  if(i == 1 || i == (N-2)){
    int i_aux = 0;
    int j_aux = 0;
    for (int ii = i-1; ii <= i+1; ii++){
      for(int jj = j-1; jj <= j+1; jj++){
	neighborhood[6+(3*i_aux)+j_aux] = h_old[ii][jj];
	j_aux++;
      }
      j_aux--;
      i_aux++;
    }
    if (i == 1){
      if(j==1){
	neighborhood[0] = h_old[N-2][N-2];
	for (int u = 1; u <= 4; u++){
	  neighborhood [u*5] = h_old[u-1][N-2];
	  neighborhood [(u*5)+4] = h_old [u-1][3];
	  neighborhood [20+u] = h_old[3][u-1];
	  neighborhood [u] = h_old [N-2] [u-1];
	}
      }
      else if (j==(N-2)){
	neighborhood[4] = h_old[1][1];
	for (int u = 1; u <= 4; u++){
	  neighborhood[u*5] = h_old[u-1][N-4];
	  neighborhood [(u*5)+4] = h_old [u-1][1];
	  neighborhood [19+u] = h_old[3][(N-5)+u];
	  neighborhood [u-1] = h_old[N-2][(N-5)+u];
	}
      }
      else{
	for (int u = 0; u < 5; u++){
	  neighborhood [u] = h_old[N-2][(j-2)+u];
	  neighborhood [5+u] = h_old[i-1][(j-2)+u];
	  neighborhood [10+u] = h_old [i][(j-2)+u];
	  neighborhood [15+u] = h_old [i+1][(j-2)+u];
	  neighborhood [20+u] = h_old [i+2][(j-2)+u];
	}
      }
    }
    else if(i == (N-2)){
      if (j==1){
	neighborhood[20] = h_old [1][N-2];
	for (int u = 1; u <= 4; u++){
	  neighborhood [((u-1)*5)] = h_old [(N-5)+u][N-2];
	  neighborhood [((u-1)*5)+4] = h_old [(N-5)+u][3];
	  neighborhood [u] = h_old [N-4][u-1];
	  neighborhood [20+u] = h_old [1][u-1];
	}
      }
      else if (j==(N-2)){
	neighborhood[24]=h_old[1][1];
	for(int u = 1; u<=4; u++){
	  neighborhood [u-1] = h_old[N-4][(N-5)+u];
	  neighborhood [(u-1)*5] = h_old [(N-5)+u][N-4];
	  neighborhood [((u-1)*5)+4] = h_old [(N-5)+u][1];
	  neighborhood [19+u] = h_old [1][(N-5)+u];
	}
      }
    } 
  }
  else if ((j==1||j==(N-2)) && (i != 1 || j != (N-2))){
    if (j==1){
      for (int u = 0; u<5 ; u++){
	neighborhood [(u*5)] = h_old[(i-2)+u][N-2];
	neighborhood [(u*5)+1] = h_old [(i-2)+u][j-1];
	neighborhood [(u*5)+2] = h_old [(i-2)+u][j];
	neighborhood [(u*5)+3] = h_old [(i-2)+u][j+1];
	neighborhood [(u*5)+4] = h_old [(i-2)+u][j+2];
      }
    }
    else if (j==(N-2)){
      for (int u = 0; u<5; u++){
	neighborhood [(u*5)] = h_old [(i-2)+u][j-2];
	neighborhood [(u*5)+1] = h_old [(i-2)+u][j-1];
	neighborhood [(u*5)+2] = h_old [(i-2)+u][j];
	neighborhood [(u*5)+3] = h_old [(i-2)+u][j+1];
	neighborhood [(u*5)+4] = h_old[(i-2)+u][1];
      }
    }
  }
  //Frontera exterior
  else if (i == 0 || i == (N-1)){ 
    if (i==0){
      int i_aux = 0;
      int j_aux = 0;
      if(j==0){
	for (int ii = ; ii <= i+1; ii++){
	  for(int jj = j-1; jj <= j+1; jj++){
	    neighborhood[6+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
      }
    }
  }
  else{
    for(int ii = i-2; ii <= i+2; ii++){
      for(int jj = j-2; jj <= j+2; jj++){
	neighborhood[n_index] = h_old[ii][jj];
	n_index++;
      }
    }
  }
  
}

int Material::neighbor_get (int ii){
  return neighborhood[ii];
}

void Material::print_array(const char * Arreglo){
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      MiArchivo<<h_old [i][j]<<"  ";
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
   granos.print_array("Arreglo1.txt");
   granos.neighbor_def(N/2,N-2);
   int neighbor [25];
   for (int ii = 0; ii<25 ; ii++){
     neighbor[ii] = granos.neighbor_get (ii);
   }
   for(int i=0;i<5;i++){
    for(int j=0;j<5;j++){
      std::cout<<neighbor[(5*i)+j]<<'\t';
    }
    std::cout<<std::endl;
   }
   
  return 0;
  
}
