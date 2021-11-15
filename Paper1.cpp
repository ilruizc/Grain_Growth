#include <stdlib.h>     
#include <iostream>
#include <fstream>
#include <random>

const int N = 20; //Tamaño del arreglo
const int T = 200; // Tiempo máximo
const int Q = 10; //Número de estados

class Material
{
  int h_old[N][N];
  int h_new[N][N];
  int areas[Q];

public:
  void fill (void);
  void fill_circle(void);
  void evolution(void);
  void evolution_change (int a, int b, int c, int d, int i, int j);
  //void size_count (int &areas [Q], &int h_new[N][N];
  void print(const char * Arreglo);
};
  
void Material::fill (void){
  srand(N);
  std::mt19937 gen(N);
  std::binomial_distribution <> d(1, 0.4);
  std::uniform_int_distribution<> rand(1, Q);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      h_old[i][j] = 1;
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
	if(d(gen) == 1){
	  h_old[i][j]= rand(gen);
	}
      }
    }
  }
}
void Material::fill_circle(void){
  srand(N);
  for(int i=N; i<N; i++){
    for(int j=N; j<N; j++){
      h_old[i][j] = 1;
    }
  }
  for(int i=((N/2)-20); i<((N/2)+20); i++){
    for(int j=((N/2)-20); j<((N/2)+20); j++){
      h_old[i][j] = rand()%2;
    }
  }

}
void Material::evolution_change (int a, int b, int c, int d, int i, int j){
   
  if((a == b) && (b == c)){//Prueba si a,b,c son iguales
    if(i==0){
      h_new[i][j] = a;
      h_new[N-1][j] = a;
    }
    else if(j==0){
      h_new[i][j] = a;
      h_new[i][j] = a;
    }
    else{
      h_new[i][j] = a;
    }
    // std::cout<<a<<std::endl;
  }
  
  else  if((b == c) && (c == d)){//Prueba si b,c,d son iguales
    if(i==0){
      h_new[i][j] = b;
      h_new[N-1][j] = b;
    }
    else if(j==0){
      h_new[i][j] = b;
      h_new[i][j] = b;
    }
    else{
      h_new[i][j] = b;
    }
    // std::cout<<b<<std::endl;
  }
  
  else if((c == d) && (d == a)){//Prueba si c,d,a son iguales
    if(i==0){
      h_new[i][j] = c;
      h_new[N-1][j] = c;
    }
    else if(j==0){
      h_new[i][j] = c;
      h_new[i][j] = c;
    }
    else{
      h_new[i][j] = c;
    }
    // std::cout<<c<<std::endl;
  }
  else if((d == a) && (a == b)){//Prueba si d,a,b son iguales
    if(i==0){
      h_new[i][j] = d;
      h_new[N-1][j] = d;
    }
    else if(j==0){
      h_new[i][j] = d;
      h_new[i][j] = d;
    }
    else{
      h_new[i][j] = d;
    }
    // std::cout<<d<<std::endl;
  }
  else {
    //  std::cout<<h_old[i][j]<<"\t";
    h_new[i][j] = h_old[i][j];
    // std::cout<<h_new[i][j]<<std::endl;
  }
}
void Material::evolution(void){
  int a,b,c,d; //variables auxiliares: a=(i,j-1) ; b = (i-1,j) ; c = (i,j+1) ; d = (i+1,j)

  for(int i=0; i<(N-1); i++){
    for(int j=0; j<(N-1); j++){
      if (i == 0){
	if (j== 0){
	  a = h_old[i][N-2];
	}
	else{
	  a = h_old[i][j-1];
	}
	b = h_old[N-2][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_change (a, b, c,  d,  i, j);
      }
      else if (j == 0){
	a = h_old[i][N-2];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_change (a, b, c,  d,  i, j);
      }
      else {
	a = h_old[i][j-1];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_change (a, b, c,  d,  i, j);
      }
      //  std::cout<<a<<" "<<b<<" "<<c<<" "<<d<<std::endl;
    }
  }
  for(int i=0; i<(N-1); i++){
    for(int j=0; j<(N-1); j++){
      h_old[i][j] = h_new[i][j];
    }
  }
}
void Material::print(const char * Arreglo){
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      MiArchivo<<h_old[i][j]<<"  ";
    }
    MiArchivo<<std::endl;
  }
  MiArchivo.close();
}

int main (void){
  Material paper1;
  paper1.fill();
  paper1.print("Arreglo1.dat");
  for(int t =0 ; t<T; t++){
    paper1.evolution();
  }
  paper1.print("Arreglo2.dat");
  
}
