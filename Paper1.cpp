#include <stdlib.h>     
#include <iostream>
#include <fstream>
#include <random>
#include <bits/stdc++.h>
const int N = 200; //Tamaño del arreglo
const int T = 200; // Tiempo máximo
const int Q = 1000; //Número de estados

class Material
{
private:
  int h_old[N][N];
  int h_new[N][N];
  int areas[Q];

public:
  void fill (void);
  void fill_circle(void);
  void evolution(int t);
  void evolution_aux (int a, int b, int c, int d, int i, int j, int prob);
  void array_change(void);
  void size_count (void);
  void print_array(const char * Arreglo);
  int getArea (int a);
  void metrics(float &mean_area, float &mean_size, const char * distribucion);
};
  
void Material::fill (void){
  std::mt19937 gen(0);
  std::binomial_distribution <> d(1, 1);
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
 
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      h_old[i][j] = 1;
    }
  }
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(((i-(N/2))*(i-(N/2)))+((j-(N/2))*(j-(N/2)))<(N/6.7)*(N/6.7)){
	h_old[i][j] = 2;
      }
    }
  }
  
}
void Material::evolution_aux (int a, int b, int c, int d, int i, int j, int prob){
 
  // a=(i,j-1) ; b = (i-1,j) ; c = (i,j+1) ; d = (i+1,j)
  if (a==b && b==c && c==d){
    h_new[i][j]=h_old[i][j];
  }

  else{
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
      
    }
    
    else if((c == d) && (d == a)){//Prueba si c,d,a son iguales
      if(i==0){
	h_new[i][j] = c;
	h_new[N-1][j] = c;
      }
      else if(j==0){
	h_new[i][j] = c;
	h_new[i][N-1] = c;
      }
      else{
	h_new[i][j] = c;
      }
      
    }
    else if((d == a) && (a == b)){//Prueba si d,a,b son iguales
      if(i==0){
	h_new[i][j] = d;
	h_new[N-1][j] = d;
      }
      else if(j==0){
	h_new[i][j] = d;
	h_new[i][N-1] = d;
      }
      else{
	h_new[i][j] = d;
      }
      
    }
    else {
      //std::cout<<prob<<std::endl;
      if(1<=prob && prob<=25){
	h_new[i][j] = a;
      }
      else if (25<prob && prob<=50){
	h_new[i][j] = b;
      }
      else if (50<prob && prob<=75){
	h_new[i][j] = c;
      }
      else if (75<prob && prob<=100) {
	h_new[i][j] = d;
      }
      else{
	h_new[i][j] = h_old[i][j];
      }
    }
  }
}
void Material::evolution(int t){
  int a,b,c,d; //variables auxiliares: a=(i,j-1) ; b = (i-1,j) ; c = (i,j+1) ; d = (i+1,j)
 std::mt19937 gen(t);
  std::uniform_int_distribution<> rand(1, 100);
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
	evolution_aux (a, b, c,  d,  i, j, rand(gen));
      }
      else if (j == 0){
	a = h_old[i][N-2];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux (a, b, c,  d,  i, j, rand(gen));
      }
      else {
	a = h_old[i][j-1];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux (a, b, c,  d,  i, j, rand(gen));
      }
    }
  } 
}
void Material::array_change(void){
 for(int i=0; i<(N-1); i++){
    for(int j=0; j<(N-1); j++){
      if (i==0){
	if(j==0){
	  h_old[i][j] = h_new[i][j];
	h_old[N-1][j] = h_new[i][j];
	h_old[i][N-1] = h_new[i][j];
	}
	else{
	  h_old[i][j] = h_new[i][j];
	  h_old[N-1][j] = h_new[i][j];
	}
      }
      else if(j==0 && i!=0){
	h_old[i][j] = h_new[i][j];
	h_old[i][N-1] = h_new[i][j];
      }
      else{
	h_old[i][j] = h_new[i][j];
      }
    }
 }
 h_old[N-1][N-1] = h_old[0][0];
}
void Material::size_count (void){
  for(int u =0; u<(Q); u++){
    areas[u] = 0;
  }
  
  for (int i=0; i<N-1; i++){
    for(int j=0; j<N-1; j++){
      areas[(h_old [i][j])-1]++;
    }
  }
}

void Material::metrics(float &mean_area, float &mean_size, const char * distribucion){
  std::ofstream MiArchivo(distribucion);
  int max_freq=0;;
  int size = Q;
  std::vector <int> areas_calc(areas,areas+(size));
  
  sort(areas_calc.begin(),areas_calc.end());
  areas_calc.erase(std::remove(areas_calc.begin(),areas_calc.end(), 0),areas_calc.end());
  areas_calc.shrink_to_fit();
  mean_area = std::reduce(areas_calc.begin(),areas_calc.end())/areas_calc.size();
  
  std::map<int, float> counts;
  for (auto v : areas_calc)
    ++counts[v];
  
  
  for (auto v : areas_calc){
    if(max_freq<counts[v]){
      max_freq = counts[v];
    }
  }
  std::vector<float> radius(counts.size());
  std::vector<float> frecuency(counts.size());
  int u = 0;
  for (auto const &p : counts ){
    radius [u] = sqrt(p.first);
    frecuency[u] = p.second/max_freq;
    u++;
  }
  mean_size = std::reduce(radius.begin(),radius.end())/radius.size();
  for (int v = 0; v<u; v++){
    radius[v] = log(radius[v]/mean_size);
    }
  for( int v= 0; v<u; v++){
    MiArchivo<<radius[v]<<'\t'<<frecuency[v]<<std::endl;
  }
  MiArchivo.close();
}


void Material::print_array(const char * Arreglo){
  std::ofstream MiArchivo(Arreglo);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      MiArchivo<<h_old[i][j]<<"  ";
    }
    MiArchivo<<std::endl;
  }
  MiArchivo.close();
}
int Material::getArea (int a){
  return areas[a];
}

int main (void){
  std::ofstream MiArchivo ("time_evolution.dat");
  /*std::ofstream Archivo("Area_circulo.dat");
  Material circulo;
  circulo.fill_circle();
  circulo.print_array("Circulo1.dat");
  int Area_max = 0;
  circulo.size_count();
  Area_max = circulo.getArea(1);
  // std::cout<<Area_max<<std::endl;
   for(int t =0 ; t<1600; t++){
    circulo.evolution(t);
    circulo.array_change();
    circulo.size_count();
    Archivo<<t<<'\t'<<(1.0*circulo.getArea(1))/(1.0*Area_max)<<std::endl;
    }
  circulo.print_array("Circulo2.dat");
  Archivo.close();*/
  float mean_area = 0;
  float mean_size = 0;
  
  Material granos;
  granos.fill();
  granos.print_array("Arreglo1.dat");
  // granos.area_distribution("area_distr.dat");
  for(int t =0 ; t<10000; t++){
    granos.evolution(t);
    granos.array_change();
    granos.size_count();
    if(t%100 == 0){
    granos.metrics (mean_area,mean_size, "area_distr.dat");
    std::cout<<t<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
    MiArchivo<<t<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
    }
  }
  granos.print_array("Arreglo2.dat");
  MiArchivo.close();
}
