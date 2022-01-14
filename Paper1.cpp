#include <stdlib.h>     
#include <iostream>
#include <fstream>
#include <random>
#include <bits/stdc++.h>


const int N = 200; //Tamaño del arreglo
const int T = 200; // Tiempo máximo
const int r = 40;  // Radio del circulo
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
  void evolution_aux1 (int a, int b, int c, int d, int i, int j);
  void evolution_aux2 (int a, int b, int c, int d, int i, int j, int prob);
  void array_change(void);
  int circle_size_count(void);
  void size_count (void);
  void print_array(const char * Arreglo);
  int getArea (int a);
  void metrics(float &mean_area, float &mean_size, bool print_distr );
};
  
void Material::fill (void){
  std::mt19937 gen(N*N);
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
      if(((i-(N/2))*(i-(N/2)))+((j-(N/2))*(j-(N/2)))<r*r){
	h_old[i][j] = 2;
      }
    }
  }
  
}
void Material::evolution_aux1 (int a, int b, int c, int d, int i, int j){
 
  // a=(i,j-1) ; b = (i-1,j) ; c = (i,j+1) ; d = (i+1,j)
  if (a==b && b==c && c==d){
    h_old[i][j] = a;
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
  }
}

void Material::evolution_aux2 (int a, int b, int c, int d, int i, int j, int prob){
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
	evolution_aux1 (a, b, c,  d,  i, j);
	evolution_aux2 (a, b, c,  d,  i, j, rand(gen));
      }
      else if (j == 0){
	a = h_old[i][N-2];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
	evolution_aux2 (a, b, c,  d,  i, j, rand(gen));
      }
      else {
	a = h_old[i][j-1];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
	evolution_aux2 (a, b, c,  d,  i, j, rand(gen));
      }
    }
  }

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
	evolution_aux1 (a, b, c,  d,  i, j);
      }
      else if (j == 0){
	a = h_old[i][N-2];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
      }
      else {
	a = h_old[i][j-1];
	b = h_old[i-1][j];
	c = h_old[i][j+1];
	d = h_old[i+1][j];
	evolution_aux1 (a, b, c,  d,  i, j);
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
int Material::circle_size_count(void){
  int count = 0;
   for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(h_old [i][j]==2){
	count++;
      }
    }
  }
   return count;
}
  
void Material::size_count (void){
  for(int u =0; u<(Q); u++){
    areas[u] = 0;
  }
  
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      areas[(h_old [i][j])-1]++;
    }
  }
}

void Material::metrics(float &mean_area, float &mean_size, bool print_distr){
  
  float max_freq=0;;
  int size = Q;
  std::vector <int> areas_calc(areas,areas+(size));
  
  sort(areas_calc.begin(),areas_calc.end());
  areas_calc.erase(std::remove(areas_calc.begin(),areas_calc.end(), 0),areas_calc.end());
  areas_calc.shrink_to_fit();
  mean_area = std::accumulate(areas_calc.begin(),areas_calc.end(),0)/areas_calc.size();
  
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
  mean_size = std::accumulate(radius.begin(),radius.end(),0)/radius.size();
  for (int v = 0; v<u; v++){
    radius[v] = log(radius[v]/mean_size);
    }
  if (print_distr == true){
    std::ofstream MiArchivo("area_distr_1.dat");
    for( int v= 0; v<u; v++){
      MiArchivo<<radius[v]<<'\t'<<frecuency[v]<<std::endl;
    }
    MiArchivo.close();
  }
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
  std::ofstream MiArchivo ("size_evolve_1.dat");
  std::ofstream Archivo("Area_circulo.dat");
  Material circulo;
  circulo.fill_circle();
  circulo.print_array("Circulo1.dat");
  float area_circ = 0;
  int area_0 = 0;
  area_0 = circulo.circle_size_count();
  std::cout<<0<<'\t'<<area_0/area_0<<std::endl;
  Archivo<<0<<'\t'<<area_0/area_0<<std::endl;
   for(int t =1 ; t<1600; t++){
    circulo.evolution(t);
    circulo.array_change();
    area_circ = circulo.circle_size_count();
    Archivo<<t<<'\t'<<area_circ/area_0<<std::endl;
    std::cout<<t<<'\t'<<area_circ/area_0<<std::endl;
    if(t == 400){
      circulo.print_array("Circulo2.dat");
    }
    else if (t == 800){
      circulo.print_array("Circulo3.dat");   
    }
    else if (t == 1200){
      circulo.print_array("Circulo4.dat");   
    }
   } 
   circulo.print_array("Circulo5.dat"); 
     
  float mean_area = 0;
  float mean_size = 0;
  
  Material granos;
  granos.fill();
  granos.print_array("C1A1.dat");
  for(int t =0 ; t<=4000; t++){
    granos.evolution(t);
    granos.array_change();
    if(t%100 == 0){
      if (t == 200){
	granos.metrics (mean_area,mean_size, true);
	granos.print_array("C1A2.dat");
      }
      else if( t == 1000){
	granos.print_array("C1A3.dat");
      }
       else if( t == 4000){
	granos.print_array("C1A4.dat");
      }
      granos.size_count();
      granos.metrics (mean_area,mean_size, false);
      std::cout<<t<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
      MiArchivo<<t<<'\t'<<mean_area<<'\t'<<mean_size<<std::endl;
    } 
  }
  MiArchivo.close();
}
