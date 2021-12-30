#include <stdlib.h>     
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <time.h>

const short int N = 200; //Tamaño del arreglo
const short int T = 200; // Tiempo máximo
const short int Q = 100; //Número de estados
const short int r = 5; //Radio en celdas de la región caliente
const float E_A = 2; // Energía de activación, como multiplo de la constante de Boltzmann
const float R = 8.31446261; // Constante del Gas ideal

class Material
{
private:
  short int h_old[N][N];
  short int h_new[N][N];
  int Delta_E_B [25];
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
  bool is_boundary (int i, int j);
  long int Boundary_energy (void);
  void Delta_E_min (int i, int j);
  void evolve (void);
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
	h_old[i][j]= rand(gen1);
	h_new[i][j]=rand(gen1);
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
void Material::neighbor_def (int i, int j){
  int n_index=0; //indice que reccore los vecinos
  for (n_index = 0; n_index<25; n_index++){
    neighborhood[n_index] = 0;
  }
  // Frontera interior
  if((i == 1 || i == (N-2)) && (j>=1 && j<=(N-2))){
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
  else if ((j==1||j==(N-2)) && (i > 1 && i < (N-2))){
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
    if (i==0){//terminado
      if(j==0){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j; jj <= j+2; jj++){
	    neighborhood[12+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u < 3; u++){
	  neighborhood [7+u] = h_old [N-2][u];
	  neighborhood [2+u] = h_old [N-3][u];
	  neighborhood [(5*u)+11] = h_old [u][N-2];
	  neighborhood [(5*u)+10] = h_old [u][N-3];
	}
	neighborhood [0] = h_old [N-3][N-3];
	neighborhood [1] = h_old [N-3][N-2];
	neighborhood [5] = h_old [N-2][N-3];
	neighborhood [6] = h_old [N-2][N-2];
      }
      else if (j==1){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-1; jj <= j+2; jj++){
	    neighborhood[11+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for(int u = 1; u < 5; u++){
	  neighborhood [u] = h_old [N-3][u-1];
	  neighborhood [5+u] = h_old [N-2][u-1];
	  neighborhood [(5*u)+5] = h_old [u-1][N-2];
	}
	neighborhood [0] = h_old [N-3][N-2];
	neighborhood [5] = h_old [N-2][N-2];
      }
      
      else if (j==(N-1)){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-2; jj <= j; jj++){
	    neighborhood[10+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [5+u] = h_old [N-2][(N-3)+u];
	  neighborhood [u] = h_old [N-3][(N-3)+u];
	  neighborhood [(5*u)+13] = h_old [u][1];
	  neighborhood [(5*u)+14] = h_old [u][2];
	}	
	neighborhood [3] = h_old [N-3][1];
	neighborhood [4] = h_old [N-3][2];
	neighborhood [8] = h_old [N-2][1];
	neighborhood [9] = h_old [N-2][2];
      }
      else if (j == (N-2)){ //terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-2; jj <= j+1; jj++){
	    neighborhood[10+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for (int u = 0 ; u<4 ; u++){
	  neighborhood [u] = h_old [N-3][(N-4)+u];
	  neighborhood [5+u] = h_old [N-2][(N-4)+u];
	  neighborhood [(5*u)+14] = h_old [u][1];
	}
	neighborhood [4] = h_old [N-3][1];
	neighborhood [9] = h_old [N-2][1];
      }
      
      else {// terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i; ii <= i+2; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[10+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
        i_aux = 0;
	j_aux = 0;
	for (int ii = 0; ii < 2; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[N-3+ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
      }
    }
    
    else if (i == (N-1)){
      if (j==0){ // terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j; jj <= j+2; jj++){
	    neighborhood[2+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [5*u] = h_old [(N-3)+u][N-3];
	  neighborhood [(5*u)+1] = h_old [(N-3)+u][N-2];
	  neighborhood [17+u] = h_old [1][u];
	  neighborhood [22+u] = h_old [2][u];
	}
	neighborhood [15] = h_old [1][N-3];
	neighborhood [16] = h_old [1][N-2];
	neighborhood [20] = h_old [2][N-3];
	neighborhood [21] = h_old [2][N-2];
      }
      else if (j==1){//terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-1; jj <= j+2; jj++){
	    neighborhood[1+(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for (int u =0; u < 4; u++){
	  neighborhood [16+u] = h_old [1][u];
	  neighborhood [21+u] = h_old [2][u];
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [5*u] = h_old [(N-3)+u][N-2];
	}
	neighborhood [15] = h_old [1][N-2];
	neighborhood [20] = h_old [2][N-2];
      }
      else if (j==N-1){//terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-2; jj <= j; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux--;
	  i_aux++;
	}
	for (int u = 0; u < 3; u++){
	  neighborhood [15+u] = h_old [1][(N-3)+u];
	  neighborhood [20+u] = h_old [2][(N-3)+u];
	  neighborhood [(5*u)+3] = h_old [(N-3)+u][1];
	  neighborhood [(5*u)+4] = h_old [(N-3)+u][2];
	}
	neighborhood [18] = h_old [1][1];
	neighborhood [19] = h_old [1][2];
	neighborhood [23] = h_old [2][1];
	neighborhood [24] = h_old [2][2];
      }
      else if (j == N-2){//terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-2; jj <= j+1; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=2;
	  i_aux++;
	}
	for (int u =0; u < 4; u++){
	  neighborhood [15+u] = h_old [1][(N-4)+u];
	  neighborhood [20+u] = h_old [2][(N-4)+u];
	}
	for (int u = 0; u<3 ; u++){
	  neighborhood [(5*u)+4] = h_old [(N-3)+u][1];
	}
	neighborhood [19] = h_old [1][1];
	neighborhood [24] = h_old [2][1];
      }
      else { // terminado
	int i_aux = 0;
	int j_aux = 0;
	for (int ii = i-2; ii <= i; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[(3*i_aux)+j_aux] = h_old[ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
        i_aux = 0;
	j_aux = 0;
	for (int ii = 0; ii < 2; ii++){
	  for(int jj = j-2; jj <= j+2; jj++){
	    neighborhood[15+(3*i_aux)+j_aux] = h_old[1+ii][jj];
	    j_aux++;
	  }
	  j_aux-=3;
	  i_aux++;
	}
      }
    }
  }
   else  if ((i == 1 || i == N-2) && (j == 0 || j == N-1)){
     if ( i == 1 && j==0){ //terminado
      int i_aux = 0;
      int j_aux = 0;
      for (int ii = i-1; ii <= i+2; ii++){
	for(int jj = j; jj <= j+2; jj++){
	  neighborhood[7+(4*i_aux)+j_aux] = h_old[ii][jj];
	  j_aux++;
	}
	j_aux-=2;
	i_aux++;
      }
      for (int u = 0; u<4 ; u++){
	neighborhood [(5*u)+5] = h_old [u][N-3];
	neighborhood [(5*u)+6] = h_old [u][N-2];
      }
      for(int u = 0; u<3 ; u++){
	neighborhood [2+u] = h_old [N-2][u];
      }
      neighborhood [0] = h_old [N-2][N-3];
      neighborhood [1] = h_old [N-2][N-2];
     }
     else if (i==1 && j==N-1){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-1; ii <= i+2; ii++){
	 for(int jj = j-2; jj <= j; jj++){
	   neighborhood[5+(4*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=2;
	i_aux++;
       }
       for (int u = 0; u<4 ; u++){
	 neighborhood [(5*u)+8] = h_old [u][1];
	 neighborhood [(5*u)+9] = h_old [u][2];
       }
       for(int u = 0; u<3 ; u++){
	 neighborhood [u] = h_old [N-2][(N-3)+u];
       }
       neighborhood [3] = h_old [N-2][1];
       neighborhood [4] = h_old [N-2][2];
     }
     else if (i==(N-2) && j == 0){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-2; ii <= i+1; ii++){
	 for(int jj = j; jj <= j+2; jj++){
	   neighborhood[2+(4*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
       for (int u = 0; u<4 ; u++){
	 neighborhood [(5*u)] = h_old [(N-4)+u][N-3];
	 neighborhood [(5*u)+1] = h_old [(N-4)+u][N-2];
       }
        for(int u = 0; u<3 ; u++){
	 neighborhood [22+u] = h_old [1][u];
       }
       neighborhood [20] = h_old [1][N-3];
       neighborhood [21] = h_old [1][N-2];
     }
     else if (i == (N-2) && j== (N-1)){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-2; ii <= i+1; ii++){
	 for(int jj = j-2; jj <= j; jj++){
	   neighborhood[(4*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
       for (int u = 0; u<4 ; u++){
	 neighborhood [(5*u)+3] = h_old [(N-4)+u][1];
	 neighborhood [(5*u)+4] = h_old [(N-4)+u][2];
       }
       for(int u = 0; u<3 ; u++){
	 neighborhood [20+u] = h_old [1][(N-3)+u];
       }
       neighborhood [23] = h_old [1][1];
       neighborhood [24] = h_old [1][2];
     }
   }
   else if ((j==0 || j==(N-1)) && (i>1 && i <N-2)){
     if (j == 0) {//terminado
       int i_aux = 0;
       int j_aux = 0;
        for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = j; jj <= j+2; jj++){
	   neighborhood[2+(5*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=3;
	 i_aux++;
       }
       i_aux = 0;
       j_aux = 0;
       for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = 0; jj < 2; jj++){
	   neighborhood[(5*i_aux)+j_aux] = h_old[ii][N-3+jj];
	    j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
     }
     else if (j == (N-1)){//terminado
       int i_aux = 0;
       int j_aux = 0;
       for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = j-2; jj <= j; jj++){
	   neighborhood[(5*i_aux)+j_aux] = h_old[ii][jj];
	   j_aux++;
	 }
	 j_aux-=3;
	 i_aux++;
       }
       i_aux = 0;
       j_aux = 0;
       for (int ii = i-2; ii <= i+2; ii++){
	 for(int jj = 0; jj < 2; jj++){
	   neighborhood[3+(5*i_aux)+j_aux] = h_old[ii][1+jj];
	    j_aux++;
	 }
	 j_aux-=2;
	 i_aux++;
       }
     }
   }
   
   else{
     n_index=0;
     for(int ii = i-2; ii <= i+2; ii++){
       for(int jj = j-2; jj <= j+2; jj++){
	 neighborhood[n_index] = h_old[ii][jj];
	 n_index++;
       }
     }
   }
}

bool Material::is_boundary (int i, int j){
  int actual = h_old [i][j];
  int test = 0;
  int i_aux =  0;
  int j_aux = 0;
  neighbor_def (i,j);
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj<3; jj++){
      test = neighborhood [6+(3*i_aux)+j_aux];
      if (test != actual){
	return true;
      }
      j_aux++;
    }
    j_aux--;
    i_aux++; 
  }
  return false;
}

long int Material::Boundary_energy (void){
 long int Energy = 0;
  for (int ii=0;ii<N;ii++){
    for(int jj=0;jj<N;jj++){
      for(int u=0;u<25;u++){
	if (h_old[ii][jj]!=neighborhood[u]){
	  Energy++;
	  //  std::cout<<Energy<<std::endl;
	} 
      }
    }
  }
  return Energy;
}

void Material::Delta_E_min (int i, int j){
  
  long int EB_0 = 0;
  long int EB_f = 0;
  int EB = 0;
  int E_aux = 0;
  int value_aux = h_old[i][j];
  int count = 0;
  int rand_part = 0;
  float Delta_ET = 0.0;
  

  std::mt19937 gen(time(NULL));
  std::uniform_real_distribution<> dis(0, 1);
  
  Delta_ET = -R*tempt [i][j]*log (dis(gen));
  EB_0 = Boundary_energy ();
  
  for (int u = 0; u < 25 ; u++){
    if(u==12){
      neighborhood [12] = value_aux;
      h_old [i][j] = value_aux;
    }
    else{
      neighborhood [12] = neighborhood [u];
      h_old [i][j] = neighborhood [u];
    }
    EB_f = Boundary_energy ();
    EB = EB_f-EB_0;
    if (EB<E_aux){
      E_aux=EB;
      h_new[i][j] = h_old[i][j];
    }
    Delta_E_B [u] = EB;
  }
  
  for(int u = 0; u<25; u++){
    if(Delta_E_B[u] == E_aux){
      count++;
      }
  }
  // 
  if (count > 1){
    float part = 100.0/count;
    float rand_num = dis(gen) *100.0;
    rand_part = round(rand_num/part);
    
    if (rand_part==0){
      rand_part++;
    }
    count = 0;
    for(int u = 0; u<25; u++){
      if(Delta_E_B[u] == E_aux){
	count++;
	if(count == rand_part){
	  
	  h_new [i][j] = neighborhood [u];
	  break;
	}
      }  
    }
  }
  Delta_E [i][j] = EB-Delta_ET;
  h_old[i][j] = value_aux;
}

void Material::evolve (void){
  float D_Emin = 0.0;
  float M_max = exp ((-E_A)/tempt [N/2][N/2]);

  std::mt19937 gen(time(NULL));
  
  
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (is_boundary(i,j)==true){
	Delta_E_min(i,j);
      }
      else {
	Delta_E [i][j] = 0;
      }
    }
  }
  D_Emin = Delta_E [0][0];
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){  
      if(Delta_E [i][j] < D_Emin){
	D_Emin = Delta_E[i][j];
      }
    }
  }
  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (Delta_E [i][j] != 0){
	std::bernoulli_distribution prob(/*probability (i, j, M_max, D_Emin)*/1);
	//std::cout<<probability (i, j, M_max, D_Emin)<<std::endl;
	if (prob(gen) == false){
	  //std::cout<<"aqui"<<std::endl;
	  h_new[i][j] = h_old [i][j];
	}
      }
    }
  }

  for (int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if (j == N-1){
	h_old [i][j] = h_new[i][0];
      }
      else if (i== N-1){
	h_old [i][j] = h_new [0][j];
      }
      else {
	h_old [i][j] = h_new [i][j];
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
  int a = 0;
  Material granos;
  granos.fill_tempt(1000.0,1000.0);
  granos.fill();
  granos.print_array("Arreglo1.dat");
  for (int t = 0; t <2; t++){
    granos.evolve();
  }
  granos.print_array("Arreglo2.dat");
  
  return 0;
  
}
