#pragma omp parallel reduction(+:Energy)
 {
   int thid = omp_get_thread_num();
   int nth = omp_get_num_threads();
   int localsize = N/nth;
   int Lmin = thid*localsize;
   int Lmax = Lmin+localsize;
   int L_energy = 0;
     for (int ii= Lmin ;ii<Lmax;ii++){
       for(int jj=Lmin;jj<Lmax;jj++){
	 for(int u=0;u<25;u++){
	   if (h_old[ii][jj]!=neighborhood[u]){
	     L_energy++;
	   } 
	 }
       }
     }
     std::cout<<L_energy<<std::endl;
 }

 for (int ii= 0 ;ii<N;ii++){
    for(int jj=0;jj<N;jj++){
      for(int u=0;u<25;u++){
	if (h_old[ii][jj]!=neighborhood[u]){
	  Energy++;
	} 
      }
    }


 
  long int EB_0 = 0;
  long int EB_f = 0;
  int EB = 0;
  int E_aux = 0;
  int value_aux = h_old[i][j];
  int count = 0;
  int rand_part = 0;
  float Delta_ET = 0.0;
  std::set <int> S;
  std::set <int> Delta_EB;

  std::set<int>::iterator itr;
  
  std::mt19937 gen(time(NULL));
  std::uniform_real_distribution<> dis(0, 1);
  
  Delta_ET = -R*tempt [i][j]*log (dis(gen));
  EB_0 = Boundary_energy ();
  
  for (int u = 0; u < 25; u++) {
    S.insert(neighborhood[u]);
  }
  
  
  for (itr = S.begin(); itr != S.end() ; itr++){
    if (*itr == neighborhood [12]){
      continue;
    }
    else{
      neighborhood [12] = *itr;
      h_old [i][j] = *itr;
    }
    EB_f = Boundary_energy ();
    EB = EB_f-EB_0;
    if (EB<E_aux){
      E_aux=EB;
      h_new[i][j] = h_old[i][j];
    }
  }
  
  if (S.Delta_EB.size() > 1){
    float part = 100.0/Delta_EB.size();
    float rand_num = dis(gen) *100.0;
    rand_part = round(rand_num/part);
    itr = S.begin();
    for (int u = 0; u<rand_part; u++){ 
      itr++;
    }
    h_new[i][j] = *itr;
  }
  
    Delta_E [i][j] = EB-Delta_ET;
    h_old[i][j] = value_aux;
  
}



 
  long int EB_0 = 0;
  long int EB_f = 0;
  int EB = 0;
  int E_aux = 0;
  int value_aux = h_old[i][j];
  int count = 0;
  int rand_part = 0;
  float Delta_ET = 0.0;
  std::vector <int> S(neighborhood,neighborhood+n);
  std::vector <int> Delta_EB;

  std::mt19937 gen(time(NULL));
  std::uniform_real_distribution<> dis(0, 1);
  
  Delta_ET = -R*tempt [i][j]*log (dis(gen));
  EB_0 = Boundary_energy ();
  
  sort(S.begin(),S.end());
  neighborhood [12] = S[0];
  h_old[i][j] = S[0];
  EB_f = Boundary_energy ();
  EB = EB_f-EB_0;
  if (EB<E_aux){
      E_aux=EB;
      h_new[i][j] = h_old[i][j];
    }
  
  for (int u = 1; u<25; u++){
    if(S[u-1]==S[u]){
      
      continue;
    }
    else{
      if (S[u] == value_aux){
	continue;
      }
      else{
	neighborhood [12] = *itr;
	h_old [i][j] = *itr;
      }
      EB_f = Boundary_energy ();
      EB = EB_f-EB_0;
      if (EB<E_aux){
	E_aux=EB;
	h_new[i][j] = h_old[i][j];
      }
      Delta_EB.insert(EB);
    }
  }
  
  if (Delta_EB.size() > 1){
    float part = 100.0/Delta_EB.size();
    float rand_num = dis(gen) *100.0;
    rand_part = round(rand_num/part);
    itr = S.begin();
    for (int u = 0; u<rand_part; u++){ 
      itr++;
    }
    h_new[i][j] = *itr;
  }
  
    Delta_E [i][j] = EB-Delta_ET;
    h_old[i][j] = value_aux;
  
}
