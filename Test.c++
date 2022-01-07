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
  
