#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <iostream> /* time */

int main ()
{
  int iSecret;
  int count [10]={0};
  /* initialize random seed: */
  srand (time(NULL));


  for (int i =0; i<100000 ; i++){
    iSecret = rand() % 10 ;
    count [iSecret]++;
    // std::cout<< iSecret<<'\n';
  }
  for(int i = 0; i<10; i++){
     std::cout<<i<<'\t'<<count[i]<<'\n';
  }
  return 0;
}
