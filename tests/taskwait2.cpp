// CORRECT example of using tasks with dependency.


#include<iostream>

using namespace std;

int main(int argc, char**argv){
  int dep[100];

#pragma omp parallel
  {
#pragma omp single nowait
    {
      for(int i = 0; i < 100; i++){
        if( i == 0 ){
#pragma omp task  depend(out: dep[i])
          {
            cout << "i = " <<  i << endl;
          }
        }
        else{
#pragma omp task depend(in: dep[i-1]) depend(out: dep[i])
          {
            cout << "i = " <<  i << endl;
          }
        }
      } // for
    }
#pragma omp taskwait
  }
  return 0;
}
