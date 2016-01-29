// Wrong example of using tasks without dependency.


#include<iostream>

using namespace std;

int main(int argc, char**argv){
#pragma omp parallel
  {
#pragma omp single nowait
    {
      for(int i = 0; i < 100; i++){
#pragma omp task
        {
          cout << "i = " <<  i << endl;
        }
      }
    }
#pragma omp taskwait
  }
  return 0;
}
