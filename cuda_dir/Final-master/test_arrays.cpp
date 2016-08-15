#include <iostream>

enum Containers {matrix, vector};


int* function(int* pt, int i, Containers container){
/*int c=0;
for(int i=0; i<5; i++){
c+=ARRAY[i];
}*/
if(container == vector){
std::cout<<"it worked"<<std::endl;
}
int* p;
p=&pt[i];
return p;
//return c;
}

int main(){

int N=2;
int M=3;
int L=4;

double foo [N];

for(int j=0; j<N; j++){

foo[j]=2;
}

for(int i=0; i<N; i++){
std::cout<<foo[i]<<std::endl;
}

double Foo[N][M][L];

for(int k=0; (k < N); k++)
   {
      
      for(int l=0; (l < M); l++)
      {
       
        for(int m=0; m<L; m++){
                Foo[k][l][m] = k*l*m;
                }
        }
   }



/*
for(int x=0; (x < N); x++)
   {
std::cout<<"x"<<std::endl;
      for(int y=0; (y < M); y++)
      {
std::cout<<"y"<<std::endl;
        for(int z=0; z<L; z++){
                std::cout<<Foo[x][y][z]<<std::endl;
                }
        }
   }
*/



int* pointer;
int* p;
//pointer= new int[5];

int ARRAY[5];
p=ARRAY;
for(int it=0; it<5; it++){
p[it]=it;
}

Containers vec;
vec = vector;
*function(p, 4, vec )=100;

int test;

test=(*function(p, 3, vec)) * (*function(p, 4, vec));
//int c= function(ARRAY);

std::cout<<sizeof(ARRAY)/sizeof(*ARRAY)<<std::endl;

//delete[] pointer;

return 0;
}
