#include<stdio.h>
#include"lapacke.h"
#include"cblas.h"

void init(double* A,int n)
{
for (int j = 0; j < 2; ++j) {
        A[j]=-1;
    }
    A[0]=abs(A[1])*2;

    for (int j = 3; j < n; ++j) {
        A[j]=0;

    }

    int deb=0;
    for (int i = 1; i < n-1; ++i) {
        for (int j = deb; j < deb+3; ++j) {
            A[i*n+j]=-1;
        }
        A[i*n+i]=abs(A[i*n + i-1]+A[i*n + i+1]);



        for (int j = deb+3; j < n; ++j) {
            A[i*n+j]=0;
        }
        deb++;


    }


    for (int j = n-2; j < n; ++j) {
        A[(n-1)*n+j]=-1;
    }
    A[(n-1)*n+n-1]=abs(A[(n-1)*n+n-1])*2;
    for (int j = 0; j < n-2; ++j) {
        A[(n-1)*n+j]=0;

    }
    
}
void set_GB_operator_colMajor_poisson1D(double* AB, int lab, int la, int kv){
  int ii, jj, kk;
  for (jj=0;jj<(la);jj++){
    kk = jj*(lab);
    
    AB[kk+ 0]=-1.0;
    AB[kk+ 1]=2.0;
    AB[kk+ 2]=-1.0;
  }
  AB[0]=0.0;
  
  
  AB[(lab)*(la)-1]=0.0;
}

int main(){
double *AB,*X,*Y,*A ; 
X=(double *) malloc(sizeof(double)*10);
AB = (double *) malloc(sizeof(double)*10*3);
Y=(double *) malloc(sizeof(double)*10);
A= (double *) malloc(sizeof(double)*10*10);
for(int i=0;i<10;i++){
    X[i]=i;
}
init(A,10);
set_GB_operator_colMajor_poisson1D(AB,  3, 10, 1);
cblas_dgbmv(CblasColMajor,CblasTrans,10,3,1,1,1.0,AB,10,X,1,0.0,Y,1);
for(int i=0;i<10;i++){
    X[i]=i;
    printf("%f ",Y[i]);
}


printf("AB=\n");
    for (int i = 0; i < 10; ++i) {
        printf("[");
        for (int j = 0; j < 10; ++j) {
            printf("%f,",A[i*10+j]);
        }
        printf("%f]\n",A[i*10+9]);
    }

return 0;
}
