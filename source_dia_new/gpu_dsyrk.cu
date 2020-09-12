//nvcc -gencode=arch=compute_35,code=sm_35 -lcublas -ccbin "/work1/soft/intel/composer_xe_2013.0.079/bin/intel64/icc" -c this.cu
#include<stdio.h>
#include<stdlib.h>
#include<cuda_runtime.h>
#include<cublas_v2.h>
// u = 'u', t = 't', A(nt,nw), C(nw,nw), C=A'*A
extern "C" void gpu_dsyrk_(char *u, char *n, int *nw, int *nt, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc)
{
  int N1,N2;
  double *dA,*dC;
  cublasHandle_t handle;
  N1=*nt;
  N2=*nw;
  cudaMalloc((void **)&dA,N1*N2*sizeof(double));
  cudaMalloc((void **)&dC,N2*N2*sizeof(double));
  cublasCreate(&handle);
  cudaMemcpy(dA,a,N1*N2*sizeof(double),cudaMemcpyHostToDevice);
  cublasDsyrk(handle,CUBLAS_FILL_MODE_UPPER,CUBLAS_OP_N,N2,N1,alpha,dA,N2,beta,dC,N2);
  cudaMemcpy(c,dC,N2*N2*sizeof(double),cudaMemcpyDeviceToHost);
  cudaFree(dA);
  cudaFree(dC);
  cublasDestroy(handle);
}
