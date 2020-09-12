//nvcc -gencode=arch=compute_35,code=sm_35 -lcublas -ccbin "/work1/soft/intel/composer_xe_2013.0.079/bin/intel64/icc" -c this.cu
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<cuda_runtime.h>
#include<cublas_v2.h>
#include<cusolverDn.h>
extern "C" void cuda_cholesky_(int *n, double *a, double *b, int *info)
{
    int lwork;
    int n1,n2;
    double *d_work = NULL;
    double *d_A = NULL;
    double *d_B = NULL;
    int *devInfo = NULL;
    cusolverDnHandle_t cusolverH = NULL;

    n1=*n;
    n2=n1;

    cusolverDnCreate(&cusolverH);
    cudaMalloc((void**)&devInfo, sizeof(int));
    cudaMalloc((void**)&d_A , sizeof(double)*n1*n2);
    cudaMalloc((void**)&d_B , sizeof(double)*n1);
    cudaMemcpy(d_A, a, sizeof(double)*n1*n2, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, b, sizeof(double)*n1, cudaMemcpyHostToDevice);

    //query working space
    cusolverDnDpotrf_bufferSize(cusolverH,CUBLAS_FILL_MODE_UPPER,n1,d_A,n2,&lwork);
    cudaMalloc((void**)&d_work, sizeof(double)*lwork);

    cusolverDnDpotrf(cusolverH,CUBLAS_FILL_MODE_UPPER,n1,d_A,n2,d_work,lwork,devInfo);
    cudaMemcpy(info,devInfo,sizeof(int),cudaMemcpyDeviceToHost);

    if(*info == 0){
        cusolverDnDpotrs(cusolverH,CUBLAS_FILL_MODE_UPPER,n1,1,d_A,n2,d_B,n1,devInfo);
        cudaMemcpy(info,devInfo,sizeof(int),cudaMemcpyDeviceToHost);
        if(*info == 0) cudaMemcpy(b,d_B,n1*sizeof(double),cudaMemcpyDeviceToHost);
    }

    if(d_A) cudaFree(d_A);
    if(d_B) cudaFree(d_B);
    if(d_work) cudaFree(d_work);
    if(devInfo) cudaFree(devInfo);
    if(cusolverH) cusolverDnDestroy(cusolverH);
}
