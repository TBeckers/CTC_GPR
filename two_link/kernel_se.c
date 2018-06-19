/* kernel
* Copyright (C) Technical University of Munich
* Written by Thomas Beckers, ITR, Technical University of Munich, 2018
 */
 
#include "mex.h"
#include <math.h>

void Product(double *x, double *y, double *z, mwSize n,mwSize m,mwSize dim,double sd,double *l)
{
    mwSize i,j,k;
    mwSize count=0;
    double dis;

    for (i=0; i<n; i++) {
    	for (j=0; j<m; j++) {
            dis=0;
            for (k=0; k<dim; k++) {
                dis+=pow((x[i*dim+k]- y[j*dim+k])/l[k],2);
                /* printf("x: %f, y: %f, k: %d\n",x[i*dim+k],y[j*dim+k],k); */
            }
            z[count] = sd*exp(-1*dis/2);
            count++;
        }
    }
}

void Product_same(double *x, double *y, double *z, mwSize n,mwSize m,mwSize dim,double sd,double *l)
{
    mwSize i,j,k;
    mwSize count=0;
    double dis;

    for (i=0; i<n; i++) {
    	for (j=i; j<m; j++) {
            dis=0;
            for (k=0; k<dim; k++) {
                dis+=pow((x[i*dim+k]- y[j*dim+k])/l[k],2);
                /* printf("x: %f, y: %f, k: %d\n",x[i*dim+k],y[j*dim+k],k); */
            }
            z[count] = sd*exp(-1*dis/2);
            z[(count%n)*n+(int)(count/n)]=z[count];
            count++;
        }
        count+=i+1;
    }
}

void Derivation(double *x, double *y, double *z, mwSize n,mwSize m,mwSize dim,double sd,double *l,int d)
{
    mwSize i,j,k;
    mwSize count=0;
    double dis,dis1;

    for (i=0; i<n; i++) {
    	for (j=0; j<m; j++) {
            dis=0;
            for (k=0; k<dim; k++) {
                dis+=pow((x[i*dim+k]- y[j*dim+k])/l[k],2);
                /* printf("x: %f, y: %f, k: %d\n",x[i*dim+k],y[j*dim+k],k); */
            }
            if (d<dim)
                dis1=pow((x[i*dim+d]- y[j*dim+d])/l[d],2);
            else
                dis1=2;
            z[count] = sd*exp(-1*dis/2)*dis1;
            count++;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *inMatrix1,*inMatrix2,*outMatrix,*l;               
    size_t ncols1,ncols2,nrows1,nrows2,nrowsl; 
    double sd;
    int d;
    
    inMatrix1 = mxGetPr(prhs[0]);
    ncols1 = mxGetN(prhs[0]);
    nrows1 = mxGetM(prhs[0]);
    
    inMatrix2 = mxGetPr(prhs[1]);
    ncols2 = mxGetN(prhs[1]);
    nrows2 = mxGetM(prhs[1]);
    
    sd = pow(mxGetScalar(prhs[2]),2);
    l = mxGetPr(prhs[3]);
    nrowsl = mxGetM(prhs[3]);
    
    if ((nrows1!=nrows2) || (nrowsl!=nrows1))
    {
        mexErrMsgTxt("Matrix have not the same dimension");
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix((mwSize)ncols2,(mwSize)ncols1,mxREAL);
    
        outMatrix = mxGetPr(plhs[0]);
        if (nrhs==5)
        {
            d=mxGetScalar(prhs[4])-1;
            if (d<0)
                    Product_same(inMatrix1,inMatrix2,outMatrix,(mwSize)ncols1,(mwSize)ncols2,(mwSize)nrows1,sd,l);
            else
            Derivation(inMatrix1,inMatrix2,outMatrix,(mwSize)ncols1,(mwSize)ncols2,(mwSize)nrows1,sd,l,d);
        }
        else
        {
            Product(inMatrix1,inMatrix2,outMatrix,(mwSize)ncols1,(mwSize)ncols2,(mwSize)nrows1,sd,l);
        }
    }
}
