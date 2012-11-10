#include <math.h>
#include "mex.h"
#define	T_IN	prhs[0]
#define	Y_IN	prhs[1]
#define	YP_OUT	plhs[0]

#if !defined(MAX)
	#define MAX(A, B)       ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
	#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#endif




static void yprime(double yp[],double *t,double y[])
{
	yp[0]=y[1];
	yp[1]=-(2943*y[2]*cos(y[0]))/(300*y[2]*y[2]*cos(y[0])*cos(y[0]) + 300*y[2]*y[2]*sin(y[0])*sin(y[0]) + 325);
	yp[2]=y[3];
	yp[3]=(100*y[1]*y[1]*y[2]*cos(y[0])*cos(y[0]) - 981*sin(y[0]) + 100*y[1]*y[1]*y[2]*sin(y[0])*sin(y[0]))/(100*cos(y[0])*cos(y[0]) + 100*sin(y[0])*sin(y[0]));
}



void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
	double *yp;
	double *t,*y;
	unsigned int m,n;

	if (nrhs != 2) {
		mexErrMsgTxt("Two input arguments required.");
	} else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments.");
	}
	/* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
	m = mxGetM(Y_IN);
	n = mxGetN(Y_IN);
	if (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) ||
		(MAX(m,n) != 4) || (MIN(m,n) != 1)) {
		mexErrMsgTxt("mexODEFun requires that Y be a 4 x 1vector.");
	}
	/* Create a matrix for the return argument */ 
	YP_OUT = mxCreateDoubleMatrix(m, n, mxREAL);
	/* Assign pointers to the various parameters */ 
	yp = mxGetPr(YP_OUT);
	t = mxGetPr(T_IN);
	y = mxGetPr(Y_IN);
	/* Do the actual computations in a subroutine */ 
	yprime(yp,t,y);
	/* Return YP_OUT to MATLAB environment */ 
	return;
}
