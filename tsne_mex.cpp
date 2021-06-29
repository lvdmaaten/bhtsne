// Mex function that runs Laurens van der Maaten's Barnes-Hut implementation of t-SNE
// (C) Patrick Fletcher 2018

#include <vector>

#include "mex.h"
#include "tsne.h"


/* Input Arguments */
#define	DATA	prhs[0]
#define	NO_DIMS	prhs[1]
#define	THETA	prhs[2]
#define	PERPLEXITY	prhs[3]
#define	MAX_ITER	prhs[4]
#define	RAND_SEED	prhs[5]
#define	INIT_Y	prhs[6]
#define	STOP_LYING_ITER	prhs[7]
#define	MOM_SWITCH_ITER	prhs[8]

/* Output Arguments */
#define	MAPPED_DATA	plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
	
    /* Check for proper number of arguments */
    if (nrhs != 9) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidNumInputs",
                "9 input arguments required."); 
    }
	
	if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:maxlhs",
                "Too many output arguments."); 
    }
	
	/* Input Type checking */
	if (!mxIsSingle(DATA) && !mxIsDouble(DATA)) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidT",
                "bh_tsne requires that DATA be single or double"); 
    } 
    /* Check the dimensions of NO_DIMS, THETA, PERPLEXITY, MAX_ITER, and RAND_SEED.  They all should be scalar */ 
	if (!mxIsScalar(NO_DIMS) || (!mxIsSingle(NO_DIMS) && !mxIsDouble(NO_DIMS)) ) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidInput",
                "BH_TSNE requires that NO_DIMS be a scalar single or double"); 
    } 
	if (!mxIsScalar(THETA) || (!mxIsSingle(THETA) && !mxIsDouble(THETA)) ) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidInput",
                "BH_TSNE requires that THETA be a scalar single or double"); 
    } 
	if (!mxIsScalar(PERPLEXITY) || (!mxIsSingle(PERPLEXITY) && !mxIsDouble(PERPLEXITY)) ) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidInput",
                "BH_TSNE requires that PERPLEXITY be a scalar single or double"); 
    } 
	if (!mxIsScalar(MAX_ITER) || (!mxIsSingle(MAX_ITER) && !mxIsDouble(MAX_ITER)) ) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidInput",
                "BH_TSNE requires that MAX_ITER be a scalar single or double"); 
    } 
	if (!mxIsScalar(RAND_SEED) || (!mxIsSingle(RAND_SEED) && !mxIsDouble(RAND_SEED)) ) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidInput",
                "BH_TSNE requires that RAND_SEED be a scalar single or double");
    } 
	if (!mxIsSingle(INIT_Y) && !mxIsDouble(INIT_Y)) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidT",
                "bh_tsne requires that INIT_Y be single or double"); 
    } 
	if (!mxIsScalar(STOP_LYING_ITER) || (!mxIsSingle(STOP_LYING_ITER) && !mxIsDouble(STOP_LYING_ITER)) ) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidInput",
                "BH_TSNE requires that STOP_LYING_ITER be a scalar single or double");
    } 
	if (!mxIsScalar(MOM_SWITCH_ITER) || (!mxIsSingle(MOM_SWITCH_ITER) && !mxIsDouble(MOM_SWITCH_ITER)) ) { 
	    mexErrMsgIdAndTxt( "MATLAB:bh_tsne:invalidInput",
                "BH_TSNE requires that MOM_SWITCH_ITER be a scalar single or double");
    } 
	
	//dimensions of data: we are passing in X' (rows=vars, cols=observations) to get row major order in the data vector
    int D = (int)mxGetM(DATA); 
    int N = (int)mxGetN(DATA);
	
    mexPrintf("N: %d, D: %d\n", N, D);
	
	std::vector<double> X ( static_cast<double *>(mxGetData(DATA)),  static_cast<double *>(mxGetData(DATA)) + mxGetNumberOfElements(DATA) );  
	int no_dims=static_cast<int>(mxGetScalar(NO_DIMS));
	double theta=static_cast<double>(mxGetScalar(THETA));
	double perplexity=static_cast<double>(mxGetScalar(PERPLEXITY));
	int max_iter=static_cast<int>(mxGetScalar(MAX_ITER));
	int rand_seed=static_cast<int>(mxGetScalar(RAND_SEED));
	int stop_lying_iter=static_cast<int>(mxGetScalar(STOP_LYING_ITER));
	int mom_switch_iter=static_cast<int>(mxGetScalar(MOM_SWITCH_ITER));
	
	size_t y_length = no_dims * N;
	std::vector<double> Y(y_length);

	bool skip_random_init = !mxIsEmpty(INIT_Y);
	if (skip_random_init) {
		std::copy(static_cast<double *>(mxGetData(INIT_Y)),  static_cast<double *>(mxGetData(INIT_Y)) + mxGetNumberOfElements(INIT_Y), Y.begin() );
	}
	
	//run tsne
	//TSNE tsne();
    TSNE* tsne = new TSNE();
	tsne->run(X.data(), N, D, Y.data(), no_dims, perplexity, theta, rand_seed, skip_random_init, max_iter, stop_lying_iter, mom_switch_iter);
    //void run(double* X, int N, int D, double* Y, int no_dims, double perplexity, double theta, int rand_seed,
    //         bool skip_random_init, int max_iter=1000, int stop_lying_iter=250, int mom_switch_iter=250);
	
	//send the result back to Matlab
    MAPPED_DATA = mxCreateDoubleMatrix( (mwSize)no_dims, (mwSize)N, mxREAL); 
	std::copy(Y.begin(), Y.end(), (double *)mxGetData(MAPPED_DATA));

    delete(tsne);
	return;
}

