cimport pytsne
cimport numpy as np


cdef class tsne:
    @staticmethod
    def run(np.ndarray[double_t, ndim=2, mode="c"] X,
        int no_dims, # number of output dims
        double perplexity,
        double theta,
        int rand_seed,
        int skip_random_init,
        int max_iter=1000,
        int stop_lying_iter=250,
        int mom_switch_iter=250):

        N = X.shape[0]
        D = X.shape[1]
        cdef np.ndarray[double_t, ndim=2, mode="c"] Y = np.ndarray((N, no_dims))

        tsne_run(&X[0, 0], N, D, &Y[0, 0], no_dims, perplexity, theta, rand_seed,
            skip_random_init, max_iter, stop_lying_iter, mom_switch_iter)

        return Y
