cdef extern:
    ctypedef double double_t
    void tsne_run(double* X, int N, int D, double* Y,
      int no_dims, double perplexity, double theta, int rand_seed,
      int skip_random_init, int max_iter, int stop_lying_iter,
      int mom_switch_iter) nogil
