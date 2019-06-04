#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include "tsne.h"

// Function that runs the Barnes-Hut implementation of t-SNE
int main() {

    // Define some variables
	int origN, N, D, no_dims, max_iter;
	double perplexity, theta, *data;
    int rand_seed = -1;

    // Read the parameters and the dataset
	if(TSNE::load_data(&data, &origN, &D, &no_dims, &theta, &perplexity, &rand_seed, &max_iter)) {

		// Make dummy landmarks
        N = origN;
        int* landmarks = (int*) malloc(N * sizeof(int));
        if(landmarks == NULL) { printf("Memory allocation failed!\n"); exit(1); }
        for(int n = 0; n < N; n++) landmarks[n] = n;

		// Now fire up the SNE implementation
		double* Y = (double*) malloc(N * no_dims * sizeof(double));
		double* costs = (double*) calloc(N, sizeof(double));
        if(Y == NULL || costs == NULL) { printf("Memory allocation failed!\n"); exit(1); }
		TSNE::run(data, N, D, Y, no_dims, perplexity, theta, rand_seed, false, max_iter, 250, 250);

		// Save the results
		TSNE::save_data(Y, landmarks, costs, N, no_dims);

        // Clean up the memory
		free(data); data = NULL;
		free(Y); Y = NULL;
		free(costs); costs = NULL;
		free(landmarks); landmarks = NULL;
    }
}
