#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>

#include <bhtsne/tsne.h>
#include "CommandlineOptions.h"

using namespace bhtsne;

// Function that runs the Barnes-Hut implementation of t-SNE
int main(int argc, char * argv[]) {

    // Define some variables
	int D, no_dims, max_iter;
	double perplexity, theta, *data;
    int rand_seed = -1;
    TSNE* tsne = new TSNE();

    applyCommandlineOptions(*tsne, argc, argv);

    // Read the parameters and the dataset
	//if (tsne->loadLegacy("C:/Users/dajg1/Documents/HPI/cpp/bhtsne/build/build/Releasedata.dat")) {
	if (tsne->load_data(&data, &D, &no_dims, &theta, &perplexity, &rand_seed, &max_iter)) {
		//tsne->setNumberOfSamples(5000);
		/*tsne->run();
		tsne->saveLegacy();/**/

		rand_seed = 0;
		perplexity = 50.0;
		theta = 0.2;
		max_iter = 1000;
		// Make dummy landmarks
        int* landmarks = (int*) malloc(tsne->getNumberOfSamples() * sizeof(int));
        if(landmarks == nullptr) {
			printf("Memory allocation failed!\n");
			exit(1);
		}
        for(unsigned int n = 0; n < tsne->getNumberOfSamples(); n++)
			landmarks[n] = n;

		// Now fire up the SNE implementation
		double* Y = (double*) malloc(tsne->getNumberOfSamples() * no_dims * sizeof(double));
		double* costs = (double*) calloc(tsne->getNumberOfSamples(), sizeof(double));
        if(Y == nullptr || costs == nullptr) {
			printf("Memory allocation failed!\n");
			exit(1);
		}
		tsne->run(data, D, Y, no_dims, perplexity, theta, rand_seed, false, max_iter);

		// Save the results
		tsne->save_data(Y, landmarks, costs, tsne->getNumberOfSamples(), no_dims);

        // Clean up the memory
		//free(data); data = nullptr;
		free(Y); Y = nullptr;
		free(costs); costs = nullptr;
		free(landmarks); landmarks = nullptr;/**/
    }
    delete(tsne);
}
