/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */


#include <bhtsne/tsne.h>


#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cassert>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <stdexcept>

#include <bhtsne/sptree.h>
#include <bhtsne/vptree.h>

#include <bhtsne/bhtsne-version.h> // includes BHTSNE_VERSION macro


using namespace bhtsne;


TSNE::TSNE()
    : m_perplexity(50.0)
    , m_gradientAccuracy(0.2)
    , m_iterations(1000)
    , m_outputDimensions(2)
    , m_inputDimensions(0)
    , m_dataSize(0)
    , m_outputFile("result")
    , m_gen(0) //default seed
{
}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
void TSNE::computeGradient(unsigned int* inp_row_P, unsigned int* inp_col_P,
    double* inp_val_P, double* Y, int D, double* dC, double theta) {

    // Construct space-partitioning tree on current map
    auto tree = SPTree(D, Y, m_dataSize);

    // Compute all terms required for t-SNE gradient
    double sum_Q = .0;
    double* pos_f = (double*) calloc(m_dataSize * D, sizeof(double));
    double* neg_f = (double*) calloc(m_dataSize * D, sizeof(double));
    if(pos_f == nullptr || neg_f == nullptr) { printf("Memory allocation failed!\n"); exit(1); }
    tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, m_dataSize, pos_f);
    for(int n = 0; n < m_dataSize; n++)
        tree.computeNonEdgeForces(n, theta, neg_f + n * D, &sum_Q);

    // Compute final t-SNE gradient
    for(int i = 0; i < m_dataSize * D; i++) {
        dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
    }
    free(pos_f);
    free(neg_f);
}

// Compute gradient of the t-SNE cost function (exact)
void TSNE::computeExactGradient(double* P, double* Y /*m_result*/, std::vector<double> & gradients)
{
    assert(gradients.size() == m_dataSize * m_outputDimensions);

    auto D = m_outputDimensions;
	// Make sure the current gradient contains zeros
    std::fill(gradients.begin(), gradients.end(), 0.0);

    // Compute the squared Euclidean distance matrix
    auto distances = computeSquaredEuclideanDistance(Y); //no params; return distances
    assert(distances.size() == m_dataSize * m_dataSize);

    // Compute Q-matrix and normalization sum
    // https://en.wikipedia.org/wiki/Q-matrix
    auto Q = std::vector<double>(m_dataSize * m_dataSize);
    double sum_Q = .0;
    for(int n = 0; n < m_dataSize; n++) {
    	for(int m = 0; m < m_dataSize; m++) {
            auto matrixIndex = n*m_dataSize + m;
            if(n != m) {
                Q[matrixIndex] = 1 / (1 + distances[matrixIndex]);
                sum_Q += Q[matrixIndex];
            }
        }
    }

	// Perform the computation of the gradient
	for(int n = 0; n < m_dataSize; n++) {
    	for(int m = 0; m < m_dataSize; m++) {
            if(n != m) {
                auto matrixIndex = n*m_dataSize + m;
                double mult = (P[matrixIndex] - (Q[matrixIndex] / sum_Q)) * Q[matrixIndex];
                for(int d = 0; d < D; d++) {
                    gradients[n*m_outputDimensions + d] += (Y[n*m_outputDimensions + d] - Y[m*m_outputDimensions + d]) * mult;
                }
            }
		}
	}
}


// Evaluate t-SNE cost function (exactly)
double TSNE::evaluateError(double* P, double* Y) 
{
    // Compute the squared Euclidean distance matrix
    auto Q = std::vector<double>(m_dataSize * m_dataSize);
    auto distances = computeSquaredEuclideanDistance(Y);
    auto DD = distances.data();

    // Compute Q-matrix and normalization sum
    double sum_Q = DBL_MIN;
    for(int n = 0; n < m_dataSize; n++) {
    	for(int m = 0; m < m_dataSize; m++) {
            auto matrixIndex = n*m_dataSize + m;
            if(n != m) {
                Q[matrixIndex] = 1 / (1 + DD[matrixIndex]);
                sum_Q += Q[matrixIndex];
            }
            else Q[matrixIndex] = DBL_MIN;
        }
    }
    for(auto& each : Q)
        each /= sum_Q;

    // Sum t-SNE error
    double C = .0;
	for(int n = 0; n < m_dataSize * m_dataSize; n++) {
        C += P[n] * log((P[n] + FLT_MIN) / (Q[n] + FLT_MIN));
	}

	return C;
}

// Evaluate t-SNE cost function (approximately)
double TSNE::evaluateError(unsigned int* row_P, unsigned int* col_P,
    double* val_P, double* Y, int D, double theta) {

    // Get estimate of normalization term
    auto tree = SPTree(D, Y, m_dataSize);
    double* buff = (double*) calloc(D, sizeof(double));
    double sum_Q = .0;
    for(int n = 0; n < m_dataSize; n++)
        tree.computeNonEdgeForces(n, theta, buff, &sum_Q);

    // Loop over all edges to compute t-SNE error
    int ind1, ind2;
    double C = .0, Q;
    for(int n = 0; n < m_dataSize; n++) {
        ind1 = n * D;
        for(int i = row_P[n]; i < row_P[n + 1]; i++) {
            Q = .0;
            ind2 = col_P[i] * D;
            for(int d = 0; d < D; d++) buff[d]  = Y[ind1 + d];
            for(int d = 0; d < D; d++) buff[d] -= Y[ind2 + d];
            for(int d = 0; d < D; d++) Q += buff[d] * buff[d];
            Q = (1.0 / (1.0 + Q)) / sum_Q;
            C += val_P[i] * log((val_P[i] + FLT_MIN) / (Q + FLT_MIN));
        }
    }

    // Clean up memory
    free(buff);
    return C;
}

// Symmetrizes a sparse matrix
void TSNE::symmetrizeMatrix(unsigned int** _row_P, unsigned int** _col_P, double** _val_P, int N) {

    // Get sparse matrix
    unsigned int* row_P = *_row_P;
    unsigned int* col_P = *_col_P;
    double* val_P = *_val_P;

    // Count number of elements and row counts of symmetric matrix
    auto row_counts = std::vector<int>(N, 0);
    for(int n = 0; n < N; n++) {
        for(int i = row_P[n]; i < row_P[n + 1]; i++) {

            // Check whether element (col_P[i], n) is present
            bool present = false;
            for(int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
                if(col_P[m] == n) present = true;
            }
            if(present) row_counts[n]++;
            else {
                row_counts[n]++;
                row_counts[col_P[i]]++;
            }
        }
    }
    int no_elem = 0;
    for(int n = 0; n < N; n++) no_elem += row_counts[n];

    // Allocate memory for symmetrized matrix
    unsigned int* sym_row_P = (unsigned int*) malloc((N + 1) * sizeof(unsigned int));
    unsigned int* sym_col_P = (unsigned int*) malloc(no_elem * sizeof(unsigned int));
    double* sym_val_P = (double*) malloc(no_elem * sizeof(double));
    if(sym_row_P == nullptr || sym_col_P == nullptr || sym_val_P == nullptr) {
        printf("Memory allocation failed!\n");
        exit(1);
    }

    // Construct new row indices for symmetric matrix
    sym_row_P[0] = 0;
    for(int n = 0; n < N; n++) sym_row_P[n + 1] = sym_row_P[n] + (unsigned int) row_counts[n];

    // Fill the result matrix
    auto offset = std::vector<int>(N, 0);
    for(int n = 0; n < N; n++) {
        for(unsigned int i = row_P[n]; i < row_P[n + 1]; i++) { // considering element(n, col_P[i])

            // Check whether element (col_P[i], n) is present
            bool present = false;
            for(unsigned int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
                if(col_P[m] == n) {
                    present = true;
                    if(n <= col_P[i]) { // make sure we do not add elements twice
                        sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                        sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                        sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
                        sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
                    }
                }
            }

            // If (col_P[i], n) is not present, there is no addition involved
            if(!present) {
                sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
                sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
            }

            // Update offsets
            if(!present || (present && n <= col_P[i])) {
                offset[n]++;
                if(col_P[i] != n)
                    offset[col_P[i]]++;
            }
        }
    }

    // Divide the result by two
    for(int i = 0; i < no_elem; i++) sym_val_P[i] /= 2.0;

    // Return symmetrized matrices
    free(*_row_P); *_row_P = sym_row_P;
    free(*_col_P); *_col_P = sym_col_P;
    free(*_val_P); *_val_P = sym_val_P;
}

// with mean zero and standard deviation one
double bhtsne::TSNE::gaussNumber()
{
    // Knuth, Art of Computer Programming vol v2, Section 3.4.1, Algorithm P (p.117)

    double S, V1, V2;
    // executed 1.27 times on the average
    do
    {
        // V1, V2 uniformly distributed between -1 and +l.
        V1 = 2.0 * static_cast<double>(m_gen()) / m_gen.max() - 1.0;
        V2 = 2.0 * static_cast<double>(m_gen()) / m_gen.max() - 1.0;
        S = V1*V1 + V2*V2;
    } while (S >= 1);

    auto X1 = V1 * std::sqrt(-2 * std::log(S) / S);
    //same can be done for X2

    return X1;
}


// Makes data zero-mean
void TSNE::zeroMean(double* X, int N, int D) {

	// Compute data mean
	auto mean = std::vector<double>(D, 0.0);
    int nD = 0;
	for(int n = 0; n < N; n++) {
		for(int d = 0; d < D; d++) {
			mean[d] += X[nD + d];
		}
        nD += D;
	}
	for(int d = 0; d < D; d++) {
		mean[d] /= (double) N;
	}

	// Subtract data mean
    nD = 0;
	for(int n = 0; n < N; n++) {
		for(int d = 0; d < D; d++) {
			X[nD + d] -= mean[d];
		}
        nD += D;
	}
}


void TSNE::setRandomSeed(int seed)
{
    if (seed >= 0) {
        std::cout << "Using random seed: " << seed << std::endl;
        m_gen.seed(seed);

    }
    else {
        std::cout << "Using current time as random seed..." << std::endl;
        m_gen.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }
}

double TSNE::perplexity() const
{
	return m_perplexity;
}

void TSNE::setPerplexity(double perplexity)
{
	m_perplexity = perplexity;
    if (m_perplexity < 2.0)
    {
        std::cerr << "perplexity has to be at least 2.0, setting perplexity to 2.0" << std::endl;
        m_perplexity = 2.0;
    }
}

double TSNE::gradientAccuracy() const
{
	return m_gradientAccuracy;
}

void TSNE::setGradientAccuracy(double accuracy)
{
	m_gradientAccuracy = accuracy;
}

unsigned int TSNE::iterations() const
{
	return m_iterations;
}

void TSNE::setIterations(unsigned int iterations)
{
	m_iterations = iterations;
}

unsigned int TSNE::outputDimensions() const
{
	return m_outputDimensions;
}

void TSNE::setOutputDimensions(unsigned int dimensions)
{
	m_outputDimensions = dimensions;
}

unsigned int TSNE::inputDimensions() const
{
	return m_inputDimensions;
}

unsigned int TSNE::dataSize() const
{
	return m_dataSize;
}

void TSNE::setDataSize(unsigned int value)
{
	m_dataSize = value;
}

std::string TSNE::outputFile() const
{
	return m_outputFile;
}

void TSNE::setOutputFile(const std::string& file)
{
	m_outputFile = file;
}


//load methods--------------------------------------------------------------------------------------

bool bhtsne::TSNE::loadFromStream(std::istream & stream)
{
    auto separator = ',';

    //read data points
    auto line = std::string();
    bool first = true;
    while (std::getline(stream, line))
    {
        auto iss = std::istringstream(line);
        auto element = std::string();

        auto point = std::vector<double>();
        point.reserve(m_inputDimensions);

        //read values of data point
        while (std::getline(iss, element, separator))
        {
            point.push_back(std::stod(element));
        }

        //set dimensionality
        if (first)
        {
            first = false;
            m_inputDimensions = point.size();
        }

        //fail if inconsistent dimensionality
        if (m_inputDimensions != point.size() || m_inputDimensions == 0)
        {
            m_inputDimensions = 0;
            m_data.clear();
            return false;
        }

        m_data.emplace_back(std::move(point));
    }

    if (m_data.empty())
        return false;

    m_dataSize = m_data.size();

    return true;
}

bool TSNE::loadLegacy(const std::string & file)
{
	auto f = std::ifstream(file, std::ios::binary);
	if (!f.is_open())
    {
        std::cerr << "Could not open " << file << std::endl;
        return false;
    }

    //read params
	f.read(reinterpret_cast<char*>(&m_dataSize), sizeof(m_dataSize));
	f.read(reinterpret_cast<char*>(&m_inputDimensions), sizeof(m_inputDimensions));
	f.read(reinterpret_cast<char*>(&m_gradientAccuracy), sizeof(m_gradientAccuracy));
	f.read(reinterpret_cast<char*>(&m_perplexity), sizeof(m_perplexity));
	f.read(reinterpret_cast<char*>(&m_outputDimensions), sizeof(m_outputDimensions));
	f.read(reinterpret_cast<char*>(&m_iterations), sizeof(m_iterations));

    //read data
	for (size_t i = 0; i < m_dataSize; ++i) {
		auto point = std::vector<double>(m_inputDimensions);
		f.read(reinterpret_cast<char*>(point.data()), sizeof(double) * m_inputDimensions);
        m_data.emplace_back(std::move(point));
	}

    //read seed
	if (!f.eof()) {
        int seed;
		f.read(reinterpret_cast<char*>(&seed), sizeof(seed));
        setRandomSeed(seed);
	}

	return true;
}

bool TSNE::loadCSV(const std::string & file)
{
	auto f = std::ifstream(file);
	if (!f.is_open())
    {
        std::cerr << "Could not open " << file << std::endl;
        return false;
    }

    return loadFromStream(f);
}

bool TSNE::loadTSNE(const std::string & file)
{
	auto f = std::ifstream(file, std::ios::binary);
	if (!f.is_open())
    {
        std::cerr << "Could not open " << file << std::endl;
        return false;
    }

	f.read(reinterpret_cast<char*>(&m_dataSize), sizeof(m_dataSize));
	f.read(reinterpret_cast<char*>(&m_inputDimensions), sizeof(m_inputDimensions));

	for (size_t i = 0; i < m_dataSize; ++i)
    {
		auto point = std::vector<double>(m_inputDimensions);
		f.read(reinterpret_cast<char*>(point.data()), sizeof(double) * m_inputDimensions);
		m_data.emplace_back(std::move(point));
	}

	return true;
}

bool TSNE::loadCin()
{
    return loadFromStream(std::cin);
}


//run method---------------------------------------------------------------------------------------
void TSNE::run()
{
    if (m_dataSize - 1 < 3 * m_perplexity) {
        auto message = "perplexity (perplexity=" + std::to_string(m_perplexity) +
            ") has to be smaller than a third of the dataSize (dataSize=" + std::to_string(m_dataSize) + ")";
        std::cerr << message << std::endl;
        throw std::invalid_argument(message);
    }

    std::cout << "Using: " << std::endl
        << "data size " << m_dataSize << std::endl
        << "in dimensions " << m_inputDimensions << std::endl
        << "out dimensions " << m_outputDimensions << std::endl
        << "perplexity " << m_perplexity << std::endl
        << "gradient accuracy " << m_gradientAccuracy << std::endl;

    if (m_gradientAccuracy == 0.0)
        runExact();
    else
        runApproximation();
}


void TSNE::runApproximation()
{
	int stop_lying_iter = 250;
	int mom_switch_iter = 250;

	// Set learning parameters
	double momentum = .5, final_momentum = .8;
	double eta = 200.0;

	// Allocate some memory
	double* dY = (double*)malloc(m_dataSize * m_outputDimensions * sizeof(double));
	auto uY = std::vector<double>(m_dataSize * m_outputDimensions, 0.0);
	auto gains = std::vector<double>(m_dataSize * m_outputDimensions, 1.0);
	if (dY == nullptr) {
		std::cout << "Memory allocation failed!" << std::endl;
		exit(1);
	}

	// Normalize input data (to prevent numerical problems)
	std::cout << "Computing input similarities..." << std::endl;
	zeroMean(m_data, m_inputDimensions);
	// TODO: extract normalization function for vector
	double max_X = .0;
	for (auto & point : m_data)
	{
		for (auto & val : point)
		{
			if (std::fabs(val) > max_X) {
				max_X = std::fabs(val);
			}
		}
	}
	for (auto & point : m_data)
	{
		for (auto & val : point)
		{
			val /= max_X;
		}
	}

	// Compute input similarities for exact t-SNE
	unsigned int* row_P; unsigned int* col_P;
    std::vector<double> modern_val_P;

	// Compute asymmetric pairwise input similarities
	computeGaussianPerplexity(&row_P, &col_P, modern_val_P);
    auto val_P = modern_val_P.data();

	// Symmetrize input similarities
	symmetrizeMatrix(&row_P, &col_P, &val_P, m_dataSize);
	//normalize val_P so that sum of all val = 1
	double sum_P = .0;
	for (int i = 0; i < row_P[m_dataSize]; i++) sum_P += val_P[i];
	for (int i = 0; i < row_P[m_dataSize]; i++) val_P[i] /= sum_P;


	// Lie about the P-values
	for (int i = 0; i < row_P[m_dataSize]; i++)
		val_P[i] *= 12.0;

	double* Y = (double*)malloc(m_dataSize * m_outputDimensions * sizeof(double));
	// Initialize solution (randomly)
	for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
	    Y[i] = gaussNumber() * .0001;

	// Perform main training loop
    std::cout << " Input similarities computed. Learning embedding..." << std::endl;

	for (int iter = 0; iter < m_iterations; iter++) {

		// Compute approximate gradient
        computeGradient(row_P, col_P, val_P, Y, m_outputDimensions, dY, m_gradientAccuracy);

		// Update gains
		for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
			gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
		for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
			gains[i] = std::max(gains[i], 0.1);

		// Perform gradient update (with momentum and gains)
		for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
			uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
		for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
			Y[i] = Y[i] + uY[i];

		// Make solution zero-mean
		zeroMean(Y, m_dataSize, m_outputDimensions);

		// Stop lying about the P-values after a while, and switch momentum
		if (iter == stop_lying_iter) {
			for (int i = 0; i < row_P[m_dataSize]; i++)
				val_P[i] /= 12.0;
		}
		if (iter == mom_switch_iter) momentum = final_momentum;

		// Print out progress
		if (iter > 0 && (iter % 50 == 0 || iter == m_iterations - 1)) {
			double C = .0;
			// doing approximate computation here!
			C = evaluateError(row_P, col_P, val_P, Y, m_outputDimensions, m_gradientAccuracy);

			std::cout << "Iteration " << (iter + 1) << ": error is " << C << std::endl;
		}
	}

	// Clean up memory
	free(dY);
	free(row_P); row_P = nullptr;
	free(col_P); col_P = nullptr;

	size_t offset = 0;
	for (size_t i = 0; i < m_dataSize; ++i) {
		auto point = std::vector<double>();
		for (size_t j = 0; j < m_outputDimensions; ++j) {
			point.push_back(Y[offset++]);
		}
		m_result.push_back(point);
	}
}


void TSNE::runExact()
{
    int stop_lying_iter = 250;
    int mom_switch_iter = 250;

    // Set learning parameters
    double momentum = .5, final_momentum = .8;
    double eta = 200.0;

    // Normalize input data (to prevent numerical problems)
    std::cout << "Computing input similarities..." << std::endl;
    zeroMean(m_data, m_inputDimensions);
    // TODO: extract normalization function for vector
    double max_X = .0;
    for (auto & point : m_data)
    {
        for (auto & val : point)
        {
            if (std::fabs(val) > max_X) {
                max_X = std::fabs(val);
            }
        }
    }
    for (auto & point : m_data)
    {
        for (auto & val : point)
        {
            val /= max_X;
        }
    }

    // Compute input similarities for exact t-SNE
    double* P = (double*)malloc(m_dataSize * m_dataSize * sizeof(double));
    if (P == nullptr)
    {
        std::cout << "Memory allocation failed!" << std::endl;
        exit(1);
    }
    computeGaussianPerplexity(P);

    // Symmetrize input similarities
    std::cout << "Symmetrizing..." << std::endl;
    int nN = 0;
    for (int n = 0; n < m_dataSize; n++) {
        int mN = (n + 1) * m_dataSize;
        for (int m = n + 1; m < m_dataSize; m++) {
            P[nN + m] += P[mN + n];
            P[mN + n] = P[nN + m];
            mN += m_dataSize;
        }
        nN += m_dataSize;
    }
    double sum_P = .0;
    for (int i = 0; i < m_dataSize * m_dataSize; i++) sum_P += P[i];
    for (int i = 0; i < m_dataSize * m_dataSize; i++) P[i] /= sum_P;
    

    // Lie about the P-values
    for (int i = 0; i < m_dataSize * m_dataSize; i++)
        P[i] *= 12.0;

    double* Y = (double*)malloc(m_dataSize * m_outputDimensions * sizeof(double));
    // Initialize solution (randomly)
    for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
        Y[i] = gaussNumber() * .0001;

    // Perform main training loop
    std::cout << "Input similarities computed. Learning embedding..." << std::endl;

    auto gradients = std::vector<double>(m_dataSize * m_outputDimensions, 0.0);
    auto uY = std::vector<double>(m_dataSize * m_outputDimensions, 0.0);
    auto gains = std::vector<double>(m_dataSize * m_outputDimensions, 1.0);
    for (int iter = 0; iter < m_iterations; iter++) 
    {
        // Compute exact gradient
        computeExactGradient(P, Y, gradients);

        // Update gains
        for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
            gains[i] = (sign(gradients[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
        for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
            gains[i] = std::max(gains[i], 0.1);

        // Perform gradient update (with momentum and gains)
        for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
            uY[i] = momentum * uY[i] - eta * gains[i] * gradients[i];
        for (int i = 0; i < m_dataSize * m_outputDimensions; i++)
            Y[i] = Y[i] + uY[i];

        // Make solution zero-mean
        zeroMean(Y, m_dataSize, m_outputDimensions);

        // Stop lying about the P-values after a while, and switch momentum
        if (iter == stop_lying_iter) 
        {
            for (int i = 0; i < m_dataSize * m_dataSize; i++)
                P[i] /= 12.0;
        }

        if (iter == mom_switch_iter) momentum = final_momentum;

        // Print out progress
        if (iter > 0 && (iter % 50 == 0 || iter == m_iterations - 1)) {
            double C = evaluateError(P, Y);
            std::cout << "Iteration " << (iter + 1) << ": error is " << C << std::endl;
        }
    }

    // Clean up memory
    free(P);

    size_t offset = 0;
    for (size_t i = 0; i < m_dataSize; ++i) {
        auto point = std::vector<double>();
        for (size_t j = 0; j < m_outputDimensions; ++j) {
            point.push_back(Y[offset++]);
        }
        m_result.push_back(point);
    }
    free(Y);
}



//save methods--------------------------------------------------------------------------------------

void TSNE::saveToStream(std::ostream & stream)
{
	for (size_t i = 0; i < m_dataSize; ++i)
	{
		for (size_t j = 0; j < m_outputDimensions; ++j)
		{
			stream << m_result[i][j];
			if (j < m_outputDimensions - 1)
			{
				stream << ',';
			}
		}
		stream << '\n';
	}

    stream.flush();
}

void TSNE::saveToCout()
{
    saveToStream(std::cout);
}

void TSNE::saveCSV()
{
    auto csv_fstream = std::ofstream(m_outputFile + ".csv");
	if (!csv_fstream.is_open())
	{
		std::cerr << "can't open " << m_outputFile << ".csv" << std::endl;
        return;
	}

    saveToStream(csv_fstream);
}

void TSNE::saveLegacy()
{
	auto f = std::ofstream(m_outputFile + ".dat", std::ios::binary);
	if (!f.is_open()) {
		std::cerr << "can't open " << m_outputFile << ".dat" << std::endl;
		return;
	}

	auto landmarks = std::vector<int>(m_dataSize);
	for (size_t i = 0; i < m_dataSize; ++i)
	{
		landmarks[i] = i;
	}

	auto costs = std::vector<double>(m_dataSize, 0.0);

	f.write(reinterpret_cast<char*>(&m_dataSize), sizeof(m_dataSize));
	f.write(reinterpret_cast<char*>(&m_outputDimensions), sizeof(m_outputDimensions));
	for (auto& point : m_result) {
		f.write(reinterpret_cast<char*>(point.data()), m_outputDimensions * sizeof(double));
	}
	f.write(reinterpret_cast<char*>(landmarks.data()), landmarks.size() * sizeof(int));
	f.write(reinterpret_cast<char*>(costs.data()), costs.size() * sizeof(double));

	std::cout << "Wrote the " << m_dataSize << " x " << m_outputDimensions
        << " data matrix successfully!" << std::endl;
}

void TSNE::saveSVG()
{
    double extreme = 0;
    for (size_t i = 0; i < m_dataSize; ++i)
	{
        for (size_t j = 0; j < m_outputDimensions; ++j)
		{
            extreme = std::max(extreme, std::abs(m_result[i][j]));
        }
    }
    double radius = 0.5;
    double halfWidth = extreme + radius;

    //TODO: allow setting a labelFile, e.g. by command line option
    auto labelFile = std::string();
    auto labels = std::vector<uint8_t>();
    bool usingLabels = false;
	if (!labelFile.empty())
	{
		usingLabels = true;
		auto labelInput = std::ifstream(labelFile, std::ios::in | std::ios::binary);
		if (!labelInput.is_open())
		{
            std::cerr << "Could not open " << labelFile << std::endl;
			return;
		}

		uint32_t labelCount;
		labelInput.read(reinterpret_cast<char*>(&labelCount), sizeof(labelCount));
		std::cout << "Labels file contains " << labelCount << " labels." << std::endl;
		if (labelCount < m_dataSize)
		{
			std::cerr << "Not enough labels for result" << std::endl;
			return;
		}

		labelCount = std::min(labelCount, m_dataSize);
		labels.resize(labelCount);
		labelInput.read(reinterpret_cast<char*>(labels.data()), labels.size());
	}

	uint8_t maxLabel = 0;
	for (auto label : labels)
		maxLabel = std::max(label, maxLabel);
	auto colors = std::vector<std::string>();
	for (int i = 0; i <= maxLabel; ++i)
		colors.push_back("hsl(" + std::to_string(360.0 * i / (maxLabel / 2 + 1)) + ", 100%, " + (i % 2 == 0 ? "25" : "60") + "%)");

	auto f = std::ofstream(m_outputFile + ".svg", std::ios::out | std::ios::trunc);
	if (!f.is_open()) {
		std::cerr << "can't open " << m_outputFile << ".svg" << std::endl;
		return;
	}

	f << "<?xml version='1.0' encoding='UTF-8' ?>\n";
	f << "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='600' height='600' viewBox='" << -halfWidth << " " << -halfWidth << " " << 2 * halfWidth << " " << 2 * halfWidth << "'>\n";

	auto color = std::string("black");
	for (unsigned int i = 0; i < m_dataSize; i++)
	{
		if (usingLabels)
        {
            color = labels[i] < colors.size() ? colors[labels[i]] : "black";
        }

		f << "<circle "
			<< "cx='" << m_result[i][0] << "' "
			<< "cy='" << m_result[i][1] << "' "
			<< "fill='" << color << "' "
			<< "r='" << radius << "' "
			<< "stroke='none' opacity='0.5'/>\n";
	}
	f << "</svg>\n";
}

void TSNE::zeroMean(std::vector<std::vector<double>> & data, unsigned int dimensions)
{
	// Compute data mean
	auto mean = std::vector<double>(dimensions, 0.0);
	int nD = 0;
	for (int n = 0; n < m_dataSize; n++) {
		for (int d = 0; d < dimensions; d++) {
			mean[d] += data[n][d];//X[nD + d];
		}
		nD += dimensions;
	}
	for (int d = 0; d < dimensions; d++) {
		mean[d] /= (double)m_dataSize;
	}

	// Subtract data mean
	nD = 0;
	for (int n = 0; n < m_dataSize; n++) {
		for (int d = 0; d < dimensions; d++) {
			data[n][d] -= mean[d];
		}
		nD += dimensions;
	}
}

void TSNE::computeGaussianPerplexity(double* P)
{
	// Compute the squared Euclidean distance matrix
	auto distances = computeSquaredEuclideanDistance(m_data);
    auto DD = distances.data();
	// Compute the Gaussian kernel row by row
	int nN = 0;
	for (int n = 0; n < m_dataSize; n++) {

		// Initialize some variables
		bool found = false;
		double beta = 1.0;
		double min_beta = -DBL_MAX;
		double max_beta = DBL_MAX;
		double tol = 1e-5;
		double sum_P;

		// Iterate until we found a good perplexity
		int iter = 0;
		while (!found && iter < 200) {

			// Compute Gaussian kernel row
			for (int m = 0; m < m_dataSize; m++)
				P[nN + m] = exp(-beta * DD[nN + m]);
			P[nN + n] = DBL_MIN;

			// Compute entropy of current row
			sum_P = DBL_MIN;
			for (int m = 0; m < m_dataSize; m++)
				sum_P += P[nN + m];
			double H = 0.0;
			for (int m = 0; m < m_dataSize; m++)
				H += beta * (DD[nN + m] * P[nN + m]);
			H = (H / sum_P) + log(sum_P);

			// Evaluate whether the entropy is within the tolerance level
			double Hdiff = H - log(m_perplexity);
			if (Hdiff < tol && -Hdiff < tol) {
				found = true;
			}
			else {
				if (Hdiff > 0) {
					min_beta = beta;
					if (max_beta == DBL_MAX || max_beta == -DBL_MAX)
						beta *= 2.0;
					else
						beta = (beta + max_beta) / 2.0;
				}
				else {
					max_beta = beta;
					if (min_beta == -DBL_MAX || min_beta == DBL_MAX)
						beta /= 2.0;
					else
						beta = (beta + min_beta) / 2.0;
				}
			}

			// Update iteration counter
			iter++;
		}

		// Row normalize P
		for (int m = 0; m < m_dataSize; m++)
			P[nN + m] /= sum_P;
		nN += m_dataSize;
	}
}

// compute squared eucl. dist. on OUTPUT!!! maptrix
std::vector<double> TSNE::computeSquaredEuclideanDistance(double* X)
{
    auto distances = std::vector<double>(m_dataSize * m_dataSize, 0.0);
    auto D = m_outputDimensions; // this is for output matrix; must change for input matrix
    const double* XnD = X;
    for (int n = 0; n < m_dataSize; ++n, XnD += D) {
        const double* XmD = XnD + D;
        double* curr_elem = &distances[n*m_dataSize + n];
        double* curr_elem_sym = curr_elem + m_dataSize;
        for (int m = n + 1; m < m_dataSize; ++m, XmD += D, curr_elem_sym += m_dataSize) {
            ++curr_elem;
            for (int d = 0; d < D; ++d) {
                *curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]);
            }
            *curr_elem_sym = *curr_elem;
        }
    }
    return distances;
}

// Compute squared Euclidean distance matrix for INPUT!!! data
std::vector<double> TSNE::computeSquaredEuclideanDistance(std::vector<std::vector<double>> data)
{
	// hacky conversion TODO remove
    assert(data.size() == m_dataSize);
	double* XnD = (double*)malloc(data.size() * data[0].size() * sizeof(double));// = data;
	int offset = 0;
	for (auto & point : data)
	{
        assert(point.size() == m_inputDimensions);
		for (auto & val : point)
		{
			XnD[offset++] = val;
		}
	}

    auto distances = std::vector<double>(m_dataSize * m_dataSize, 0.0);
    auto D = m_inputDimensions;
	for (int n = 0; n < m_dataSize; ++n, XnD += D) {
		const double* XmD = XnD + D;
		double* curr_elem = &distances[n*m_dataSize + n];
		double* curr_elem_sym = curr_elem + m_dataSize;
		for (int m = n + 1; m < m_dataSize; ++m, XmD += D, curr_elem_sym += m_dataSize) {
			++curr_elem;
			for (int d = 0; d < D; ++d) {
				*curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]);
			}
			*curr_elem_sym = *curr_elem;
		}
	}
    return distances;
}

void TSNE::computeGaussianPerplexity(unsigned int** _row_P,	unsigned int** _col_P, std::vector<double> & _val_P) 
{
    assert(m_data.size() == m_dataSize);
    assert(m_data[0].size() == m_inputDimensions);

    //hacky conversion TODO remove
	double* X = (double*)malloc(m_dataSize * m_inputDimensions * sizeof(double));// = data;
	int offset = 0;
	for (auto & point : m_data)
	{
		for (auto & val : point)
		{
			X[offset++] = val;
		}
	}

	int K = (int)(3 * m_perplexity);

	// Allocate the memory we need
	*_row_P = (unsigned int*)malloc((m_dataSize + 1) * sizeof(unsigned int));
	*_col_P = (unsigned int*)calloc(m_dataSize * K, sizeof(unsigned int));
    _val_P.resize(m_dataSize * K, 0.0);
    auto val_P = _val_P.data();
	if (*_row_P == nullptr || *_col_P == nullptr) {
		printf("Memory allocation failed!\n");
		exit(1);
	}
	unsigned int* row_P = *_row_P;
	unsigned int* col_P = *_col_P;
	
	auto cur_P = std::vector<double>(m_dataSize - 1);
	row_P[0] = 0;
	for (int n = 0; n < m_dataSize; n++)
		row_P[n + 1] = row_P[n] + (unsigned int)K;

	// Build ball tree on data set
	auto tree = VpTree<DataPoint, euclidean_distance>();
	auto obj_X = std::vector<DataPoint>(m_dataSize, DataPoint(m_inputDimensions, -1, X));
	for (int n = 0; n < m_dataSize; n++)
		obj_X[n] = DataPoint(m_inputDimensions, n, X + n * m_inputDimensions);
	tree.create(obj_X);

	// Loop over all points to find nearest neighbors
	std::cout << "Building tree..." << std::endl;
	std::vector<DataPoint> indices;
	std::vector<double> distances;
	for (int n = 0; n < m_dataSize; n++) {

		if (n % 10000 == 0)
        {
            std::cout << " - point " << n << " of " << m_dataSize << std::endl;
        }

		// Find nearest neighbors
		indices.clear();
		distances.clear();
		tree.search(obj_X[n], K + 1, &indices, &distances);

		// Initialize some variables for binary search
		bool found = false;
		double beta = 1.0;
		double min_beta = -DBL_MAX;
		double max_beta = DBL_MAX;
		double tol = 1e-5;

		// Iterate until we found a good perplexity
		int iter = 0; double sum_P;
		while (!found && iter < 200) {

			// Compute Gaussian kernel row
			for (int m = 0; m < K; m++)
				cur_P[m] = exp(-beta * distances[m + 1] * distances[m + 1]);

			// Compute entropy of current row
			sum_P = DBL_MIN;
			for (int m = 0; m < K; m++)
				sum_P += cur_P[m];
			double H = .0;
			for (int m = 0; m < K; m++)
				H += beta * (distances[m + 1] * distances[m + 1] * cur_P[m]);
			H = (H / sum_P) + log(sum_P);

			// Evaluate whether the entropy is within the tolerance level
			double Hdiff = H - log(m_perplexity);
			if (Hdiff < tol && -Hdiff < tol) {
				found = true;
			}
			else {
				if (Hdiff > 0) {
					min_beta = beta;
					if (max_beta == DBL_MAX || max_beta == -DBL_MAX)
						beta *= 2.0;
					else
						beta = (beta + max_beta) / 2.0;
				}
				else {
					max_beta = beta;
					if (min_beta == -DBL_MAX || min_beta == DBL_MAX)
						beta /= 2.0;
					else
						beta = (beta + min_beta) / 2.0;
				}
			}

			// Update iteration counter
			iter++;
		}

		// Row-normalize current row of P and store in matrix
		for (unsigned int m = 0; m < K; m++)
			cur_P[m] /= sum_P;
		for (unsigned int m = 0; m < K; m++) {
			col_P[row_P[n] + m] = (unsigned int)indices[m + 1].index();
			val_P[row_P[n] + m] = cur_P[m];
		}
	}
}
