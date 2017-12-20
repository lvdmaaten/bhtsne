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
    , m_gen(static_cast<unsigned long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm) (approximately)
void TSNE::computeGradient(unsigned int *inp_row_P, unsigned int *inp_col_P, double *inp_val_P, double *Y, double *dC)
{

    // Construct space-partitioning tree on current map
    auto tree = SPTree(m_outputDimensions, Y, m_dataSize);

    // Compute all terms required for t-SNE gradient
    auto pos_f = std::vector<double>(m_dataSize * m_outputDimensions, 0.0);
    tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, m_dataSize, pos_f.data());

    double sum_Q = 0.0;
    auto neg_f = std::vector<double>(m_dataSize * m_outputDimensions, 0.0);
    for(unsigned n = 0; n < m_dataSize; n++)
    {
        tree.computeNonEdgeForces(n, m_gradientAccuracy, neg_f.data() + n * m_outputDimensions, &sum_Q);
    }

    // Compute final t-SNE gradient
    for(unsigned i = 0; i < m_dataSize * m_outputDimensions; i++)
    {
        dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
    }
}

// Compute gradient of the t-SNE cost function (exact)
void TSNE::computeExactGradient(double* P, double* Y /*m_result*/, std::vector<double> & gradients)
{
    assert(gradients.size() == m_dataSize * m_outputDimensions);

	// Make sure the current gradient contains zeros
    std::fill(gradients.begin(), gradients.end(), 0.0);

    // Compute the squared Euclidean distance matrix
    auto distances = computeSquaredEuclideanDistance(Y);
    assert(distances.size() == m_dataSize * m_dataSize);

    // Compute Q-matrix and normalization sum
    // https://en.wikipedia.org/wiki/Q-matrix ??
    auto Q = std::vector<double>(m_dataSize * m_dataSize);
    double sum_Q = .0;
    for(int n = 0; n < m_dataSize; n++)
    {
    	for(int m = 0; m < m_dataSize; m++)
        {
            if(n != m)
            {
                auto matrixIndex = n*m_dataSize + m;
                Q[matrixIndex] = 1.0 / (1.0 + distances[matrixIndex]);
                sum_Q += Q[matrixIndex];
            }
        }
    }

	// Perform the computation of the gradient
	for(int n = 0; n < m_dataSize; n++)
    {
    	for(int m = 0; m < m_dataSize; m++)
        {
            if(n != m)
            {
                auto matrixIndex = n*m_dataSize + m;
                double mult = (P[matrixIndex] - (Q[matrixIndex] / sum_Q)) * Q[matrixIndex];
                for(unsigned d = 0; d < m_outputDimensions; d++)
                {
                    gradients[n*m_outputDimensions + d] +=
                            (Y[n*m_outputDimensions + d] - Y[m*m_outputDimensions + d]) * mult;
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

    // Compute Q-matrix and normalization sum
    double sum_Q = std::numeric_limits<double>::min();
    for(int n = 0; n < m_dataSize; n++)
    {
    	for(int m = 0; m < m_dataSize; m++)
        {
            auto matrixIndex = n*m_dataSize + m;
            if(n != m)
            {
                Q[matrixIndex] = 1.0 / (1.0 + distances[matrixIndex]);
                sum_Q += Q[matrixIndex];
            }
            else
            {
                Q[matrixIndex] = std::numeric_limits<double>::min();
            }
        }
    }
    for(auto& each : Q)
    {
        each /= sum_Q;
    }

    // Sum t-SNE error
    double error = 0.0;
	for(int n = 0; n < m_dataSize * m_dataSize; n++)
    {
        error += P[n] * log((P[n] + std::numeric_limits<float>::min())
                            / (Q[n] + std::numeric_limits<float>::min()));
	}

    //TODO: tests still succeed even if any static value is returned...
	return error;
}

// Evaluate t-SNE cost function (approximately)
double TSNE::evaluateError(unsigned int * row_P, unsigned int * col_P, double * val_P, double * Y)
{
    // Get estimate of normalization term
    auto tree = SPTree(m_outputDimensions, Y, m_dataSize);
    auto buff = std::vector<double>(m_outputDimensions, 0.0);
    double sum_Q = 0.0;
    for(unsigned i = 0; i < m_dataSize; i++)
    {
        tree.computeNonEdgeForces(i, m_gradientAccuracy, buff.data(), &sum_Q);
    }

    // Loop over all edges to compute t-SNE error
    double error = 0.0;
    for(unsigned n = 0; n < m_dataSize; n++)
    {
        int ind1 = n * m_outputDimensions;
        for(unsigned i = row_P[n]; i < row_P[n + 1]; i++)
        {
            double Q = 0.0;
            int ind2 = col_P[i] * m_outputDimensions;
            for(unsigned d = 0; d < m_outputDimensions; d++)
            {
                buff[d] = Y[ind1 + d] - Y[ind2 + d];
                Q += buff[d] * buff[d];
            }

            Q = (1.0 / (1.0 + Q)) / sum_Q;
            error += val_P[i] * log((val_P[i] + std::numeric_limits<float>::min())
                                    / (Q + std::numeric_limits<float>::min()));
        }
    }

    //TODO: tests still succeed even if any static value is returned...
    return error;
}

// Symmetrizes a sparse matrix
void TSNE::symmetrizeMatrix(std::vector<unsigned int> & row_P, std::vector<unsigned int> & col_P, std::vector<double> & val_P)
{
    // Count number of elements and row counts of symmetric matrix
    auto row_counts = std::vector<unsigned int>(m_dataSize, 0);
    for(unsigned n = 0; n < m_dataSize; n++)
    {
        for(unsigned i = row_P[n]; i < row_P[n + 1]; i++)
        {
            // Check whether element (col_P[i], n) is present
            auto first = row_P.begin() + row_P[col_P[i]];
            auto last = row_P.begin() + row_P[col_P[i] + 1];

            if(std::find(first, last, n) == last)
            {
                row_counts[col_P[i]]++;
            }
            row_counts[n]++;
        }
    }
    unsigned no_elem = std::accumulate(row_counts.begin(), row_counts.begin() + m_dataSize, 0u);

    // Allocate memory for symmetrized matrix
    // TODO reuse the memory in row,col and val!!
    auto sym_row_P = std::vector<unsigned int>(m_dataSize + 1);
    auto sym_col_P = std::vector<unsigned int>(no_elem);
    auto sym_val_P = std::vector<double>(no_elem);

    // Construct new row indices for symmetric matrix
    sym_row_P[0] = 0;
    for(int n = 0; n < m_dataSize; n++)
    {
        sym_row_P[n + 1] = sym_row_P[n] + row_counts[n];
    }

    // Fill the result matrix
    auto offset = std::vector<int>(m_dataSize, 0);
    for(unsigned n = 0; n < m_dataSize; n++)
    {
        for(unsigned i = row_P[n]; i < row_P[n + 1]; i++)
        { // considering element(n, col_P[i])

            // Check whether element (col_P[i], n) is present
            bool present = false;
            for(unsigned m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++)
            {
                if(col_P[m] == n)
                {
                    present = true;
                    if(n <= col_P[i])
                    { // make sure we do not add elements twice
                        sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                        sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                        sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
                        sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
                    }
                }
            }

            // If (col_P[i], n) is not present, there is no addition involved
            if(!present)
            {
                sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
                sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
            }

            // Update offsets
            if(!present || n <= col_P[i])
            {
                offset[n]++;
                if(col_P[i] != n)
                {
                    offset[col_P[i]]++;
                }
            }
        }
    }

    // Divide the result by two
    for(int i = 0; i < no_elem; i++)
    {
        sym_val_P[i] /= 2.0;
    }

    // Return symmetrized matrices
    row_P = std::move(sym_row_P);
    col_P = std::move(sym_col_P);
    val_P = std::move(sym_val_P);
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

void TSNE::setRandomSeed(unsigned long seed)
{
    std::cout << "Using random seed: " << seed << std::endl;
    m_gen.seed(seed);
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
	for (size_t i = 0; i < m_dataSize; ++i)
    {
		auto point = std::vector<double>(m_inputDimensions);
		f.read(reinterpret_cast<char*>(point.data()), sizeof(double) * m_inputDimensions);
        m_data.emplace_back(std::move(point));
	}

    //read seed
	if (!f.eof())
    {
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
    if (m_dataSize - 1 < 3 * m_perplexity)
    {
        auto message = "perplexity (perplexity=" + std::to_string(m_perplexity) +
            ") has to be smaller than a third of the dataSize (dataSize=" + std::to_string(m_dataSize) + ")";
        std::cerr << message << std::endl;
        throw std::invalid_argument(message);
    }

    std::cout << "Using:"
        << "\ndata size " << m_dataSize
        << "\nin dimensions " << m_inputDimensions
        << "\nout dimensions " << m_outputDimensions
        << "\nperplexity " << m_perplexity
        << "\ngradient accuracy " << m_gradientAccuracy
        << std::endl;

    if (m_gradientAccuracy == 0.0)
    {
        runExact();
    }
    else
    {
        runApproximation();
    }
}


void TSNE::runApproximation()
{
	// Normalize input data to prevent numerical problems
	std::cout << "Computing input similarities..." << std::endl;
    zeroMean(m_data);
	// TODO: extract normalization function for vector
	double max_X = 0.0;
	for (auto & point : m_data)
	{
		for (auto & val : point)
		{
            max_X = std::max(std::fabs(val), max_X);
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
    // init sparse matrix P, TODO: create and use class SparseMatrix
    auto row_P = std::vector<unsigned int>();
    auto col_P = std::vector<unsigned int>();
    auto val_P = std::vector<double>();

	// Compute asymmetric pairwise input similarities
	computeGaussianPerplexity(row_P, col_P, val_P);

	// Symmetrize input similarities
	symmetrizeMatrix(row_P, col_P, val_P);

	//normalize val_P so that sum of all val = 1
	double sum_P = std::accumulate(val_P.begin(), val_P.end(), 0.0);
	for (auto & each : val_P)
    {
        each /= sum_P;
    	// Lie about the P-values
        each *= 12.0;
    }

	auto Y = std::vector<double>(m_dataSize * m_outputDimensions);
	// Initialize solution (randomly)
	for (auto & each : Y)
    {
        each = gaussNumber() * 0.0001;
    }



    //TODO: documentation for all these magic numbers
    int stop_lying_iter = 250;
    int mom_switch_iter = 250;

    // Set learning parameters
    double momentum = 0.5;
    double final_momentum = 0.8;
    double eta = 200.0;

    auto dY = std::vector<double>(m_dataSize * m_outputDimensions);
    auto uY = std::vector<double>(m_dataSize * m_outputDimensions, 0.0);
    auto gains = std::vector<double>(m_dataSize * m_outputDimensions, 1.0);

	// Perform main training loop
    std::cout << " Input similarities computed. Learning embedding..." << std::endl;
	for (unsigned iter = 0; iter < m_iterations; iter++)
    {
		// Compute approximate gradient
        computeGradient(row_P.data(), col_P.data(), val_P.data(), Y.data(), dY.data());

		// Update gains
		for (unsigned i = 0; i < m_dataSize * m_outputDimensions; i++)
        {
            if (sign(dY[i]) != sign(uY[i]))
            {
                gains[i] += 0.2;
            }
            else
            {
                gains[i] *= 0.8;
            }
            gains[i] = std::max(0.1, gains[i]);
        }

		// Perform gradient update (with momentum and gains)
		for (unsigned i = 0; i < m_dataSize * m_outputDimensions; i++)
        {
            uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
            Y[i] += uY[i];
        }

		// Make solution zero-mean
		zeroMean(Y.data(), m_dataSize, m_outputDimensions);

		// Stop lying about the P-values after a while, and switch momentum
		if (iter == stop_lying_iter)
        {
			for (auto & each : val_P)
            {
                each /= 12.0;
            }
		}
		if (iter == mom_switch_iter)
        {
            momentum = final_momentum;
        }

		// Print out progress
		if (iter % 50 == 49 || iter == m_iterations - 1)
        {
			// doing approximate computation here!
			double error = evaluateError(row_P.data(), col_P.data(), val_P.data(), Y.data());
			std::cout << "Iteration " << (iter + 1) << ": error is " << error << std::endl;
		}
	}

    // convert Y to vector<vector> result
	size_t offset = 0;
	for (size_t i = 0; i < m_dataSize; ++i)
    {
		auto point = std::vector<double>();
		for (size_t d = 0; d < m_outputDimensions; ++d)
        {
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
    double momentum = 0.5;
    double final_momentum = 0.8;
    double eta = 200.0;

    // Normalize input data (to prevent numerical problems)
    std::cout << "Computing input similarities..." << std::endl;
    zeroMean(m_data);
    // TODO: extract normalization function for vector
    double max_X = 0.0;
    for (auto & point : m_data)
    {
        for (auto & val : point)
        {
            max_X = std::max(std::fabs(val), max_X);
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
    for (int n = 0; n < m_dataSize; n++)
    {
        int mN = (n + 1) * m_dataSize;
        for (int m = n + 1; m < m_dataSize; m++)
        {
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
        if (iter > 0 && (iter % 50 == 0 || iter == m_iterations - 1))
        {
            double C = evaluateError(P, Y);
            std::cout << "Iteration " << (iter + 1) << ": error is " << C << std::endl;
        }
    }

    // Clean up memory
    free(P);

    size_t offset = 0;
    for (size_t i = 0; i < m_dataSize; ++i)
    {
        auto point = std::vector<double>();
        for (size_t j = 0; j < m_outputDimensions; ++j)
        {
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
	if (!f.is_open())
    {
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
	for (auto& point : m_result)
    {
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
	if (!f.is_open())
    {
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

// TODO: remove, always use vector as parameter
void TSNE::zeroMean(double* X, int N, int D)
{

    auto points = std::vector<std::vector<double>>(m_dataSize, std::vector<double>(D, 0.0));

    for (auto i = 0u; i < N; ++i)
    {
        for (auto j = 0u; j < D; ++j)
        {
            points[i][j] = X[i*D + j];
        }
    }

    zeroMean(points);

    for (auto i = 0u; i < N; ++i)
    {
        for (auto j = 0u; j < D; ++j)
        {
            X[i*D + j] = points[i][j];
        }
    }
}

//make the mean of all data points equal 0 for each dimension -> zero mean
void TSNE::zeroMean(std::vector<std::vector<double>> & points)
{
    auto dimensions = points[0].size();

    for (auto d = 0u; d < dimensions; d++)
    {
        auto mean = 0.0;
        for (auto i = 0u; i < points.size(); i++)
        {
            mean += points[i][d];
        }
        mean /= points.size();
        for (auto i = 0u; i < points.size(); i++)
        {
            points[i][d] -= mean;
        }
    }
}

void TSNE::computeGaussianPerplexity(double* P)
{
	// Compute the squared Euclidean distance matrix
	auto distances = computeSquaredEuclideanDistance(m_data);
    auto DD = distances.data();
	// Compute the Gaussian kernel row by row
	int nN = 0;
	for (int n = 0; n < m_dataSize; n++)
    {

		// Initialize some variables
		bool found = false;
		double beta = 1.0;
		double min_beta = -DBL_MAX;
		double max_beta = DBL_MAX;
		double tol = 1e-5;
		double sum_P;

		// Iterate until we found a good perplexity
		int iter = 0;
		while (!found && iter < 200)
        {

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
			if (Hdiff < tol && -Hdiff < tol)
            {
				found = true;
			}
			else {
				if (Hdiff > 0)
                {
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
//TODO: remove this, always use vector as parameter
std::vector<double> TSNE::computeSquaredEuclideanDistance(double* data)
{
    auto points = std::vector<std::vector<double>>(m_dataSize, std::vector<double>(m_outputDimensions, 0.0));

    for (auto i = 0u; i < points.size(); ++i)
    {
        for (auto j = 0u; j < points[0].size(); ++j)
        {
            points[i][j] = data[i*m_outputDimensions + j];
        }
    }

    return computeSquaredEuclideanDistance(points);
}

std::vector<double> TSNE::computeSquaredEuclideanDistance(const std::vector<std::vector<double>> & points)
{
    auto dimensions = points[0].size();

    auto distances = std::vector<double>(points.size() * points.size(), 0.0);

    for (auto i = 0u; i < points.size(); ++i)
    {
        for (auto j = i + 1; j < points.size(); ++j)
        {
            auto distance = 0.0;
            for (auto d = 0u; d < dimensions; ++d)
            {
                auto diff = points[i][d] - points[j][d];
                distance += diff * diff;
            }

            distances[i*points.size() + j] = distance;
            distances[j*points.size() + i] = distance;
        }
    }

    return distances;
}

void TSNE::computeGaussianPerplexity(std::vector<unsigned int> & row_P, std::vector<unsigned int> & col_P, std::vector<double> & val_P)
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
	row_P.resize(m_dataSize + 1);
    col_P.resize(m_dataSize * K);
    val_P.resize(m_dataSize * K, 0.0);

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
	for (int n = 0; n < m_dataSize; n++)
    {

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
		while (!found && iter < 200)
        {

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
			if (Hdiff < tol && -Hdiff < tol)
            {
				found = true;
			}
			else {
				if (Hdiff > 0)
                {
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
		for (unsigned int m = 0; m < K; m++)
        {
			col_P[row_P[n] + m] = (unsigned int)indices[m + 1].index();
			val_P[row_P[n] + m] = cur_P[m];
		}
	}
}
