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


#include <bhtsne/TSNE.h>


#include <cassert>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>

#include "SpacePartitioningTreeTemplate.h"
#include "VantagePointTree.h"

#include <bhtsne/bhtsne-version.h> // includes BHTSNE_VERSION macro


using namespace bhtsne;


TSNE::TSNE()
    : m_perplexity(50.0)
    , m_gradientAccuracy(0.2)
    , m_iterations(1000)
    , m_outputDimensions(2)
    , m_inputDimensions(0)
    , m_dataSize(0)
    , m_gen(static_cast<unsigned long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
    , m_outputFile("result")
{
}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm) (approximately)
template<unsigned int D>
Vector2D<double> TSNE::computeGradient(SparseMatrix & similarities)
{
    // Construct space-partitioning tree on current map
    auto tree = SpacePartitioningTree<D>(m_result);

    // Compute all terms required for t-SNE gradient
    auto pos_f = Vector2D<double>(m_dataSize, m_outputDimensions, 0.0);
    tree.computeEdgeForces(similarities.rows, similarities.columns, similarities.values, pos_f);

    auto neg_f = Vector2D<double>(m_dataSize, m_outputDimensions, 0.0);
    double sum_Q = 0.0;
    // omp version on windows (2.0) does only support signed loop variables, should be unsigned
    #pragma omp parallel for reduction(+:sum_Q)
    for (int n = 0; n < m_dataSize; ++n)
    {
        tree.computeNonEdgeForces(n, m_gradientAccuracy, neg_f[n], sum_Q);
    }

    auto result = Vector2D<double>(m_dataSize, m_outputDimensions);
    // Compute final t-SNE gradient
    for (unsigned int i = 0; i < m_dataSize; ++i)
    {
        for (unsigned int j = 0; j < m_outputDimensions; ++j)
        {
            result[i][j] = pos_f[i][j] - (neg_f[i][j] / sum_Q);
        }
    }
    return result;
}

// Compute gradient of the t-SNE cost function (exact)
Vector2D<double> TSNE::computeGradientExact(const Vector2D<double> & Perplexity)
{
    auto gradients = Vector2D<double>(m_dataSize, m_outputDimensions, 0.0);

    // Compute the squared Euclidean distance matrix
    auto distances = computeSquaredEuclideanDistance(m_result);
    assert(distances.height() == m_dataSize);
    assert(distances.width() == m_dataSize);

    // Compute Q-matrix and normalization sum
    // Q = similarities of low dimensional output data
    auto Q = Vector2D<double>(m_dataSize, m_dataSize);
    double sum_Q = 0.0;
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        for (unsigned int m = n + 1; m < m_dataSize; ++m)
        {
            Q[n][m] = 1.0 / (1.0 + distances[n][m]);
            Q[m][n] = Q[n][m];
            sum_Q += 2 * Q[n][m];
        }
    }

	// Perform the computation of the gradient
	for (unsigned int n = 0; n < m_dataSize; ++n)
    {
    	for (unsigned int m = 0; m < m_dataSize; ++m)
        {
            if (n == m)
            {
                continue;
            }

            double mult = (Perplexity[n][m] - (Q[n][m] / sum_Q)) * Q[n][m];
            for (unsigned int d = 0; d < m_outputDimensions; ++d)
            {
                gradients[n][d] += (m_result[n][d] - m_result[m][d]) * mult;
            }
		}
	}

    return gradients;
}


// Evaluate t-SNE cost function (exactly)
double TSNE::evaluateErrorExact(const Vector2D<double> & Perplexity)
{
    assert(Perplexity.height() == m_dataSize);
    assert(Perplexity.width() == m_dataSize);

    // Compute the squared Euclidean distance matrix
    auto Q = Vector2D<double>(m_dataSize, m_dataSize, std::numeric_limits<double>::min());
    auto distances = computeSquaredEuclideanDistance(m_result);
    assert(distances.height() == m_result.height());
    assert(distances.width() == m_result.height());

    // Compute Q-matrix and normalization sum
    //TODO init to 0 (or evaluate consequences)
    double sum_Q = std::numeric_limits<double>::min();
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        for (unsigned int m = n + 1; m < m_dataSize; ++m)
        {
            Q[n][m] = 1.0 / (1.0 + distances[n][m]);
            Q[m][n] = Q[n][m];
            sum_Q += 2 * Q[n][m];
        }
    }

    //TODO use vector normalization method
    for (auto & each : Q)
    {
        each /= sum_Q;
    }

    // Sum t-SNE error
    //TODO remove numeric limits, coz Q is init like this
    double error = 0.0;
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        for (unsigned int m = 0; m < m_dataSize; ++m)
        {
            error += Perplexity[n][m] * log((Perplexity[n][m] + std::numeric_limits<float>::min())
                / (Q[n][m] + std::numeric_limits<float>::min()));
        }
    }

	return error;
}

// Evaluate t-SNE cost function (approximately)
template<unsigned int D>
double TSNE::evaluateError(SparseMatrix & similarities)
{
    // Get estimate of normalization term
    auto tree = SpacePartitioningTree<D>(m_result);
    auto buff = std::vector<double>(m_outputDimensions, 0.0);
    double sum_Q = 0.0;
    for (unsigned int i = 0; i < m_dataSize; ++i)
    {
        tree.computeNonEdgeForces(i, m_gradientAccuracy, buff.data(), sum_Q);
    }

    // Loop over all edges to compute t-SNE error
    double error = 0.0;
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        for (unsigned int i = similarities.rows[n]; i < similarities.rows[n + 1]; ++i)
        {
            double Q = 0.0;
            for (unsigned int d = 0; d < m_outputDimensions; d++)
            {
                buff[d] = m_result[n][d] - m_result[similarities.columns[i]][d];
                Q += buff[d] * buff[d];
            }

            Q = (1.0 / (1.0 + Q)) / sum_Q;
            error += similarities.values[i] * log((similarities.values[i] + std::numeric_limits<float>::min())
                                                  / (Q + std::numeric_limits<float>::min()));
        }
    }

    //TODO: write tests for evaluate error
    return error;
}

// Symmetrizes a sparse matrix
void TSNE::symmetrizeMatrix(SparseMatrix & similarities)
{
    // Count number of elements and row counts of symmetric matrix
    auto row_counts = std::vector<unsigned int>(m_dataSize, 0);
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        for (unsigned int i = similarities.rows[n]; i < similarities.rows[n + 1]; ++i)
        {
            // Check whether element (similarities.columns[i], n) is present
            auto first = similarities.columns.begin() + similarities.rows[similarities.columns[i]];
            auto last = similarities.columns.begin() + similarities.rows[similarities.columns[i] + 1];

            if (std::find(first, last, n) == last)
            {
                row_counts[similarities.columns[i]]++;
            }
            row_counts[n]++;
        }
    }
    unsigned int numberOfElements = std::accumulate(row_counts.begin(), row_counts.begin() + m_dataSize, 0u);

    // Allocate memory for symmetrized matrix
    // TODO reuse the memory in similarities
    auto symmetrized_similarities = SparseMatrix();
    symmetrized_similarities.rows.resize(m_dataSize + 1);
    symmetrized_similarities.columns.resize(numberOfElements);
    symmetrized_similarities.values.resize(numberOfElements);

    // Construct new row indices for symmetric matrix
    symmetrized_similarities.rows[0] = 0;
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        symmetrized_similarities.rows[n + 1] = symmetrized_similarities.rows[n] + row_counts[n];
    }

    // Fill the result matrix
    auto offset = std::vector<int>(m_dataSize, 0);
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        for (unsigned int i = similarities.rows[n]; i < similarities.rows[n + 1]; ++i)
        { // considering element(n, similarities.columns[i])

            // Check whether element (similarities.columns[i], n) is present
            bool present = false;
            for (unsigned int m = similarities.rows[similarities.columns[i]]; m < similarities.rows[similarities.columns[i] + 1]; ++m)
            {
                if (similarities.columns[m] == n)
                {
                    present = true;
                    if (n <= similarities.columns[i])
                    { // make sure we do not add elements twice
                        symmetrized_similarities.columns[symmetrized_similarities.rows[n] + offset[n]] = similarities.columns[i];
                        symmetrized_similarities.columns[symmetrized_similarities.rows[similarities.columns[i]] + offset[similarities.columns[i]]] = n;
                        symmetrized_similarities.values[symmetrized_similarities.rows[n] + offset[n]] = similarities.values[i] + similarities.values[m];
                        symmetrized_similarities.values[symmetrized_similarities.rows[similarities.columns[i]] + offset[similarities.columns[i]]] =
                                similarities.values[i] + similarities.values[m];
                    }
                }
            }

            // If (similarities.columns[i], n) is not present, there is no addition involved
            if (!present)
            {
                symmetrized_similarities.columns[symmetrized_similarities.rows[n] + offset[n]] = similarities.columns[i];
                symmetrized_similarities.columns[symmetrized_similarities.rows[similarities.columns[i]] + offset[similarities.columns[i]]] = n;
                symmetrized_similarities.values[symmetrized_similarities.rows[n] + offset[n]] = similarities.values[i];
                symmetrized_similarities.values[symmetrized_similarities.rows[similarities.columns[i]] + offset[similarities.columns[i]]] = similarities.values[i];
            }

            // Update offsets
            if (!present || n <= similarities.columns[i])
            {
                offset[n]++;
                if (similarities.columns[i] != n)
                {
                    offset[similarities.columns[i]]++;
                }
            }
        }
    }

    // Divide the result by two
    for (auto & each : symmetrized_similarities.values)
    {
        each /= 2.0;
    }

    // Return symmetrized matrices
    similarities.rows = std::move(symmetrized_similarities.rows);
    similarities.columns = std::move(symmetrized_similarities.columns);
    similarities.values = std::move(symmetrized_similarities.values);
}

// with mean zero and standard deviation one
double bhtsne::TSNE::gaussNumber()
{
    // Knuth, Art of Computer Programming vol v2, Section 3.4.1, Algorithm P (p.117)

    double S, V1, V2;
    // executed 1.27 times on average
    do
    {
        // V1, V2 uniformly distributed between -1 and +l.
        V1 = 2.0 * static_cast<double>(m_gen()) / m_gen.max() - 1.0;
        V2 = 2.0 * static_cast<double>(m_gen()) / m_gen.max() - 1.0;
        S = V1 * V1 + V2 * V2;
    }
    while (S >= 1);

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

std::string TSNE::outputFile() const
{
	return m_outputFile;
}

void TSNE::setOutputFile(const std::string & file)
{
	m_outputFile = file;
}


//load methods--------------------------------------------------------------------------------------

bool bhtsne::TSNE::loadFromStream(std::istream & stream)
{
    char separator = ',';

    //read data points
    auto line = std::string();
    bool first = true;
    while (std::getline(stream, line))
    {
        std::istringstream iss(line);
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
            m_inputDimensions = static_cast<unsigned int>(point.size());
            m_data.initialize(0, point.size());
        }

        //fail if inconsistent dimensionality
        if (m_inputDimensions != point.size() || m_inputDimensions == 0)
        {
            m_inputDimensions = 0;
            m_data.initialize(0, 0);
            return false;
        }

        m_data.appendRow(point);
    }

    if (m_data.height() == 0 || m_data.width() == 0)
    {
        return false;
    }

    m_dataSize = static_cast<unsigned int>(m_data.height());

    return true;
}

bool TSNE::loadLegacy(const std::string & file)
{
    std::ifstream f(file, std::ios::binary);
	if (!f.is_open())
    {
        std::cerr << "Could not open " << file << std::endl;
        return false;
    }

    //read params
	f.read(reinterpret_cast<char *>(&m_dataSize), sizeof(m_dataSize));
	f.read(reinterpret_cast<char *>(&m_inputDimensions), sizeof(m_inputDimensions));
	f.read(reinterpret_cast<char *>(&m_gradientAccuracy), sizeof(m_gradientAccuracy));
	f.read(reinterpret_cast<char *>(&m_perplexity), sizeof(m_perplexity));
	f.read(reinterpret_cast<char *>(&m_outputDimensions), sizeof(m_outputDimensions));
	f.read(reinterpret_cast<char *>(&m_iterations), sizeof(m_iterations));

    //read data
    m_data.initialize(m_dataSize, m_inputDimensions);
	f.read(reinterpret_cast<char *>(m_data[0]), m_dataSize * sizeof(double) * m_inputDimensions);

    //read seed
	if (!f.eof())
    {
        int seed;
		f.read(reinterpret_cast<char *>(&seed), sizeof(seed));
        setRandomSeed(static_cast<unsigned long>(seed));
	}

	return true;
}

bool TSNE::loadCSV(const std::string & file)
{
    std::ifstream f(file);
	if (!f.is_open())
    {
        std::cerr << "Could not open " << file << std::endl;
        return false;
    }

    return loadFromStream(f);
}

bool TSNE::loadTSNE(const std::string & file)
{
    std::ifstream f(file, std::ios::binary);
	if (!f.is_open())
    {
        std::cerr << "Could not open " << file << std::endl;
        return false;
    }

	f.read(reinterpret_cast<char *>(&m_dataSize), sizeof(m_dataSize));
	f.read(reinterpret_cast<char *>(&m_inputDimensions), sizeof(m_inputDimensions));

    //read data
    m_data.initialize(m_dataSize, m_inputDimensions);
    f.read(reinterpret_cast<char *>(m_data[0]), m_dataSize * sizeof(double) * m_inputDimensions);

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

    m_result.initialize(m_dataSize, m_outputDimensions);
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
    normalize(m_data);

    // Compute input similarities for exact t-SNE
    auto inputSimilarities = SparseMatrix();

	// Compute asymmetric pairwise input similarities
	computeGaussianPerplexity(inputSimilarities);

	// Symmetrize input similarities
	symmetrizeMatrix(inputSimilarities);

	//normalize inputSimilarities so that sum of all values = 1
	double sum_P = std::accumulate(inputSimilarities.values.begin(), inputSimilarities.values.end(), 0.0);
	for (auto & each : inputSimilarities.values)
    {
        each /= sum_P;
    	// Lie about the inputSimilarities
        each *= 12.0;
    }

	// Initialize solution (randomly)
    for (auto & each : m_result)
    {
        each = gaussNumber() * 0.0001;
    }

    //TODO: documentation for all these magic numbers
    unsigned int stop_lying_iteration = 251;
    unsigned int momentum_switch_iteration = 251;

    // Set learning parameters
    double momentum = 0.5;
    double final_momentum = 0.8;
    double eta = 200.0;

    auto uY = Vector2D<double>(m_dataSize, m_outputDimensions);
    auto gains = Vector2D<double>(m_dataSize, m_outputDimensions, 1.0);

	// Perform main training loop
    std::cout << " Input similarities computed. Learning embedding..." << std::endl;
	for (unsigned int iteration = 1; iteration <= m_iterations; ++iteration)
    {
		// Compute approximate gradient
        auto gradients =
            (m_outputDimensions == 2) ? computeGradient<2>(inputSimilarities) :
            (m_outputDimensions == 3) ? computeGradient<3>(inputSimilarities) :
            computeGradient<0>(inputSimilarities);

		// Update gains
        for (unsigned int i = 0; i < m_dataSize; ++i)
        {
            for (unsigned int j = 0; j < m_outputDimensions; ++j)
            {
                if (sign(gradients[i][j]) != sign(uY[i][j]))
                {
                    gains[i][j] += 0.2;
                }
                else
                {
                    gains[i][j] *= 0.8;
                }
                gains[i][j] = std::max(0.1, gains[i][j]);
            }
        }

		// Perform gradient update (with momentum and gains)
        for (unsigned int i = 0; i < m_dataSize; ++i)
        {
            for (unsigned int j = 0; j < m_outputDimensions; ++j)
            {
                uY[i][j] = momentum * uY[i][j] - eta * gains[i][j] * gradients[i][j];
                m_result[i][j] += uY[i][j];
            }
        }

		// Make solution zero-mean
		zeroMean(m_result);

		// Stop lying about the inputSimilarities-values after a while, and switch momentum
		if (iteration == stop_lying_iteration)
        {
			for (auto & each : inputSimilarities.values)
            {
                each /= 12.0;
            }
		}

		if (iteration == momentum_switch_iteration)
        {
            momentum = final_momentum;
        }

		// Print out progress
		if (iteration % 50 == 0 || iteration == m_iterations)
        {
			// doing approximate computation here!
			double error =
                (m_outputDimensions == 2) ? evaluateError<2>(inputSimilarities) :
                (m_outputDimensions == 3) ? evaluateError<3>(inputSimilarities) :
                evaluateError<0>(inputSimilarities); // assert(false)
			std::cout << "Iteration " << iteration << ": error is " << error << std::endl;
		}
	}
}


void TSNE::runExact()
{
    unsigned int stop_lying_iteration = 251;
    unsigned int momentum_switch_iteration = 251;

    // Set learning parameters
    double momentum = 0.5;
    double final_momentum = 0.8;
    double eta = 200.0;

    // Normalize input data (to prevent numerical problems)
    std::cout << "Computing input similarities..." << std::endl;
    zeroMean(m_data);
    normalize(m_data);

    // Compute input similarities for exact t-SNE
    auto P = computeGaussianPerplexityExact();
    assert(P.width() == m_dataSize);
    assert(P.height() == m_dataSize);

    // Symmetrize input similarities
    std::cout << "Symmetrizing..." << std::endl;
    for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        for (unsigned int m = n + 1; m < m_dataSize; ++m)
        {
            P[n][m] += P[m][n];
            P[m][n] = P[n][m];
        }
    }

    double sum_P = std::accumulate(P.begin(), P.end(), 0.0);
    for (auto & each : P)
    {
        each /= sum_P;
    }

    // Lie about the P-values
    for (auto & each : P)
    {
        each *= 12.0;
    }

    // Initialize solution (randomly)
    for (auto & each : m_result)
    {
        each = gaussNumber() * 0.0001;
    }

    // Perform main training loop
    std::cout << "Input similarities computed. Learning embedding..." << std::endl;

    auto uY    = Vector2D<double>(m_dataSize, m_outputDimensions, 0.0);
    auto gains = Vector2D<double>(m_dataSize, m_outputDimensions, 1.0);

    for (unsigned int iteration = 1; iteration <= m_iterations; ++iteration)
    {
        // Compute exact gradient
        auto gradients = computeGradientExact(P);
        assert(gradients.height() == m_dataSize);
        assert(gradients.width() == m_outputDimensions);

        // Update gains
        for (size_t i = 0; i < m_dataSize; ++i)
        {
            for (size_t j = 0; j < m_outputDimensions; ++j)
            {
                gains[i][j] = (sign(gradients[i][j]) != sign(uY[i][j])) ? (gains[i][j] + .2) : (gains[i][j] * .8);
            }
        }

        for (auto & each : gains)
        {
            each = std::max(each, 0.1);
        }

        // Perform gradient update (with momentum and gains)
        for (size_t i = 0; i < m_dataSize; ++i)
        {
            for (size_t j = 0; j < m_outputDimensions; ++j)
            {
                uY[i][j] = momentum * uY[i][j] - eta * gains[i][j] * gradients[i][j];
            }
        }

        for (size_t i = 0; i < m_dataSize; ++i)
        {
            for (size_t j = 0; j < m_outputDimensions; ++j)
            {
                m_result[i][j] += uY[i][j];
            }
        }

        // Make solution zero-mean
        zeroMean(m_result);

        // Stop lying about the P-values after a while, and switch momentum
        if (iteration == stop_lying_iteration)
        {
            for (auto & each : P)
            {
                each /= 12.0;
            }
        }

        if (iteration == momentum_switch_iteration)
        {
            momentum = final_momentum;
        }

        // Print out progress
        if (iteration % 50 == 0 || iteration == m_iterations)
        {
            double C = evaluateErrorExact(P);
            std::cout << "Iteration " << (iteration + 1) << ": error is " << C << std::endl;
        }
    }
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
    std::ofstream csv_fstream(m_outputFile + ".csv");
	if (!csv_fstream.is_open())
	{
		std::cerr << "can't open " << m_outputFile << ".csv" << std::endl;
        return;
	}

    saveToStream(csv_fstream);
}

void TSNE::saveLegacy()
{
    std::ofstream f(m_outputFile + ".dat", std::ios::binary);
	if (!f.is_open())
    {
		std::cerr << "can't open " << m_outputFile << ".dat" << std::endl;
		return;
	}

	auto landmarks = std::vector<int>(m_dataSize);
    std::iota(landmarks.begin(), landmarks.end(), 0);

	auto costs = std::vector<double>(m_dataSize, 0.0);

	f.write(reinterpret_cast<char *>(&m_dataSize), sizeof(m_dataSize));
	f.write(reinterpret_cast<char *>(&m_outputDimensions), sizeof(m_outputDimensions));
	f.write(reinterpret_cast<char *>(m_result[0]), m_dataSize * m_outputDimensions * sizeof(double));
	f.write(reinterpret_cast<char *>(landmarks.data()), landmarks.size() * sizeof(int));
	f.write(reinterpret_cast<char *>(costs.data()), costs.size() * sizeof(double));

	std::cout << "Wrote the " << m_dataSize << " x " << m_outputDimensions
        << " data matrix successfully!" << std::endl;
}

void TSNE::saveSVG()
{
    assert(m_result.size() > 0);
    double extreme = *std::max_element(m_result.begin(), m_result.end());

    double radius = 0.5;
    double halfWidth = extreme + radius;

    //TODO: allow setting a labelFile, e.g. by command line option
    auto labelFile = std::string();
    auto labels = std::vector<uint8_t>();
    auto usingLabels = false;
	if (!labelFile.empty())
	{
		usingLabels = true;
        std::ifstream labelInput(labelFile, std::ios::in | std::ios::binary);
		if (!labelInput.is_open())
		{
            std::cerr << "Could not open " << labelFile << std::endl;
			return;
		}

		uint32_t labelCount;
		labelInput.read(reinterpret_cast<char *>(&labelCount), sizeof(labelCount));
		std::cout << "Labels file contains " << labelCount << " labels." << std::endl;
		if (labelCount < m_dataSize)
		{
			std::cerr << "Not enough labels for result" << std::endl;
			return;
		}

		labelCount = std::min(labelCount, m_dataSize);
		labels.resize(labelCount);
		labelInput.read(reinterpret_cast<char *>(labels.data()), labels.size());
	}

	uint8_t maxLabel = 0;
	for (uint8_t label : labels)
    {
        maxLabel = std::max(label, maxLabel);
    }
	auto colors = std::vector<std::string>();
	for (uint8_t i = 0; i <= maxLabel; ++i)
    {
        colors.push_back("hsl(" + std::to_string(360.0 * i / (maxLabel / 2 + 1)) + ", 100%, " + (i % 2 == 0 ? "25" : "60") + "%)");
    }

    std::ofstream f(m_outputFile + ".svg", std::ios::out | std::ios::trunc);
	if (!f.is_open())
    {
		std::cerr << "can't open " << m_outputFile << ".svg" << std::endl;
		return;
	}

	f << "<?xml version='1.0' encoding='UTF-8' ?>\n";
	f << "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='600' height='600' viewBox='" << -halfWidth << " " << -halfWidth << " " << 2 * halfWidth << " " << 2 * halfWidth << "'>\n";

	auto color = std::string("black");
	for (unsigned int i = 0; i < m_dataSize; ++i)
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

//make the mean of all data points equal 0 for each dimension -> zero mean
void TSNE::zeroMean(Vector2D<double> & points)
{
    const auto dimensions = points.width();
    const auto size = points.height();

    for (size_t d = 0; d < dimensions; d++)
    {
        auto mean = 0.0;
        for (unsigned int i = 0; i < size; ++i)
        {
            mean += points[i][d];
        }
        mean /= size;
        for (unsigned int i = 0; i < size; ++i)
        {
            points[i][d] -= mean;
        }
    }
}

void TSNE::normalize(Vector2D<double> & vec)
{
    assert(vec.size() > 0);
    double max_X = *std::max_element(vec.begin(), vec.end());
    for (auto & each : vec)
    {
        each /= max_X;
    }
}

Vector2D<double> TSNE::computeGaussianPerplexityExact()
{
	// Compute the squared Euclidean distance matrix
	auto distances = computeSquaredEuclideanDistance(m_data);
    auto P = Vector2D<double>(m_dataSize, m_dataSize);

	// Compute the Gaussian kernel row by row
	for (unsigned int n = 0; n < m_dataSize; ++n)
    {
		// Initialize some variables
		double beta = 1.0;
		double min_beta = std::numeric_limits<double>::lowest();
		double max_beta = std::numeric_limits<double>::max();
		double tolerance_threshold = 1e-5;
        double sum_P = 0.0;

		// Iterate until we found a good perplexity
		for (unsigned int iteration = 0; iteration < 200u; iteration++)
        {

			// Compute Gaussian kernel row
			for (unsigned int m = 0; m < m_dataSize; ++m)
            {
                P[n][m] = exp(-beta * distances[n][m]);
            }
			P[n][n] = std::numeric_limits<double>::min();

			// Compute entropy of current row
            sum_P = std::numeric_limits<double>::min();
			for (unsigned int m = 0; m < m_dataSize; ++m)
            {
                sum_P += P[n][m];
            }

			double H = 0.0; //TODO what is H?
			for (unsigned int m = 0; m < m_dataSize; m++)
            {
                H += beta * (distances[n][m] * P[n][m]);
            }
			H = (H / sum_P) + log(sum_P);

			// Evaluate whether the entropy is within the tolerance level
			double Hdiff = H - log(m_perplexity);
			if (std::abs(Hdiff) < tolerance_threshold)
            {
				break;
			}

			if (Hdiff > 0)
            {
                min_beta = beta;
                if (max_beta == std::numeric_limits<double>::max() || max_beta == std::numeric_limits<double>::lowest())
                {
                    beta *= 2.0;
                }
                else
                {
                    beta = (beta + max_beta) / 2.0;
                }
            }
            else
            {
                max_beta = beta;
                if (min_beta == std::numeric_limits<double>::lowest() || min_beta == std::numeric_limits<double>::max())
                {
                    beta /= 2.0;
                }
                else
                {
                    beta = (beta + min_beta) / 2.0;
                }
            }
		}

		// Row normalize P
		for (unsigned int m = 0; m < m_dataSize; ++m)
        {
            P[n][m] /= sum_P;
        }
	}

    return P;
}

Vector2D<double> TSNE::computeSquaredEuclideanDistance(const Vector2D<double> & points)
{
    auto dimensions = points.width();
    auto number = points.height();

    auto distances = Vector2D<double>(number, number, 0.0);

    for (unsigned int i = 0; i < number; ++i)
    {
        for (unsigned int j = i + 1; j < number; ++j)
        {
            double distance = 0.0;
            for (size_t d = 0; d < dimensions; ++d)
            {
                double diff = points[i][d] - points[j][d];
                distance += diff * diff;
            }

            distances[i][j] = distance;
            distances[j][i] = distance;
        }
    }

    return distances;
}

void TSNE::computeGaussianPerplexity(SparseMatrix & similarities)
{
    assert(m_data.height() == m_dataSize);
    assert(m_data.width() == m_inputDimensions);

	auto K = static_cast<unsigned int>(3 * m_perplexity);

	// Allocate the memory we need
    similarities.rows.resize(m_dataSize + 1);
    similarities.columns.resize(m_dataSize * K);
    similarities.values.resize(m_dataSize * K, 0.0);

    similarities.rows[0] = 0;
	for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        similarities.rows[n + 1] = similarities.rows[n] + K;
    }

	// Build ball tree on data set
	auto vantagePointTree = VantagePointTree();
	auto obj_X = std::vector<DataPoint>(m_dataSize, DataPoint(m_inputDimensions, 0, m_data[0]));
	for (unsigned int n = 0; n < m_dataSize; ++n)
    {
        obj_X[n].index = n;
        obj_X[n].data.assign(m_data[n], m_data[n] + m_inputDimensions);
    }
	vantagePointTree.create(obj_X);

	// Loop over all points to find nearest neighbors
	std::cout << "building vantage point tree..." << std::endl;
	auto indices = std::vector<DataPoint>();
	auto distances = std::vector<double>();
    auto cur_P = std::vector<double>(m_dataSize - 1);
	for (unsigned int n = 0; n < m_dataSize; ++n)
    {
		if (n % 10000 == 0)
        {
            std::cout << " - point " << n << " of " << m_dataSize << std::endl;
        }

		// Find nearest neighbors
		vantagePointTree.search(obj_X[n], K + 1, indices, distances);

		// Initialize some variables for binary search
		double beta = 1.0;
		double min_beta = std::numeric_limits<double>::lowest();
		double max_beta = std::numeric_limits<double>::max();
		double tolerance_threshold = 1e-5;

		// Iterate until we found a good perplexity
		double sum_P = 0.0;
		for (unsigned int iteration = 0; iteration < 200u; iteration++)
        {
			// Compute Gaussian kernel row
			for (unsigned int m = 0; m < K; ++m)
            {
                cur_P[m] = exp(-beta * distances[m + 1] * distances[m + 1]);
            }

			// Compute entropy of current row
			sum_P = std::numeric_limits<double>::min();
			for (unsigned int m = 0; m < K; ++m)
            {
                sum_P += cur_P[m];
            }

			double H = 0.0;
			for (unsigned int m = 0; m < K; ++m)
            {
                H += beta * (distances[m + 1] * distances[m + 1] * cur_P[m]);
            }
			H = (H / sum_P) + log(sum_P);

			// Evaluate whether the entropy is within the tolerance level
			double Hdiff = H - log(m_perplexity);
			if (std::abs(Hdiff) < tolerance_threshold)
            {
				break;
			}
            if (Hdiff > 0)
            {
                min_beta = beta;
                if (max_beta == std::numeric_limits<double>::max()
                    || max_beta == std::numeric_limits<double>::lowest())
                {
                    beta *= 2.0;
                }
                else
                {
                    beta = (beta + max_beta) / 2.0;
                }
            }
            else
            {
                max_beta = beta;
                if (min_beta == std::numeric_limits<double>::lowest()
                    || min_beta == std::numeric_limits<double>::max())
                {
                    beta /= 2.0;
                }
                else
                {
                    beta = (beta + min_beta) / 2.0;
                }
            }
		}

		// Row-normalize current row of P and store in matrix
		for (unsigned int m = 0; m < K; ++m)
        {
            cur_P[m] /= sum_P;
            similarities.columns[similarities.rows[n] + m] = indices[m + 1].index;
            similarities.values[similarities.rows[n] + m] = cur_P[m];
		}
	}
}
