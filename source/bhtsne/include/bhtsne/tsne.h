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

#pragma once


#include <string>
#include <vector>

#include <bhtsne/bhtsne_api.h> // generated header for export macros


static inline double sign(double x) { return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0)); }


namespace bhtsne
{


class BHTSNE_API TSNE
{
public:

    /**
    *  @brief
    *    Constructor. Sets resonable defaults.
    *    TODO: mention default values.
    */
    TSNE();

    /**
    *  @brief
    *    Get seed for random generator
    *
    *  @return
    *    Currently used random seed.
    */
    int randomSeed() const;

    /**
    *  @brief
    *    Set option value
    *
    *  @param[in] seed
    *    an int
    */
    void setRandomSeed(int seed);

    double perplexity() const;
    void setPerplexity(double perplexity);

    /**
    *  @brief
    *    Get gradient accuracy
    *
    *  @remarks
    *    Called theta in other implementations
    */
    double gradientAccuracy() const;

    /**
    *  @brief
    *    Set gradient accuracy
    *
    *  @remarks
    *    Called theta in other implementations
    */
    void setGradientAccuracy(double accuracy);

    unsigned int iterations() const;
    void setIterations(unsigned int iterations);

    unsigned int outputDimensions() const;
    void setOutputDimensions(unsigned int dimensions);

    //no setter coz from dataset deduced
    unsigned int inputDimensions() const;

    unsigned int getNumberOfSamples() const; //TODO: rename (w/o get)
    void setNumberOfSamples(unsigned int value);

    std::string outputFile() const;
    void setOutputFile(const std::string& file);


    bool loadLegacy(std::string file);
    bool loadCSV(std::string file);
    bool loadTSNE(std::string file);

    void run();

    void saveLegacy();
    void saveCSV();
    void saveSVG();

    //TODO: remove legacy stuff
    void run(double* X, int D, double* Y, int no_dims, double perplexity, double theta, int rand_seed,
             bool skip_random_init, int max_iter=1000, int stop_lying_iter=250, int mom_switch_iter=250);
    bool load_data(double** data, int* d, int* no_dims, double* theta, double* perplexity, int* rand_seed, int* max_iter);
    void save_data(double* data, int* landmarks, double* costs, int n, int d);
    static void symmetrizeMatrix(unsigned int** row_P, unsigned int** col_P, double** val_P, int N); // should be static!



private:
    void computeGradient(unsigned int* inp_row_P, unsigned int* inp_col_P, double* inp_val_P, double* Y, int D, double* dC, double theta);
    void computeExactGradient(double* P, double* Y, int D, double* dC);
    double evaluateError(double* P, double* Y, int D);
    double evaluateError(unsigned int* row_P, unsigned int* col_P, double* val_P, double* Y, int D, double theta);
    void zeroMean(double* X, int N, int D); //static?
    void computeGaussianPerplexity(double* X, int N, int D, double* P, double perplexity);
    void computeGaussianPerplexity(double* X, int N, int D, unsigned int** _row_P, unsigned int** _col_P, double** _val_P, double perplexity, int K);
    void computeSquaredEuclideanDistance(double* X, int N, int D, double* DD); //static?
    double randn();


    //params
    int          m_randomSeed;      ///< TODO comment
    double       m_perplexity;      ///< TODO comment
    double       m_gradientAccuracy;///< TODO comment
    unsigned int m_iterations;      ///< TODO comment

    //dataset
    unsigned int m_outputDimensions;    ///< TODO comment
    unsigned int m_inputDimensions;     ///< TODO comment
    unsigned int m_numberOfSamples;     ///< TODO comment and rename
	std::vector<std::vector<double>> m_data;

    std::string  m_outputFile;          ///< TODO comment
};


}; //bhtsne
