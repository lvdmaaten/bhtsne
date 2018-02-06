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
#include <random>

#include <bhtsne/Vector2D.h>
#include <bhtsne/SparseMatrix.h>
#include <bhtsne/bhtsne_api.h> // generated header for export macros


static inline double sign(double x) { return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0)); }


namespace bhtsne
{

/**
*  @brief
*    Representation of the Barnes-Hut approximation for
*    the t-distributed stochastic neighbor embedding algorithm
*
*    Default parameters are:
*    - randomSeed          random
*    - perplexity          50
*    - gradientAccuracy    0.2
*    - iterations          1000
*    - outputDimensions    2
*    - outputFile          "./result"
*
*    This class follows the method object pattern (SmalltalkBestPracticePatterns, page 34-37).
*/
class BHTSNE_API TSNE
{
public:
    /**
    *  @brief
    *    Constructor with resonable defaults
    *
    *  @see TSNE
    */
    TSNE();

    /**
    *  @brief
    *    Set random seed
    *
    *  @param[in] seed
    *    Random seed
    */
    void setRandomSeed(unsigned long seed);

    /**
    *  @brief
    *    Get perplexity
    *
    *  @return
    *    Perplexity (0..inputDimensions)
    *
    *  @remarks
    *    Perplexity controls how to balance attention between local and global aspects of the
    *    data. The parameter is a guess about the number of close neighbors each point has. Typical
    *    values are between 5 and 50.
    */
    double perplexity() const;

    /**
    *  @brief
    *    Set perplexity
    *
    *  @param[in] perplexity
    *    Perplexity (0..inputDimensions)
    *
    *  @see perplexity()
    */
    void setPerplexity(double perplexity);

    /**
    *  @brief
    *    Get gradient accuracy
    *
    *  @return
    *    Gradient accuracy (0..1)
    *
    *  @remarks
    *    Called theta in other implementations. This value
    *    is used as the width for the gauss sampling kernel.
    */
    double gradientAccuracy() const;

    /**
    *  @brief
    *    Set gradient accuracy
    *
    *  @param[in] accuracy
    *    Gradient accuracy (0..1)
    *
    *  @see gradientAccuracy()
    */
    void setGradientAccuracy(double accuracy);

    /**
    *  @brief
    *    Get number of iterations
    *
    *  @return
    *    Number of iterations
    *
    *  @remarks
    *    This defines how many times the algorithm adjusts the result to minimize the error.
    *    A higher value can improve the result but also increases computation times.
    */
    unsigned int iterations() const;

    /**
    *  @brief
    *    Set number of iterations
    *
    *  @param[in] iterations
    *    Number of iterations (usually larger than 250)
    *
    *  @see iterations()
    */
    void setIterations(unsigned int iterations);

    /**
    *  @brief
    *    Get output dimensionality
    *
    *  @return
    *    Output dimensionality (1 <= outDim < inDim)
    *
    *  @remarks
    *    This defines the dimensionality of the result. Therefore, it should be smaller than the input dimensionality.
    */
    unsigned int outputDimensions() const;

    /**
    *  @brief
    *    Set output dimensionality
    *
    *  @param[in] dimensions
    *    Output dimensionality (1 <= outDim < inDim)
    *
    *  @see outputDimensions()
    */
    void setOutputDimensions(unsigned int dimensions); //2 or 3

    /**
    *  @brief
    *    Get input dimensionality
    *
    *  @return
    *    Input dimensionality (1 <= outDim < inDim)
    *
    *  @remarks
    *    This is the dimensionality of the loaded dataset.
    *    The loading method (e.g. loadCSV()) deduces this value from the given dataset file.
    */
    unsigned int inputDimensions() const;

    /**
    *  @brief
    *    Get number of data points in input
    *
    *  @return
    *    number of data points
    *
    *  @remarks
    *    This is the number of data points in the loaded dataset.
    *    The loading method (e.g. loadCSV()) deduces this value from the given dataset file.
    */
    unsigned int dataSize() const;

    /**
    *  @brief
    *    Get output file
    *
    *  @return
    *    Output file
    *
    *  @remarks
    *    This is the path and basename that is used to create all output files.
    */
    std::string outputFile() const;

    /**
    *  @brief
    *    Set output file
    *
    *  @param[in] file
    *    Output file
    *
    *  @see outputFile()
    */
    void setOutputFile(const std::string & file);


    // load methods---------------------------------------------------------------------------------

    /**
    *  @brief
    *    Loads a csv dataset from a stream
    *
    *  @param[in] file
    *    Input file path
    *
    *  @post
    *    If this fuction returns true, the dataset was loaded and run() can be called.
    *
    *  @remarks
    *    Loads csv data from a stream, expects the same format as loadCSV().
    *    The data should contain x lines with y numerical values each.
    *    The number of samples is set to x and the input dimensionality to y.
    */
    bool loadFromStream(std::istream & stream);


    /**
    *  @brief
    *    Loads a dataset from std::cin
    *
    *  @remarks
    *    Expects the same format as loadCSV()
    */
    bool loadCin();


    /**
    *  @brief
    *    Loads a dataset from a ".dat" file
    *
    *  @param[in] file
    *    Input file path
    *
    *  @post
    *    If this fuction returns true, the dataset was loaded and run() can be called.
    *
    *  @remarks
    *    Loads files with the extension ".dat" as used in the original implementation (https://github.com/lvdmaaten/bhtsne).
    *    In addition to the dataset parameters also
    *    sets (and overrides) gradient accuracy, perplexity, output dimensionality, iterations, and possibly the seed.
    *    The file must contain the following binary information:
    *    - int         number of datapoints
    *    - int         input dimensionality
    *    - double      gradient accuracy
    *    - double      perplexity
    *    - int         output dimensionality
    *    - int         number of iterations
    *    - double...   the data as an interleaved buffer (number of datapoints * input dimensionality)
    *    - int         random seed (optional)
    */
    bool loadLegacy(const std::string & file);

    /**
    *  @brief
    *    Loads a dataset from a ".csv" file
    *
    *  @param[in] file
    *    Input file path
    *
    *  @post
    *    If this fuction returns true, the dataset was loaded and run() can be called.
    *
    *  @remarks
    *    Loads files with the extension ".csv".
    *    The file should contain x lines with y numerical values each.
    *    The number of samples is set to x and the input dimensionality to y.
    */
    bool loadCSV(const std::string & file);

    /**
    *  @brief
    *    Loads a dataset from a ".tsne" file
    *
    *  @param[in] file
    *    Input file path
    *
    *  @post
    *    If this fuction returns true, the dataset was loaded and run() can be called.
    *
    *  @remarks
    *    Loads files with the extension ".tsne".
    *    The file must contain the following binary information:
    *    - int         number of datapoints
    *    - int         input dimensionality
    *    - double...   the data as an interleaved buffer (number of datapoints * input dimensionality)
    */
    bool loadTSNE(const std::string & file);


    //run method------------------------------------------------------------------------------------

    /**
    *  @brief
    *    Runs the algorithm
    *
    *  @pre
    *    A dataset must be loaded by any of the "load" functions.
    *
    *  @post
    *    The result is computed and ready to be saved by any of the "save" functions.
    *
    *  @remarks
    *    This function runs the barnes hut implementation of the T-SNE algorithm.
    */
    void run();


    //save methods----------------------------------------------------------------------------------

    /**
    *  @brief
    *    Pushes csv result to stream
    *
    *  @pre
    *    The algorithm must have ran (i.e. run() was called).
    *
    *  @remarks
    *    Pushes the csv result to the given stream.
    */
    void saveToStream(std::ostream & stream);

    /**
    *  @brief
    *    Pushes the csv result to std::cout
    *
    *  @pre
    *    The algorithm must have ran (i.e. run() was called).
    *
    *  @remarks
    *    Pushes the csv result to std::cout.
    */
    void saveToCout();

    /**
    *  @brief
    *    Saves the result in a ".csv" file
    *
    *  @pre
    *    The algorithm must have ran (i.e. run() was called).
    *
    *  @remarks
    *    Saves the result of the computation to "<outputFile>.csv".
    */
    void saveCSV();

    /**
    *  @brief
    *    Saves the result in a ".dat" file
    *
    *  @pre
    *    The algorithm must have ran (i.e. run() was called).
    *
    *  @remarks
    *    Saves the result of the computation to "<outputFile>.dat".
    *    This format was used in the original implementation (https://github.com/lvdmaaten/bhtsne).
    */
    void saveLegacy();

    /**
    *  @brief
    *    Saves the result in a ".svg" file
    *
    *  @pre
    *    - the algorithm must have ran (i.e. run() was called)
    *    - output dimensionality = 2
    *
    *  @remarks
    *    Saves the result of the computation to "<outputFile>.svg".
    */
    void saveSVG();

protected:
    void runApproximation();
    void runExact();

    template<unsigned int D>
    Vector2D<double> computeGradient(SparseMatrix & similarities);
    Vector2D<double> computeGradientExact(const Vector2D<double> & Perplexity);
    template<unsigned int D>
    double evaluateError(SparseMatrix & similarities);
    double evaluateErrorExact(const Vector2D<double> & Perplexity);
    void computeGaussianPerplexity(SparseMatrix & similarities);
    Vector2D<double> computeGaussianPerplexityExact();

    // params
    double       m_perplexity;         ///< balance local/global data aspects, see documentation of perplexity()
    double       m_gradientAccuracy;   ///< used as the width for the gauss sampling kernel
    unsigned int m_iterations;         ///< defines how many iterations the algorithm does in run()

    // dataset
    unsigned int m_outputDimensions;   ///< dimensionality of the result
    unsigned int m_inputDimensions;    ///< dimensionality of the input; set during load
    unsigned int m_dataSize;           ///< size of data; set during load
	Vector2D<double> m_data;           ///< loaded data
    std::mt19937 m_gen;                ///< random number generator

    // output
    std::string  m_outputFile;         ///< path and basename used to create output files
	Vector2D<double> m_result;         ///< computation results



    //helper
    static Vector2D<double> computeSquaredEuclideanDistance(const Vector2D<double> & points);
    void symmetrizeMatrix(SparseMatrix & similarities);
    static void zeroMean(Vector2D<double>& points);
    static void normalize(Vector2D<double>& vec);
    double gaussNumber();
};


} // bhtsne
