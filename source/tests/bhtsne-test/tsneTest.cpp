#include <bhtsne/tsne.h>
#include <gtest/gtest.h>
#include <string>
#include <cstdio>
#include <fstream>
#include <initializer_list>

class PublicTSNE : public bhtsne::TSNE
{
public:
    std::vector<std::vector<double>> data() 
    {
        return m_data;
    };

    void setData(std::vector<std::vector<double>> data)
    {
        m_data = data;
    };

    std::vector<std::vector<double>> result()
    {
        return m_result;
    };

    void setResult(std::vector<std::vector<double>> result)
    {
        m_result = result;
    };

    void setInputDimensions(int dimensions)
    {
        m_inputDimensions = dimensions;
    }

    double firstGaussNumber()
    {
        return gaussNumber();
    }
};

class BinaryWriter
{
private:
    std::ofstream* filestream;

public:
    BinaryWriter() = default;

    explicit BinaryWriter(std::ofstream& f)
    {
        filestream = &f;
    }

    template <typename T>
    BinaryWriter& operator<<(const T& value)
    {
        filestream->write(reinterpret_cast<const char*>(&value), sizeof(value));
        return *this;
    }
};

class TsneTest : public testing::Test
{
protected:
    TsneTest()
    {
        m_tsne = PublicTSNE();
        tempfile = std::tmpnam(nullptr);
    }

    void createTempfile()
    {
        filestream.open(tempfile, std::ios::out | std::ios::binary | std::ios::trunc);
        writer = BinaryWriter(filestream);
    }

    void removeTempfile()
    {
        remove(tempfile.c_str());
    }

    double firstGaussNumber(int seed)
    {
        auto gen = std::mt19937(seed);
        auto dist = std::normal_distribution<>(0,2);
        return dist(gen);
    }

    PublicTSNE m_tsne;
    std::string tempfile;
    std::ofstream filestream;
    BinaryWriter writer;
};

TEST(SanityChecks, Equality)
{
    EXPECT_EQ((unsigned int) 0, 0);
    EXPECT_EQ((unsigned int) 1, 1);
}   

TEST_F(TsneTest, DefaultValues)
{
    EXPECT_EQ(firstGaussNumber(0), m_tsne.firstGaussNumber());
    EXPECT_EQ(50.0, m_tsne.perplexity());
    EXPECT_EQ(0.2, m_tsne.gradientAccuracy());
    EXPECT_EQ(1000, m_tsne.iterations());
    EXPECT_EQ(2, m_tsne.outputDimensions());
    EXPECT_EQ(0, m_tsne.dataSize());
    EXPECT_EQ("result", m_tsne.outputFile());
}

TEST_F(TsneTest, RandomSeed)
{
    m_tsne.setRandomSeed(1);
    EXPECT_EQ(firstGaussNumber(1), m_tsne.firstGaussNumber());
    m_tsne.setRandomSeed(2);
    EXPECT_EQ(firstGaussNumber(2), m_tsne.firstGaussNumber());
}

TEST_F(TsneTest, Perplexity)
{
    m_tsne.setPerplexity(1.0);
    EXPECT_EQ(1.0, m_tsne.perplexity());
    m_tsne.setPerplexity(2.0);
    EXPECT_EQ(2.0, m_tsne.perplexity());
}

TEST_F(TsneTest, GradientAccuracy)
{
    m_tsne.setGradientAccuracy(1.0);
    EXPECT_EQ(1.0, m_tsne.gradientAccuracy());
    m_tsne.setGradientAccuracy(2.0);
    EXPECT_EQ(2.0, m_tsne.gradientAccuracy());
}

TEST_F(TsneTest, Iterations)
{
    m_tsne.setIterations(1);
    EXPECT_EQ(1, m_tsne.iterations());
    m_tsne.setIterations(2);
    EXPECT_EQ(2, m_tsne.iterations());
}

TEST_F(TsneTest, OutputDimensions)
{
    m_tsne.setOutputDimensions(1);
    EXPECT_EQ(1, m_tsne.outputDimensions());
    m_tsne.setOutputDimensions(2);
    EXPECT_EQ(2, m_tsne.outputDimensions());
}

TEST_F(TsneTest, InputDimensions)
{
    // only test accessor - setting is tested alongside the respective load method
    m_tsne.setInputDimensions(1);
    EXPECT_EQ(1, m_tsne.inputDimensions());
    m_tsne.setInputDimensions(2);
    EXPECT_EQ(2, m_tsne.inputDimensions());
}

TEST_F(TsneTest, DataSize)
{
    m_tsne.setDataSize(1);
    EXPECT_EQ(1, m_tsne.dataSize());
    m_tsne.setDataSize(2);
    EXPECT_EQ(2, m_tsne.dataSize());
}

TEST_F(TsneTest, OutputFile)
{
    m_tsne.setOutputFile("foo.dat");
    EXPECT_EQ("foo.dat", m_tsne.outputFile());
    m_tsne.setOutputFile("bar.dat");
    EXPECT_EQ("bar.dat", m_tsne.outputFile());
}

TEST_F(TsneTest, LoadLegacy)
{
    int dataSize = 1;
    int inputDimensions = 1;
    double gradientAccuracy = 0.5;
    double perplexity = 25.0;
    int outputDimensions = 1;
    int iterations = 100;
    int randomSeed = 42;
    double data = 42.0;

    createTempfile();
    writer << dataSize << inputDimensions << gradientAccuracy << perplexity << outputDimensions << iterations;
    writer << data;
    writer << randomSeed;
    filestream.flush();
    filestream.close();

    EXPECT_TRUE(m_tsne.loadLegacy(tempfile));
    EXPECT_EQ(dataSize, m_tsne.dataSize());
    EXPECT_EQ(inputDimensions, m_tsne.inputDimensions());
    EXPECT_EQ(gradientAccuracy, m_tsne.gradientAccuracy());
    EXPECT_EQ(perplexity, m_tsne.perplexity());
    EXPECT_EQ(outputDimensions, m_tsne.outputDimensions());
    EXPECT_EQ(iterations, m_tsne.iterations());
    EXPECT_EQ(firstGaussNumber(randomSeed), m_tsne.firstGaussNumber());
    EXPECT_EQ(data, m_tsne.data()[0][0]);

    removeTempfile();
}

TEST_F(TsneTest, LoadCSV)
{
    auto data = std::vector<std::vector<double>>{ { 1.0, 2.0, 3.0 },{ 4.0, 5.0, 6.0 } };

    createTempfile();
    for (auto sample : data)
    {
        for (auto value : sample)
        {
            filestream << value << ',';
        }
        filestream.seekp(-1, std::ios_base::cur);
        filestream << std::endl;
    }

    filestream.flush();
    filestream.close();

    EXPECT_TRUE(m_tsne.loadCSV(tempfile));
    EXPECT_EQ(2, m_tsne.dataSize());
    EXPECT_EQ(3, m_tsne.inputDimensions());
    for (auto i = 0; i < 2; i++)
    {
        for (auto j = 0; j < 3; j++)
        {
            EXPECT_EQ(data[i][j], m_tsne.data()[i][j]);
        }
    }

    removeTempfile();
}

TEST_F(TsneTest, LoadTSNE)
{
    int dataSize = 1;
    int inputDimensions = 1;
    double data = 42.0;

    createTempfile();
    writer << dataSize << inputDimensions;
    writer << data;
    filestream.flush();
    filestream.close();

    EXPECT_TRUE(m_tsne.loadTSNE(tempfile));
    EXPECT_EQ(1, m_tsne.dataSize());
    EXPECT_EQ(1, m_tsne.inputDimensions());
    EXPECT_EQ(data, m_tsne.data()[0][0]);

    removeTempfile();
}

TEST_F(TsneTest, Run)
{
    m_tsne.setDataSize(2);
    m_tsne.setInputDimensions(1);
    m_tsne.setGradientAccuracy(0.2);
    m_tsne.setPerplexity(0.1);
    m_tsne.setOutputDimensions(1);
    m_tsne.setIterations(100);
    m_tsne.setRandomSeed(42);
    m_tsne.setData(std::vector<std::vector<double>>{ { 42.0 }, { 17.0 } });

    EXPECT_NO_THROW(m_tsne.run());
    std::vector<unsigned long long> expected = { 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    for (auto i = 0; i < 2; i++)
    {
        EXPECT_EQ(*reinterpret_cast<double*>(&expected[i]), m_tsne.result()[i][0]);
    }
}

TEST_F(TsneTest, SaveLegacy)
{
    std::vector<unsigned long long> expected = { 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    std::vector<double> expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });
    
    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setOutputFile(tempfile); 
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0] }, { expectedDouble[1] }});
    EXPECT_NO_THROW(m_tsne.saveLegacy());
    // check file exists and has right size
    std::ifstream result;
    EXPECT_NO_THROW(result.open(tempfile + ".dat", std::ios::in | std::ios::binary | std::ios::ate));
    EXPECT_TRUE(result.is_open());
    EXPECT_EQ(48, result.tellg());
    // check values in file
    result.seekg(0);
    int dataSize;
    int outputDimensions;
    std::vector<double> data = std::vector<double>(2);
    std::vector<int> landmarks = std::vector<int>(2);
    std::vector<double> costs = std::vector<double>(2);
    result.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));
    result.read(reinterpret_cast<char*>(&outputDimensions), sizeof(outputDimensions));
    result.read(reinterpret_cast<char*>(data.data()), 2 * sizeof(double));
    result.read(reinterpret_cast<char*>(landmarks.data()), 2 * sizeof(int));
    result.read(reinterpret_cast<char*>(costs.data()), 2 * sizeof(double));
    EXPECT_EQ(2, dataSize);
    EXPECT_EQ(1, outputDimensions);
    for (auto i = 0; i < 2; i++)
    {
        EXPECT_EQ(expectedDouble[i], data[i]);
        EXPECT_EQ(i, landmarks[i]);
        EXPECT_EQ(0, costs[i]);
    }

    result.close();
    remove((tempfile + ".dat").c_str());
}

TEST_F(TsneTest, SaveToStream)
{
    std::vector<unsigned long long> expected = { 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    std::vector<double> expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });

    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0] }, { expectedDouble[1] }});
    std::ostringstream result;
    EXPECT_NO_THROW(m_tsne.saveToStream(result));
    // check values in stream
    std::ostringstream expectedOut;
    expectedOut << std::setprecision(6);
    expectedOut << expectedDouble[0] << std::endl << expectedDouble[1] << std::endl;

    EXPECT_EQ(expectedOut.str(), result.str());
}

TEST_F(TsneTest, SaveCSV)
{
    std::vector<unsigned long long> expected = { 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    std::vector<double> expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });

    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0] }, { expectedDouble[1] }});
    m_tsne.setOutputFile(tempfile);
    EXPECT_NO_THROW(m_tsne.saveCSV());
    // check file exists and has right size
    std::ifstream result;
    EXPECT_NO_THROW(result.open(tempfile + ".csv", std::ios::in | std::ios::ate));
    EXPECT_TRUE(result.is_open());
    // file size may change based on current os (\n vs \r\n)
    EXPECT_LE(17, result.tellg());
    EXPECT_GE(19, result.tellg());
    // check values in file
    std::ostringstream expectedOut;
    expectedOut << std::setprecision(6);
    expectedOut << expectedDouble[0] << std::endl << expectedDouble[1] << std::endl;
    std::istringstream expectedIn(expectedOut.str());

    result.seekg(0);
    std::string actualLine;
    std::string expectedLine;
    while ((bool)std::getline(result, actualLine) & (bool)std::getline(expectedIn, expectedLine))
    {
        EXPECT_EQ(expectedLine, actualLine);
    }

    EXPECT_TRUE(result.eof());
    EXPECT_TRUE(expectedIn.eof());

    result.close();
    EXPECT_EQ(0, remove((tempfile + ".csv").c_str()));
}

TEST_F(TsneTest, SaveSVG)
{
    std::vector<unsigned long long> expected = { 0xC07DAC741DC680D4, 0xC07DAC741DC680D4, 0x407DAC741DC680D4, 0x407DAC741DC680D4 };
    std::vector<double> expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });

    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0], expectedDouble[1] }, { expectedDouble[2], expectedDouble[3] }});
    m_tsne.setOutputFile(tempfile);
    EXPECT_NO_THROW(m_tsne.saveSVG());
    // check file exists
    std::ifstream result;
    EXPECT_NO_THROW(result.open(tempfile + ".svg", std::ios::in));
    EXPECT_TRUE(result.is_open());
    // check values in file
    double expectedRadius = 0.5;
    double expectedViewBoxMin = expectedDouble[0] - expectedRadius;
    double expectedViewBoxSize = 2 * -expectedViewBoxMin;
    std::ostringstream expectedOut;
    expectedOut << std::setprecision(6);
    expectedOut << "<?xml version='1.0' encoding='UTF-8' ?>\n";
    expectedOut << "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='600' height='600' viewBox='" << expectedViewBoxMin << " " << expectedViewBoxMin << " " << expectedViewBoxSize << " " << expectedViewBoxSize << "'>\n";
    expectedOut << "<circle cx='" << expectedDouble[0] << "' cy='" << expectedDouble[1] << "' fill='black' r='" << expectedRadius << "' stroke='none' opacity='0.5'/>\n";
    expectedOut << "<circle cx='" << expectedDouble[2] << "' cy='" << expectedDouble[3] << "' fill='black' r='" << expectedRadius << "' stroke='none' opacity='0.5'/>\n";
    expectedOut << "</svg>\n";
    std::istringstream expectedIn(expectedOut.str());

    std::string actualLine;
    std::string expectedLine;
    while ((bool)std::getline(result, actualLine) & (bool)std::getline(expectedIn, expectedLine))
    {
        EXPECT_EQ(expectedLine, actualLine);
    }

    EXPECT_TRUE(result.eof());
    EXPECT_TRUE(expectedIn.eof());

    result.close();
    EXPECT_EQ(0, remove((tempfile + ".svg").c_str()));
}
