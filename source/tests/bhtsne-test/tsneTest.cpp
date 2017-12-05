#include <bhtsne/tsne.h>
#include <gtest/gtest.h>
#include <string>
#include <cstdio>
#include <fstream>
#include <initializer_list>

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
        m_tsne = bhtsne::TSNE();
    }

    void createTempfile()
    {
        tempfile = std::tmpnam(nullptr);
        std::cout << tempfile << std::endl;
        filestream.open(tempfile, std::ios::out | std::ios::binary | std::ios::trunc);
        writer = BinaryWriter(filestream);
    }

    void createTempfileLegacy(int dataSize, int inputDimensions, double gradientAccuracy, double perplexity, int outputDimensions, int iterations, int randomSeed, std::initializer_list<double> data)
    {
        createTempfile();
        writer << dataSize << inputDimensions << gradientAccuracy << perplexity << outputDimensions << iterations;
        for (auto value : data)
        {
            writer << value;
        }
        writer << randomSeed;
        filestream.flush();
        filestream.close();
    }

    void createTempfileTSNE(int dataSize, int inputDimensions, std::initializer_list<double> data)
    {
        createTempfile();
        writer << dataSize << inputDimensions;
        for (auto value : data)
        {
            writer << value;
        }
        filestream.flush();
        filestream.close();
    }

    void createTempfileCSV(int inputDimensions, std::initializer_list<double> data)
    {
        createTempfile();
        int dimension = 0;
        for (auto value : data)
        {
            filestream << value;
            if (++dimension % inputDimensions == 0)
            {
                filestream << std::endl;
            }
            else
            {
                filestream << ',';
            }
        }
        filestream.flush();
        filestream.close();
    }

    void removeTempfile()
    {
        remove(tempfile.c_str());
    }

    bhtsne::TSNE m_tsne;
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
    EXPECT_EQ(0, m_tsne.randomSeed());
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
    EXPECT_EQ(1, m_tsne.randomSeed());
    m_tsne.setRandomSeed(2);
    EXPECT_EQ(2, m_tsne.randomSeed());
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
    createTempfileLegacy(1, 1, 0.5, 25.0, 1, 100, 42, { 42.0 });

    EXPECT_TRUE(m_tsne.loadLegacy(tempfile));
    EXPECT_EQ(1, m_tsne.inputDimensions());

    removeTempfile();

    createTempfileLegacy(1, 3, 0.5, 25.0, 1, 100, 42, { 42.0, 0.0, 0.0 });

    EXPECT_TRUE(m_tsne.loadLegacy(tempfile));
    EXPECT_EQ(3, m_tsne.inputDimensions());

    removeTempfile();
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
    createTempfileLegacy(1, 1, 0.5, 25.0, 1, 100, 42, { 42.0 });

    EXPECT_TRUE(m_tsne.loadLegacy(tempfile));
    EXPECT_EQ(1, m_tsne.dataSize());
    EXPECT_EQ(1, m_tsne.inputDimensions());
    EXPECT_EQ(0.5, m_tsne.gradientAccuracy());
    EXPECT_EQ(25.0, m_tsne.perplexity());
    EXPECT_EQ(1, m_tsne.outputDimensions());
    EXPECT_EQ(100, m_tsne.iterations());
    EXPECT_EQ(42, m_tsne.randomSeed());

    removeTempfile();
}

TEST_F(TsneTest, LoadCSV)
{
    createTempfileCSV(3, { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });

    EXPECT_TRUE(m_tsne.loadCSV(tempfile));
    EXPECT_EQ(2, m_tsne.dataSize());
    EXPECT_EQ(3, m_tsne.inputDimensions());

    removeTempfile();
}

TEST_F(TsneTest, LoadTSNE)
{
    createTempfileTSNE(1, 1, { 42.0 });

    EXPECT_TRUE(m_tsne.loadTSNE(tempfile));
    EXPECT_EQ(1, m_tsne.dataSize());
    EXPECT_EQ(1, m_tsne.inputDimensions());

    removeTempfile();
}

TEST_F(TsneTest, Run)
{
    createTempfileLegacy(2, 1, 0.1, 0.01, 1, 100, 42, { 42.0, 17.0 });

    EXPECT_TRUE(m_tsne.loadLegacy(tempfile));
    EXPECT_NO_THROW(m_tsne.run());

    removeTempfile();
}

TEST_F(TsneTest, SaveLegacy)
{
    createTempfileLegacy(2, 1, 0.1, 0.01, 1, 100, 42, { 42.0, 17.0 });
    // test load/run/save
    EXPECT_TRUE(m_tsne.loadLegacy(tempfile));
    EXPECT_NO_THROW(m_tsne.run());
    m_tsne.setOutputFile(tempfile);
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
    std::vector<unsigned long long> expected = { 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    for (auto i = 0; i < 2; i++)
    {
        EXPECT_EQ(*reinterpret_cast<double*>(&expected[i]), data[i]);
        EXPECT_EQ(i, landmarks[i]);
        EXPECT_EQ(0, costs[i]);
    }

    removeTempfile();
    remove((tempfile + ".dat").c_str());
}

TEST_F(TsneTest, SaveCSV)
{
    FAIL();
}

TEST_F(TsneTest, SaveSVG)
{
    FAIL();
}
