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
    BinaryWriter() {}

    BinaryWriter(std::ofstream& f)
    {
        filestream = &f;
    }

    template <typename T>
    BinaryWriter& operator<<(const T& value)
    {
        filestream->write(reinterpret_cast<const char*>(&value), sizeof(value));
        return *this;
    };
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
        filestream;
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
    FAIL();
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
    FAIL();
}

TEST_F(TsneTest, LoadTSNE)
{
    FAIL();
}

TEST_F(TsneTest, Run)
{
    createTempfileLegacy(3, 1, 0.2, 0.1, 1, 100, 42, { 42.0, 17.0, 1.0 });

    EXPECT_TRUE(m_tsne.loadLegacy(tempfile));
    EXPECT_NO_THROW(m_tsne.run());

    removeTempfile();
}

TEST_F(TsneTest, SaveLegacy)
{
    FAIL();
}

TEST_F(TsneTest, SaveCSV)
{
    FAIL();
}

TEST_F(TsneTest, SaveSVG)
{
    FAIL();
}
