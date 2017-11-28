#include <bhtsne/tsne.h>
#include <gtest/gtest.h>
#include <string>
#include <cstdio>
#include <fstream>

class TsneTest : public testing::Test
{
protected:
    TsneTest()
    {
        m_tnse = bhtsne::TSNE();
    }
    /*virtual void SetUp()
    {
    m_tnse = bhtsne::TSNE();
    }*/

    void createTempfile()
    {
        std::string tempfile = std::tmpnam(nullptr);
        std::cout << tempfile << std::endl;
        std::ofstream filestream;
        filestream.open(tempfile, std::ios::out | std::ios::binary | std::ios::trunc);
    }

    void removeTempfile()
    {
        filestream.close();
        remove(tempfile.c_str());
    }

    bhtsne::TSNE m_tnse;
    std::string tempfile;
    std::ofstream filestream;
};

TEST(SanityChecks, Equality)
{
    EXPECT_EQ((unsigned int) 0, 0);
    EXPECT_EQ((unsigned int) 1, 1);
}   

TEST_F(TsneTest, DefaultValues)
{
    FAIL();
}

TEST_F(TsneTest, RandomSeed)
{
    m_tnse.setRandomSeed(1);
    EXPECT_EQ(m_tnse.randomSeed(), 1);
    m_tnse.setRandomSeed(2);
    EXPECT_EQ(m_tnse.randomSeed(), 2);
}

TEST_F(TsneTest, Perplexity)
{
    m_tnse.setPerplexity(1.0);
    EXPECT_EQ(m_tnse.perplexity(), 1.0);
    m_tnse.setPerplexity(2.0);
    EXPECT_EQ(m_tnse.perplexity(), 2.0);
}

TEST_F(TsneTest, GradientAccuracy)
{
    m_tnse.setGradientAccuracy(1.0);
    EXPECT_EQ(m_tnse.gradientAccuracy(), 1.0);
    m_tnse.setGradientAccuracy(2.0);
    EXPECT_EQ(m_tnse.gradientAccuracy(), 2.0);
}

TEST_F(TsneTest, Iterations)
{
    m_tnse.setIterations(1);
    EXPECT_EQ(m_tnse.iterations(), 1);
    m_tnse.setIterations(2);
    EXPECT_EQ(m_tnse.iterations(), 2);
}

TEST_F(TsneTest, OutputDimensions)
{
    m_tnse.setOutputDimensions(1);
    EXPECT_EQ(m_tnse.outputDimensions(), 1);
    m_tnse.setOutputDimensions(2);
    EXPECT_EQ(m_tnse.outputDimensions(), 2);
}

TEST_F(TsneTest, InputDimensions)
{
    FAIL();
}

TEST_F(TsneTest, NumberOfSamples)
{
    m_tnse.setNumberOfSamples(1);
    EXPECT_EQ(m_tnse.getNumberOfSamples(), 1);
    m_tnse.setNumberOfSamples(2);
    EXPECT_EQ(m_tnse.getNumberOfSamples(), 2);
}

TEST_F(TsneTest, OutputFile)
{
    m_tnse.setOutputFile("foo.dat");
    EXPECT_EQ(m_tnse.outputFile(), "foo.dat");
    m_tnse.setOutputFile("bar.dat");
    EXPECT_EQ(m_tnse.outputFile(), "bar.dat");
}

TEST_F(TsneTest, LoadLegacy)
{
    //createTempfile();
    //filestream << "text in tempfile";
    // run load and chceck if read data is correct
    //removeTempfile();
    FAIL();
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
    FAIL();
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
