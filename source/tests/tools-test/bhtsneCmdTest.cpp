#include <bhtsne/tsne.h>
#include <gtest/gtest.h>
#include <cstring>
#include "../../tools/bhtsne_cmd/ArgumentParser.cpp" //needed because we are not linking to a ArgumentParser definition
#include "../../tools/bhtsne_cmd/CommandlineOptions.cpp"

class bhtsneCmdTest : public testing::Test
{
protected:
    bhtsneCmdTest()
    {
        m_tsne = bhtsne::TSNE();
    }
    /*virtual void SetUp()
    {
    m_tsne = bhtsne::TSNE();
    }*/

    bhtsne::TSNE m_tsne;
};

TEST(SanityChecks, Equality)
{
    EXPECT_EQ((unsigned int) 0, 0);
    EXPECT_EQ((unsigned int) 1, 1);
}

TEST_F(bhtsneCmdTest, CommandLineArgsParsing)
{
    m_tsne.setPerplexity(30.123);
    m_tsne.setGradientAccuracy(1.123);
    m_tsne.setIterations(5123);
    m_tsne.setNumberOfSamples(2123);
    m_tsne.setOutputDimensions(3);
    m_tsne.setOutputFile("result123.csv");
    m_tsne.setRandomSeed(123);

    char argString[] = "./bhtsne_cmd --perplexity 40.123 --gradient-accuracy 2.123 --iterations 4123 --number-of-samples 3123 --output-dimensions 2 --output-file resultbla123.dat --random-seed 321 input_file.dat";
    auto argv = std::vector<char *>();
    auto token = strtok(argString, " ");
    while (token != nullptr) {
        argv.push_back(token);
        token = strtok(nullptr, " ");
    }

    applyCommandlineOptions(m_tsne, static_cast<int>(argv.size()), argv.data());
    free(argString);
    EXPECT_EQ(40.123, m_tsne.perplexity()) << "perplexity was not set correctly via commandline option";
    EXPECT_EQ(2.123, m_tsne.gradientAccuracy()) << "gradient-accuracy was not set correctly via commandline option";
    EXPECT_EQ(4123, m_tsne.iterations()) << "iterations was not set correctly via commandline option";
    EXPECT_EQ(3123, m_tsne.getNumberOfSamples()) << "number-of-samples was not set correctly via commandline option";
    EXPECT_EQ(2, m_tsne.outputDimensions()) << "output-dimensions was not set correctly via commandline option";
    EXPECT_EQ("resultbla123.dat", m_tsne.outputFile()) << "output-file was not set correctly via commandline option";
    EXPECT_EQ(321, m_tsne.randomSeed()) << "random-seed was not set correctly via commandline option";
}
