#include <bhtsne/TSNE.h>
#include <gtest/gtest.h>
#include <cstring>
#include "../../tools/bhtsne_cmd/ArgumentParser.cpp" //needed because we are not linking to an ArgumentParser definition
#include "../../tools/bhtsne_cmd/CommandlineOptions.cpp"

void parseArguments(cppassist::ArgumentParser & parsedArguments, std::string argString) {
    auto argv = std::vector<char *>();
    auto token = strtok(&argString[0], " ");
    while (token != nullptr) {
        argv.push_back(token);
        token = strtok(nullptr, " ");
    }
    parsedArguments.parse(static_cast<int>(argv.size()), argv.data());
}

class BhtsneCmdTest : public testing::Test
{
protected:
    BhtsneCmdTest()
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

TEST_F(BhtsneCmdTest, SettingCommandLineParameters)
{
    m_tsne.setPerplexity(30.123);
    m_tsne.setGradientAccuracy(1.123);
    m_tsne.setIterations(5123);
    // m_tsne.setDataSize(2123);
    m_tsne.setOutputDimensions(3);
    m_tsne.setOutputFile("result123.csv");
    m_tsne.setRandomSeed(123);

    auto parsedArguments = cppassist::ArgumentParser();
    parseArguments(parsedArguments, "./bhtsne_cmd "
                           "--perplexity 40.123 "
                           "--gradient-accuracy 2.123 "
                           "--iterations 4123 "
                           // "--data-size 3123 "
                           "--output-dimensions 2 "
                           "--output-file another_result.dat "
                           "--random-seed 321 input_file.dat");

    applyCommandlineOptions(m_tsne, parsedArguments.options());

    EXPECT_EQ(40.123, m_tsne.perplexity()) << "perplexity was not set correctly via commandline option";
    EXPECT_EQ(2.123, m_tsne.gradientAccuracy()) << "gradient-accuracy was not set correctly via commandline option";
    EXPECT_EQ(4123, m_tsne.iterations()) << "iterations was not set correctly via commandline option";
    // EXPECT_EQ(3123, m_tsne.dataSize()) << "number-of-samples was not set correctly via commandline option";
    EXPECT_EQ(2, m_tsne.outputDimensions()) << "output-dimensions was not set correctly via commandline option";
    EXPECT_EQ("another_result.dat", m_tsne.outputFile()) << "output-file was not set correctly via commandline option";
}

TEST_F(BhtsneCmdTest, SettingCommandLineOptions)
{
    auto parsedArguments = cppassist::ArgumentParser();
    parseArguments(parsedArguments, "./bhtsne_cmd");

    EXPECT_TRUE(!parsedArguments.isSet("-legacy")
                && !parsedArguments.isSet("-svg")
                && !parsedArguments.isSet("-csv")) << "did not recognize that no output option was set";

    parseArguments(parsedArguments, "./bhtsne_cmd -svg -csv -legacy");

    EXPECT_TRUE(parsedArguments.isSet("-svg")) << "svg output was not set correctly via commandline option";
    EXPECT_TRUE(parsedArguments.isSet("-csv")) << "csv output was not set correctly via commandline option";
    EXPECT_TRUE(parsedArguments.isSet("-legacy")) << "legacy output was not set correctly via commandline option";
}
