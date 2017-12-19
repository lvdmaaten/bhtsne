#include <bhtsne/tsne.h>
#include <gtest/gtest.h>
#include <string>
#include <cstdio>
#include <fstream>
#include <initializer_list>

inline bool doubleEquals(double a, double b, double epsilon = 0.0001) {
    return std::fabs(a - b) < epsilon;
}

class PublicTSNE : public bhtsne::TSNE
{
public:
    auto data()
    {
        return m_data;
    };

    auto setData(std::vector<std::vector<double>> data)
    {
        m_data = data;
    };

    auto result()
    {
        return m_result;
    };

    auto setResult(std::vector<std::vector<double>> result)
    {
        m_result = result;
    };

    auto setInputDimensions(unsigned int dimensions)
    {
        m_inputDimensions = dimensions;
    }

    auto gaussNumber()
    {
        return TSNE::gaussNumber();
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
    auto& operator<<(const T& value)
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

    auto createTempfile()
    {
        filestream.open(tempfile, std::ios::out | std::ios::binary | std::ios::trunc);
        writer = BinaryWriter(filestream);
    }

    auto removeTempfile()
    {
        if (filestream.is_open())
        {
            filestream.close();
        }
        EXPECT_EQ(0, remove(tempfile.c_str()));
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
    ASSERT_DOUBLE_EQ(1.1630780958763871, m_tsne.gaussNumber());
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
    ASSERT_DOUBLE_EQ(0.15606557998386178, m_tsne.gaussNumber());
    m_tsne.setRandomSeed(0);
    ASSERT_DOUBLE_EQ(1.1630780958763871, m_tsne.gaussNumber());
}

TEST_F(TsneTest, Perplexity)
{
    m_tsne.setPerplexity(5.0);
    EXPECT_EQ(5.0, m_tsne.perplexity());
    m_tsne.setPerplexity(25.0);
    EXPECT_EQ(25.0, m_tsne.perplexity());
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
    auto dataSize = 1;
    auto inputDimensions = 1;
    auto gradientAccuracy = 0.5;
    auto perplexity = 25.0;
    auto outputDimensions = 1;
    auto iterations = 100;
    auto randomSeed = 42;
    auto data = 42.0;

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
    ASSERT_DOUBLE_EQ(-0.51696416431811998, m_tsne.gaussNumber());
    EXPECT_EQ(data, m_tsne.data()[0][0]);

    removeTempfile();
}

TEST_F(TsneTest, LoadFromStream)
{
    auto data = std::vector<std::vector<double>>{ { 1.0, 2.0, 3.0 },{ 4.0, 5.0, 6.0 } };

    auto out = std::ostringstream();
    for (auto sample : data)
    {
        for (auto value : sample)
        {
            out << value << ',';
        }
        out.seekp(-1, std::ios_base::cur);
        out << std::endl;
    }

    auto in = std::istringstream(out.str());
    EXPECT_TRUE(m_tsne.loadFromStream(in));
    EXPECT_EQ(2, m_tsne.dataSize());
    EXPECT_EQ(3, m_tsne.inputDimensions());
    for (auto i = 0; i < 2; i++)
    {
        for (auto j = 0; j < 3; j++)
        {
            EXPECT_EQ(data[i][j], m_tsne.data()[i][j]);
        }
    }
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
    auto dataSize = 1;
    auto inputDimensions = 1;
    auto data = 42.0;

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
    auto data = std::vector<std::vector<double>>{
            { 0,56,19,80,58 },
            { 47,35,89,82,74 },
            { 17,85,71,51,30 },
            { 1,9,36,14,16 },
            { 98,44,11,0,0 },
            { 37,53,57,60,60 },
            { 16,66,45,35,5 },
            { 60,78,80,51,30 },
            { 87,72,95,92,53 },
            { 14,46,23,86,20 }
    };
    m_tsne.setDataSize(data.size());
    m_tsne.setInputDimensions(data[0].size());
    m_tsne.setGradientAccuracy(0.2);
    m_tsne.setPerplexity(3.0);
    m_tsne.setOutputDimensions(2);
    m_tsne.setIterations(100);
    m_tsne.setRandomSeed(42);
    m_tsne.setData(data);

    try {
        m_tsne.run();
    } catch (std::exception & e) {
        FAIL() << "run method exception: " << e.what();
    }

    auto result = m_tsne.result();

    auto expected = std::vector<std::vector<double>>{
        { -37.969682880670781, 13.53176098293458 },
        { 17.747120647448977, 11.769259684044014 },
        { 31.823026258301418, -9.5819351822503673 },
        { -9.325729922489689, -18.726748023156155 },
        { -20.281048368449309, -1.4635189684119658 },
        { 13.980657219661596, 15.40660294520521 },
        { -5.1885801679382269, -12.140850978402652 },
        { 34.539709042863528, -3.5143701898564852 },
        { 15.971308234165035, 17.035942809420174 },
        { -41.296780062892559, -12.316143079526352 }
    };

    for (size_t i = 0; i < m_tsne.dataSize(); i++)
    {
        for (size_t j = 0; j < m_tsne.outputDimensions(); j++)
        {
            EXPECT_TRUE(doubleEquals(expected[i][j], result[i][j]))
                << "expected != result at [" << i << "][" << j <<"]: "
                << expected[i][j] << " != " << result[i][j] << std::endl;
        }
    }
}

TEST_F(TsneTest, SaveLegacy)
{
    auto expected = std::vector<unsigned long long>{ 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    auto expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });

    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setOutputFile(tempfile);
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0] }, { expectedDouble[1] }});
    EXPECT_NO_THROW(m_tsne.saveLegacy());
    // check file exists and has right size
    auto result = std::ifstream();
    EXPECT_NO_THROW(result.open(tempfile + ".dat", std::ios::in | std::ios::binary | std::ios::ate));
    EXPECT_TRUE(result.is_open());
    EXPECT_EQ(48, result.tellg());
    // check values in file
    result.seekg(0);
    int dataSize;
    int outputDimensions;
    auto data = std::vector<double>(2);
    auto landmarks = std::vector<int>(2);
    auto costs = std::vector<double>(2);
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
    auto expected = std::vector<unsigned long long>{ 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    auto expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });

    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0] }, { expectedDouble[1] }});
    auto result = std::ostringstream();
    EXPECT_NO_THROW(m_tsne.saveToStream(result));
    // check values in stream
    auto expectedOut = std::ostringstream();
    expectedOut << std::setprecision(6);
    expectedOut << expectedDouble[0] << std::endl << expectedDouble[1] << std::endl;

    EXPECT_EQ(expectedOut.str(), result.str());
}

TEST_F(TsneTest, SaveCSV)
{
    auto expected = std::vector<unsigned long long>{ 0xC07DAC741DC680D4, 0x407DAC741DC680D4 };
    auto expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });

    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0] }, { expectedDouble[1] }});
    m_tsne.setOutputFile(tempfile);
    EXPECT_NO_THROW(m_tsne.saveCSV());
    // check file exists and has right size
    auto result = std::ifstream();
    EXPECT_NO_THROW(result.open(tempfile + ".csv", std::ios::in | std::ios::ate));
    EXPECT_TRUE(result.is_open());
    // file size may change based on current os (\n vs \r\n)
    EXPECT_LE(17, result.tellg());
    EXPECT_GE(19, result.tellg());
    // check values in file
    auto expectedOut = std::ostringstream();
    expectedOut << std::setprecision(6);
    expectedOut << expectedDouble[0] << std::endl << expectedDouble[1] << std::endl;
    auto expectedIn = std::istringstream(expectedOut.str());

    result.seekg(0);
    auto actualLine = std::string();
    auto expectedLine = std::string();
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
    auto expected = std::vector<unsigned long long>{ 0xC07DAC741DC680D4, 0xC07DAC741DC680D4, 0x407DAC741DC680D4, 0x407DAC741DC680D4 };
    auto expectedDouble = std::vector<double>();
    std::for_each(expected.begin(), expected.end(), [&](auto& v) { expectedDouble.push_back(*reinterpret_cast<double*>(&v)); });

    // test save
    m_tsne.setDataSize(2);
    m_tsne.setOutputDimensions(1);
    m_tsne.setResult(std::vector<std::vector<double>>{ { expectedDouble[0], expectedDouble[1] }, { expectedDouble[2], expectedDouble[3] }});
    m_tsne.setOutputFile(tempfile);
    EXPECT_NO_THROW(m_tsne.saveSVG());
    // check file exists
    auto result = std::ifstream();
    EXPECT_NO_THROW(result.open(tempfile + ".svg", std::ios::in));
    EXPECT_TRUE(result.is_open());
    // check values in file
    auto expectedRadius = 0.5;
    auto expectedViewBoxMin = expectedDouble[0] - expectedRadius;
    auto expectedViewBoxSize = 2 * -expectedViewBoxMin;
    auto expectedOut = std::ostringstream();
    expectedOut << std::setprecision(6);
    expectedOut << "<?xml version='1.0' encoding='UTF-8' ?>\n";
    expectedOut << "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='600' height='600' viewBox='" << expectedViewBoxMin << " " << expectedViewBoxMin << " " << expectedViewBoxSize << " " << expectedViewBoxSize << "'>\n";
    expectedOut << "<circle cx='" << expectedDouble[0] << "' cy='" << expectedDouble[1] << "' fill='black' r='" << expectedRadius << "' stroke='none' opacity='0.5'/>\n";
    expectedOut << "<circle cx='" << expectedDouble[2] << "' cy='" << expectedDouble[3] << "' fill='black' r='" << expectedRadius << "' stroke='none' opacity='0.5'/>\n";
    expectedOut << "</svg>\n";
    auto expectedIn = std::istringstream(expectedOut.str());

    auto actualLine = std::string();
    auto expectedLine = std::string();
    while ((bool)std::getline(result, actualLine) & (bool)std::getline(expectedIn, expectedLine))
    {
        EXPECT_EQ(expectedLine, actualLine);
    }

    EXPECT_TRUE(result.eof());
    EXPECT_TRUE(expectedIn.eof());

    result.close();
    EXPECT_EQ(0, remove((tempfile + ".svg").c_str()));
}

TEST_F(TsneTest, ResultConsistency)
{
    std::string input =
R"(0,56,19,80,58
47,35,89,82,74
17,85,71,51,30
1,9,36,14,16
98,44,11,0,0
37,53,57,60,60
16,66,45,35,5
60,78,80,51,30
87,72,95,92,53
14,46,23,86,20)";

    std::string expected =
R"(3.16374,-91.2264
81.7362,0.506676
-9.78514,23.3906
-81.2272,-9.16138
-73.7557,63.4705
51.512,-12.2191
-46.786,18.4669
23.8261,40.8484
77.3819,48.292
-26.0659,-82.3683
)";

    m_tsne.setGradientAccuracy(0.2);
    m_tsne.setPerplexity(3.0);
    m_tsne.setIterations(2000);
    m_tsne.setOutputDimensions(2);
    m_tsne.setRandomSeed(1337);

    auto iss = std::istringstream(input);
    EXPECT_TRUE(m_tsne.loadFromStream(iss));

    EXPECT_NO_THROW(m_tsne.run());

    auto oss = std::ostringstream();
    EXPECT_NO_THROW(m_tsne.saveToStream(oss));

    EXPECT_EQ(expected, oss.str());
}

TEST_F(TsneTest, ExactRun)
{
    std::string input =
R"(0,56,19,80,58
47,35,89,82,74
17,85,71,51,30
1,9,36,14,16
98,44,11,0,0
37,53,57,60,60
16,66,45,35,5
60,78,80,51,30
87,72,95,92,53
14,46,23,86,20)";

    std::string expected =
R"(-54.8785,-96.3901
-87.6806,46.6628
23.7272,19.651
82.25,-55.6907
115.188,26.6933
-62.33,15.9121
60.7718,-6.52401
-2.59574,57.4055
-56.0326,95.6331
-18.4195,-103.353
)";

    m_tsne.setGradientAccuracy(0.0);
    m_tsne.setPerplexity(3.0);
    m_tsne.setIterations(2000);
    m_tsne.setOutputDimensions(2);
    m_tsne.setRandomSeed(1337);

    auto iss = std::istringstream(input);
    EXPECT_TRUE(m_tsne.loadFromStream(iss));

    EXPECT_NO_THROW(m_tsne.runExact());

    auto oss = std::ostringstream();
    EXPECT_NO_THROW(m_tsne.saveToStream(oss));

    EXPECT_EQ(expected, oss.str());
}
