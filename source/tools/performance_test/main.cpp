#include <vector>
#include <chrono>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>

#include <bhtsne/TSNE.h>


struct MeasurementResult
{
    using rep = std::chrono::high_resolution_clock::rep;

    rep preparation_time;
    rep execution_time;
    rep save_time;
};

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Please specify the warmup and test iterations. Example: performance_test.exe 3 5" << std::endl;
        return 0;
    }

    const auto testSizes = std::vector<int>{ 250, 500, 750, 1000 };
    const auto iterationTimes = std::vector<int>{ 250, 500, 750, 1000 };
    const auto warmupIterations = std::atoi(argv[1]);
    const auto testIterations = std::atoi(argv[2]);

    auto runtimes = std::vector<std::vector<MeasurementResult>>(testSizes.size());
    for (auto & each : runtimes)
        each.resize(iterationTimes.size());

    // hide library output
    auto coutBuf = std::cout.rdbuf();
    std::ostream newCout(coutBuf);
    std::stringstream redirectStream;
    std::cout.rdbuf(redirectStream.rdbuf());

    for (size_t i = 0; i < testSizes.size(); ++i)
    {
        for (size_t j = 0; j < iterationTimes.size(); ++j)
        {
            auto testSize = testSizes[i];
            auto iterations = iterationTimes[j];

            auto & current_result = runtimes[i][j];
            current_result.preparation_time = 0;
            current_result.execution_time = 0;
            current_result.save_time = 0;

            for (auto k = 0; k < testIterations + warmupIterations; ++k)
            {
                newCout << "Running " << (k < warmupIterations ? "warmup" : "test") << " iteration "
                        << (k + 1) << "/" << (testIterations + warmupIterations)
                        << " for " << testSize << " samples and " << iterations << " iterations." << std::endl;
                auto start_prepare = std::chrono::high_resolution_clock::now();

                auto tsne = bhtsne::TSNE();

                tsne.loadLegacy("./data/data_s" + std::to_string(testSize) + "_i1000.dat");

                tsne.setOutputDimensions(2);
                tsne.setPerplexity(50);
                tsne.setGradientAccuracy(0.2);
                tsne.setRandomSeed(0);
                tsne.setOutputFile("./result");
                tsne.setIterations(iterations);

                auto start_execute = std::chrono::high_resolution_clock::now();

                tsne.run();

                auto end_execute = std::chrono::high_resolution_clock::now();

                tsne.saveLegacy();

                auto end_save = std::chrono::high_resolution_clock::now();

                if (k >= warmupIterations)
                {
                    current_result.preparation_time += (start_execute - start_prepare).count();
                    current_result.execution_time += (end_execute - start_execute).count();
                    current_result.save_time += (end_save - end_execute).count();
                }

                remove("result.dat");
            }

            current_result.preparation_time /= testIterations;
            current_result.execution_time /= testIterations;
            current_result.save_time /= testIterations;
        }
    }

    newCout << "All tests done." << std::endl;

    auto saveTime = std::chrono::system_clock::now();
    auto saveTime_t = std::chrono::system_clock::to_time_t(saveTime);
    char* timeString = new char[20];
    std::strftime(timeString, 20, "%Y-%m-%d-%H-%M-%S", std::localtime(&saveTime_t));
    auto fileName = "performance_" + std::string(timeString) + (argc > 3 ? "_" + std::string(argv[3]) : "") + ".csv";
    std::ofstream file(fileName);
    if (!file.is_open())
    {
        newCout << "Could not open output file " << fileName << "." << std::endl;
    }

    // write header
    file << "testsize;iterations;preparation_time;execution_time;save_time" << std::endl;

    // write contents
    for (size_t i = 0; i < testSizes.size(); ++i)
    {
        for (size_t j = 0; j < iterationTimes.size(); ++j)
        {
            auto testSize = testSizes[i];
            auto iterations = iterationTimes[j];

            const auto & result = runtimes[i][j];
            file
                << testSize << ";"
                << iterations << ";"
                << result.preparation_time << ";"
                << result.execution_time << ";"
                << result.save_time << std::endl;
        }
    }

    file.close();
}
