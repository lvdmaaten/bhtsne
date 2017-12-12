
#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include <bhtsne/tsne.h>


struct MeasurementResult
{
    using rep = std::chrono::high_resolution_clock::rep;

    rep preparation_time;
    rep execution_time;
    rep save_time;
};


bool fileIsPresent()
{
    auto file = std::ifstream("data.dat");
    return file.good();
}

int main(int argc, char* argv[])
{
    auto testSizes = std::vector<int>{ 250, 500, 750, 1000 };
    auto iteration_times = std::vector<int>{ 250, 500, 750, 1000 };

    std::vector<std::vector<MeasurementResult>> runtimes;
    runtimes.resize(testSizes.size());
    for (auto& each : runtimes)
        each.resize(iteration_times.size());

    for (auto i = 0; i < testSizes.size(); ++i)
    {
        for (auto j = 0; j < iteration_times.size(); ++j)
        {
            auto testSize = testSizes[i];
            auto iterations = iteration_times[j];

            //TODO update to new interface
            /*
            auto start_prepare = std::chrono::high_resolution_clock::now();

            int input_dimension = 784; //magic
            int output_dimension = 2;
            double perplexity = 50;
            double gradient_accuracy = 0.5;
            double* data;
            int rand_seed = 0;
            bhtsne::TSNE tsne;
            if (!fileIsPresent())
            {
                std::cout << "data.dat not found" << std::endl;
            }
            int temp;
            double temp2;
            tsne.load_data(&data, &temp, &temp, &temp2, &temp2, &temp, &temp);
            tsne.setDataSize(testSize);
            int* landmarks = (int*)malloc(testSize * sizeof(int));
            if (landmarks == NULL) { printf("Memory allocation failed!\n"); exit(1); }
            for (int n = 0; n < testSize; n++) landmarks[n] = n;
            double* result = (double*)malloc(testSize * output_dimension * sizeof(double));
            double* costs = (double*)calloc(testSize, sizeof(double));
            if (result == NULL || costs == NULL) { printf("Memory allocation failed!\n"); exit(1); }

            auto start_execute = std::chrono::high_resolution_clock::now();

            tsne.run(data, input_dimension, result, output_dimension, perplexity, gradient_accuracy, rand_seed, false, iterations);

            auto end_execute = std::chrono::high_resolution_clock::now();

            tsne.save_data(result, landmarks, costs, testSize, output_dimension);
            free(data);
            free(result);
            free(costs);
            free(landmarks);

            auto end_save = std::chrono::high_resolution_clock::now();

            auto& current_result = runtimes[i][j];
            current_result.preparation_time = (start_execute - start_prepare).count();
            current_result.execution_time = (end_execute - start_execute).count();
            current_result.save_time = (end_save - end_execute).count();

            */
        }
    }

    // write header
    std::cout << "testsize;iterations;preparation_time;execution_time;save_time" << std::endl;

    // write contents
    for (auto i = 0; i < testSizes.size(); ++i)
    {
        for (auto j = 0; j < iteration_times.size(); ++j)
        {
            auto testSize = testSizes[i];
            auto iterations = iteration_times[j];

            const auto& result = runtimes[i][j];
            std::cout
                << testSize << ";"
                << iterations << ";"
                << result.preparation_time << ";"
                << result.execution_time << ";"
                << result.save_time << std::endl;
        }
    }
}
