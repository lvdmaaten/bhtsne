#include "CommandlineOptions.h"
#include "ArgumentParser.h"

namespace bhtsne {
    void applyCommandlineOptions(TSNE &tsne, int argc, char *argv[]) {
        auto parser = cppassist::ArgumentParser();
        parser.parse(argc, argv);

        //set parameter values
        for (auto optionValuePair : parser.options()) {
            if (optionValuePair.first == "--perplexity") {
                tsne.setPerplexity(std::stod(optionValuePair.second));
            } else if (optionValuePair.first == "--gradient-accuracy") {
                tsne.setGradientAccuracy(std::stod(optionValuePair.second));
            } else if (optionValuePair.first == "--iterations") {
                tsne.setIterations(static_cast<unsigned int>(std::stol(optionValuePair.second)));
            } else if (optionValuePair.first == "--number-of-samples") {
                tsne.setNumberOfSamples(static_cast<unsigned int>(std::stol(optionValuePair.second)));
            } else if (optionValuePair.first == "--output-dimensions") {
                tsne.setOutputDimensions(static_cast<unsigned int>(std::stol(optionValuePair.second)));
            } else if (optionValuePair.first == "--output-file") {
                tsne.setOutputFile(optionValuePair.second);
            } else if (optionValuePair.first == "--random-seed") {
                tsne.setRandomSeed(std::stoi(optionValuePair.second));
            }
        }

        //TODO bool flags

    }
}// namespace bhtsne