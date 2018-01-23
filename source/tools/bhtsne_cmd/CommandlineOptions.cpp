#include <iostream>

#include "CommandlineOptions.h"

namespace bhtsne
{
    void applyCommandlineOptions(TSNE & tsne, const std::map<std::string, std::string> & options)
    {
        //set parameter values
        for (const auto & optionValuePair : options)
        {
            if (optionValuePair.first == "--perplexity")
            {
                tsne.setPerplexity(std::stod(optionValuePair.second));
            }
            else if (optionValuePair.first == "--gradient-accuracy")
            {
                tsne.setGradientAccuracy(std::stod(optionValuePair.second));
            }
            else if (optionValuePair.first == "--iterations")
            {
                tsne.setIterations(static_cast<unsigned int>(std::stol(optionValuePair.second)));
            }
            else if (optionValuePair.first == "--output-dimensions")
            {
                tsne.setOutputDimensions(static_cast<unsigned int>(std::stol(optionValuePair.second)));
            }
            else if (optionValuePair.first == "--output-file")
            {
                tsne.setOutputFile(optionValuePair.second);
            }
            else if (optionValuePair.first == "--random-seed")
            {
                tsne.setRandomSeed(std::stoul(optionValuePair.second));
            }
            else if (optionValuePair.first.find("--") == 0)
            {
                std::cerr << "warning: ignored unexpected command line option " << optionValuePair.first << "\n"
                    << "allowed options are: --perplexity, --gradient-accuracy, --iterations, "
                    << "--output-dimensions, --output-file, --random-seed\n";
            }
        }
    }

}// namespace bhtsne
