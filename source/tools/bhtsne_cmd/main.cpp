#include <iostream>

#include "ArgumentParser.h"
#include "CommandlineOptions.h"
#include <bhtsne/TSNE.h>


// Function that runs the Barnes-Hut implementation of t-SNE
int main(int argc, char * argv[])
{
    auto tsne = bhtsne::TSNE();

    //parse and apply command line arguments
	auto parsedArguments = cppassist::ArgumentParser();
	parsedArguments.parse(argc, argv);

    //read correct input file
    auto params = parsedArguments.params();
    if (params.size() > 1)
    {
        std::cerr << "only one input file at a time\n";
        return 2;
    }

    auto loaded = false;

    if (params.empty())
    {
        loaded = tsne.loadCin();
    }
    else
    {
        //determine extension
        auto inputfile = params.front();
        auto ext_pos = inputfile.rfind('.');
        if (ext_pos == std::string::npos)
        {
            std::cerr << "input file needs an extension\n";
            return 3;
        }
        auto ext = inputfile.substr(ext_pos);

        //load correct file
        
        if (ext == ".dat")
        {
            loaded = tsne.loadLegacy(inputfile);
        }
        else if (ext == ".csv")
        {
            loaded = tsne.loadCSV(inputfile);
        }
        else if (ext == ".tsne")
        {
            loaded = tsne.loadTSNE(inputfile);
        }
        else
        {
            std::cerr << "file extension of " << ext << " not supported\n";
            return 4;
        }
    }
    if (!loaded)
    {
        std::cerr << "failed to load input file\n";
        return 5;
    }

    applyCommandlineOptions(tsne, parsedArguments.options());

	tsne.run();

    //save in requested formats
	if (parsedArguments.isSet("-legacy"))
	{
		tsne.saveLegacy();
	}
	if (parsedArguments.isSet("-svg"))
	{
		tsne.saveSVG();
	}
	if (parsedArguments.isSet("-csv"))
	{
		tsne.saveCSV();
	}

    if (parsedArguments.isSet("-stdout"))
    {
        tsne.saveToCout();
    }
}
