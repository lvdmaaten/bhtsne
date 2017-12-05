#include <iostream>

#include "ArgumentParser.h"
#include "CommandlineOptions.h"
#include <bhtsne/tsne.h>


// Function that runs the Barnes-Hut implementation of t-SNE
int main(int argc, char * argv[])
{
    auto tsne = bhtsne::TSNE();

    //parse and apply command line arguments
	auto parsedArguments = cppassist::ArgumentParser();
	parsedArguments.parse(argc, argv);

	if (!parsedArguments.isSet("-legacy")
		&& !parsedArguments.isSet("-svg")
		&& !parsedArguments.isSet("-csv"))
	{
		std::cerr << "please specify one of the output options: -csv, -svg, -legacy\n";
		return 1;
	}

    applyCommandlineOptions(tsne, parsedArguments.options());

    //read correct input file
    auto params = parsedArguments.params();
    if (params.size() > 1)
    {
        std::cerr << "only one input file at a time\n";
        return 2;
    }

    if (params.empty())
    {
        //TODO(!): read from stdin
        std::cerr << "reading from stdin is not supported yet, please specify an input file\n";
        return 42;
    } else {
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
        bool loaded = false;
        if (ext == ".dat")
        {
            loaded = tsne.loadLegacy(inputfile);
        } else if (ext == ".csv")
        {
            loaded = tsne.loadCSV(inputfile);
        } else if (ext == ".tsne")
        {
            loaded = tsne.loadTSNE(inputfile);
        } else {
            std::cerr << "file extension of " << ext << " not supported\n";
            return 4;
        }

        if (!loaded)
        {
            std::cerr << "failed to load " << inputfile << "\n";
            return 5;
        }
    }

	tsne.run();

    //save in requested formats
	if (parsedArguments.isSet("-legacy"))
	{
		tsne.saveLegacy();
	}
	if (parsedArguments.isSet("-svg"))
	{
		tsne.saveCSV();
	}
	if (parsedArguments.isSet("-csv"))
	{
		tsne.saveSVG();
	}
}
