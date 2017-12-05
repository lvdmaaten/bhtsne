#include <iostream>

#include "ArgumentParser.h"
#include "CommandlineOptions.h"
#include <bhtsne/tsne.h>


using namespace std;


// Function that runs the Barnes-Hut implementation of t-SNE
int main(int argc, char * argv[])
{
    auto tsne = bhtsne::TSNE();

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
        cerr << "only one input file at a time" << endl;
        return 1;
    }

    if (params.empty())
    {
        //TODO
        //read from stdin
    }
    else
    {
        //determine extension
        auto inputfile = params.front();
        auto ext_pos = inputfile.rfind('.');
        if (ext_pos == std::string::npos)
        {
            cerr << "input file needs an extension" << endl;
            return 1;
        }
        auto ext = inputfile.substr(ext_pos);

        //load correct file
        bool loaded = false;
        if (ext == ".dat")
            loaded = tsne.loadLegacy(inputfile);
        else if (ext == ".csv")
            loaded = tsne.loadCSV(inputfile);
        else if (ext == ".tsne")
            loaded = tsne.loadTSNE(inputfile);
        else
        {
            cerr << "file extension of " << ext << " not supported" << endl;
            return 1;
        }

        if (!loaded)
        {
            cerr << "failed to load " << inputfile << endl;
            return 1;
        }
    }    

	tsne.run();

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
