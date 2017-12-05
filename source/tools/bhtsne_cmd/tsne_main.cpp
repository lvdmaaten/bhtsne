#include <iostream>

#include "ArgumentParser.h"
#include "CommandlineOptions.h"
#include <bhtsne/tsne.h>

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
