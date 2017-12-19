
#include <vector>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
	const auto svgFileName = "performance.svg";

	auto svg = std::ofstream();
	svg.open(svgFileName);

	if (!svg.is_open())
	{
		std::cerr << "Could not open " << svgFileName << std::endl;
		return 1;
	}

	auto width = 600.0;
	auto height = 400.0;

	svg << "<?xml version='1.0' encoding='UTF-8' ?>\n";
	svg << "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='" << width << "' height='" << height << "' ";
	svg << "viewBox='0 0 " << width << " " << height << "'>\n";


	svg << "</svg>";

	svg.close();
}
