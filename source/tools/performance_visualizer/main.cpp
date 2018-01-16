
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <functional>

struct PerformanceData
{
	std::string commitName;
	int testsize;
	int iterations;
	long long preparation_time;
	long long execution_time;
	long long save_time;
};

const auto svgFileName = "performance.svg";

const auto width = 600.0;
const auto height = 400.0;

const auto paddingLeft = 45.0;
const auto paddingRight = 180.0;
const auto paddingTop = 20.0;
const auto paddingBottom = 80.0;

const auto timeStep = 10.0;
const auto timeSubStep = 1.0;

const auto preperationColor = std::string("#f00");
const auto executionColor   = std::string("#0f0");
const auto saveColor        = std::string("#00f");

auto maxTime = 0.0;
auto xStep = 100.0;

std::ofstream svg;
auto m_data = std::vector<PerformanceData>();
int iterationsToPlot;
int testsizeToPlot;

void createSvgHeader()
{
	svg << "<?xml version='1.0' encoding='UTF-8' ?>\n"
		<< "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' "
		<< "width='" << width + paddingLeft + paddingRight << "' height='" << height + paddingTop + paddingBottom << "' "
		<< "viewBox='" << -paddingLeft << " " << -paddingTop << " " << width + paddingRight + paddingLeft << " " << height + paddingBottom + paddingTop << "'>\n"
		<< "<namedview scale-x='1'/>";
}

void drawHorizontalLines()
{
	for (auto time = 0.0; time <= maxTime; time += timeSubStep) 
	{
		auto y = height - (time * height / maxTime);
		svg << "<path style='stroke:#eee;stroke-width:1px;' "
			<< "d='M 0," << y << " H " << width << "'"
			<< "/>\n";
	}
	for (auto time = 0.0; time <= maxTime; time += timeStep)
	{
		auto y = height - (time * height / maxTime);
		svg << "<path style='stroke:#ccc;stroke-width:1px;' "
			<< "d='M 0," << y << " H " << width << "'"
			<< "/>\n";

		svg << "<path style='stroke:#000;stroke-width:1px;' "
			<< "d='M 0," << y << " H " << -5 << "'"
			<< "/>\n";

		svg << "<text style='fill:#000;font-family:sans-serif;font-size:20px;text-anchor:end;text-align:end;' "
			<< "x='" << -10 << "' "
			<< "y='" << y + 5 << "'"
			<< "><tspan>" << time << "</tspan></text>\n";
	}
}

void drawVerticalLines()
{
	auto pos = 0;
	for (auto data : m_data)
	{
		auto x = xStep * pos++;
		svg << "<path style='stroke:#ccc;stroke-width:1px;' "
			<< "d='M " << x << ",0 V " << height << "'"
			<< "/>\n";

		svg << "<path style='stroke:#000;stroke-width:1px;' "
			<< "d='M " << x << "," << height << " V " << height + 5 << "'"
			<< "/>\n";

		svg << "<text style='fill:#000;font-family:sans-serif;font-size:20px;text-anchor:start;text-align:start;' transform='rotate(45 " << x - 10 << " " << height + 5 << ")' "
			<< "x='" << x << "' "
			<< "y='" << height + 5 << "'"
			<< "><tspan>" << data.commitName.substr(0, 6) << "</tspan></text>\n";
	}
}

template<typename F>
void drawGraph(const std::string &color, const F &accaccessor,const double scaleFactor = 1.0)
{
	auto circles = std::stringstream();
	auto line = std::stringstream();
	line << "<path style='stroke:" << color << ";stroke-width:1px;fill:none;' "
		<< "d='M";
	auto pos = 0;
	for (auto &data : m_data)
	{
		auto x = xStep * pos;
		auto y = height - accaccessor(data) * 1.0E-9 / maxTime * height * scaleFactor;

		line << " " << x << "," << y;

		circles << "<circle style='fill:" << color << "' r='3' "
			<< "cx='" << x << "' "
			<< "cy='" << y << "'"
			<< "/>\n";
		pos++;
	}

	line << "'/>\n";

	svg << line.str();
	svg << circles.str();
}

void drawKey(const std::string &label, const std::string &color, const double yOffset)
{
	auto x = width + 100;
	auto y = height / 2 + yOffset;
	svg << "<text style='fill:#000;font-family:sans-serif;font-size:20px;text-anchor:middle;text-align:center;' "
		<< "x='" << x << "' "
		<< "y='" << y << "'"
		<< "><tspan>" << label << "</tspan></text>\n";

	svg << "<path style='stroke:" << color << ";stroke-width:1px;' "
		<< "d='M " << x - 10 << "," << y + 10 << " H " << x + 10 << "'"
		<< "/>\n";

	svg << "<circle style='fill:" << color << "' r='3' "
		<< "cx='" << x << "' "
		<< "cy='" << y + 10 << "'"
		<< "/>\n";
}

int main(int argc, char* argv[])
{
	if (argc < 4) {
		std::cerr << "Not enough arguments.\n" << argv[0] << " iterationToPlot testsizeToPlot filename+";
		return 1;
	}
	iterationsToPlot = atoi(argv[1]);
	testsizeToPlot = atoi(argv[2]);
	//load data
	for (auto i = size_t(3); i < argc; ++i)
	{
		auto fileName = std::string( argv[i] );

		auto posOfLastDot = fileName.find_last_of('.');
		auto posOfLastUnderscore = fileName.find_last_of('_');
		auto commitName = fileName.substr(posOfLastUnderscore + 1, posOfLastDot - posOfLastUnderscore - 1);

		auto performanceFile = std::ifstream();
		performanceFile.open(fileName);
		
		if (!performanceFile.is_open())
		{
			std::cerr << "Could not open " << fileName << std::endl;
		}

		auto line = std::string();
		bool first = true;
		while (std::getline(performanceFile, line))
		{
			if (first) {
				first = false;
				continue;
			}
			auto iss = std::istringstream(line);
			auto element = std::string();

			auto values = PerformanceData();
			values.commitName = commitName;

			std::getline(iss, element, ';');
			values.testsize = std::stoi(element);

			if (values.testsize != testsizeToPlot) continue;

			std::getline(iss, element, ';');
			values.iterations = std::stoi(element);

			if (values.iterations != iterationsToPlot) continue;

			std::getline(iss, element, ';');
			values.preparation_time = std::stoll(element);
			maxTime = std::fmax(maxTime, std::ceil(static_cast<double>(values.preparation_time) / 1'000'000));

			std::getline(iss, element, ';');
			values.execution_time = std::stoll(element);
			maxTime = std::fmax(maxTime, std::ceil(static_cast<double>(values.execution_time) / 1'000'000'000));

			std::getline(iss, element, '\n');
			values.save_time = std::stoll(element);
			maxTime = std::fmax(maxTime, std::ceil(static_cast<double>(values.save_time) / 1'000'000));

			m_data.emplace_back(std::move(values));
		}

	}

	std::cout << "maxTime: " << maxTime << std::endl;

	//save svg
	svg = std::ofstream();
	svg.open(svgFileName);

	if (!svg.is_open())
	{
		std::cerr << "Could not open " << svgFileName << std::endl;
		return 1;
	}

	xStep = width / (m_data.size() - 1);

	createSvgHeader();
	drawHorizontalLines();
	drawVerticalLines();
	drawGraph(preperationColor, [](PerformanceData d) {return d.preparation_time; }, 1000.0);
	drawGraph(executionColor, [](PerformanceData d) {return d.execution_time; });
	drawGraph(saveColor, [](PerformanceData d) {return d.save_time; }, 1000.0);

	drawKey(std::string("preparation in ms"), preperationColor, -50.0);
	drawKey(std::string("execution in s"), executionColor, 0.0);
	drawKey(std::string("save in ms"), saveColor, 50.0);

	svg << "</svg>";

	svg.close();
}
