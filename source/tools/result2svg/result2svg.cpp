#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>

using namespace std;

void fail(string reason)
{
    cout << "Error: " << reason << endl;
    exit(-1);
}

int main(int argc, char* argv[])
{
    string file = argc >= 2 ? argv[1] : "result.dat";
    ifstream input;
    input.open(file, ios::in | ios::binary);
    if(!input.is_open()) fail("Could not open input.");
    cout << "Loaded input." << endl;

    uint32_t sampleCount, dimensionCount;
    input.read(reinterpret_cast<char*>(&sampleCount), sizeof(sampleCount));
    input.read(reinterpret_cast<char*>(&dimensionCount), sizeof(dimensionCount));

    cout << "Set includes " << sampleCount << " samples with " << dimensionCount << " values each." << endl;
    cout << "Total number of values is " << (sampleCount * dimensionCount) << "." << endl;
    
    auto data = vector<double>(sampleCount * dimensionCount);
    input.read(reinterpret_cast<char*>(data.data()), data.size() * sizeof(double));
    cout << "Loaded data." << endl;

    input.close();
    cout << "Closed input." << endl;

    double extreme = 0;
    for(int i = 0; i < sampleCount * dimensionCount; i++)
        extreme = max(extreme, abs(data[i]));
    double radius = 0.5;
    double halfwidth = extreme + radius;
    string viewBox = to_string(-halfwidth) + " " + to_string(-halfwidth) + " " + to_string(2*halfwidth) + " " + to_string(2*halfwidth);

    auto labels = vector<uint8_t>();
    bool useLabels = false;
    if (argc > 2)
    {
        useLabels = true;
        ifstream labelInput;
        labelInput.open(argv[2], ios::in | ios::binary);
        if (!labelInput.is_open()) fail("Could not open labels.");

        uint32_t labelCount;
        labelInput.read(reinterpret_cast<char*>(&labelCount), sizeof(labelCount));
        cout << "Labels file contains " << labelCount << " labels." << endl;
        if (labelCount < sampleCount) fail("Not enough labels for result.");

        labelCount = min(labelCount, sampleCount);
        labels.resize(labelCount);
        labelInput.read(reinterpret_cast<char*>(labels.data()), labels.size());

        labelInput.close();
        cout << "Read labels." << endl;
    }

    uint8_t maxLabel = 0;
    for (auto label : labels)
        maxLabel = max(label, maxLabel);
    auto colors = vector<string>();
    for (int i = 0; i <= maxLabel; ++i)
        colors.push_back("hsl(" + to_string(360.0 * i / (maxLabel / 2 + 1)) + ", 100%, " + (i % 2 == 0 ? "25" : "60") + "%)");

    ofstream output;
    output.open("result.svg", ios::out | ios::trunc);
    cout << "Opened output." << endl;

    output << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
    output << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"600\" height=\"600\" viewBox=\"" << viewBox << "\">" << endl;
    cout << "Wrote header." << endl;
    
    string color = "black";
    for(int i = 0; i < sampleCount; i++)
    {
        if (useLabels)
            color = labels[i] < colors.size() ? colors[labels[i]] : "black";

        output << "<circle "
            << "cx='" << data[i * 2] << "' "
            << "cy='" << data[i * 2 + 1] << "' "
            << "fill='" << color << "' "
            << "r='" << radius << "' "
            << "stroke='none' opacity='0.5'/>" << endl;
    }
    output << "</svg>" << endl;
    cout << "Wrote data." << endl;

    output.close();
    cout << "Closed output." << endl;

    return 0;
}