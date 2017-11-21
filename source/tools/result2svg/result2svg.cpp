#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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

    uint32_t N;
    input.read((char*)&N, 4);
    cout << "Result contains " << to_string(N) << " samples." << endl;

    uint32_t no_dims;
    input.read((char*)&no_dims, 4);
    cout << "Number of dimensions is " << to_string(no_dims) << "." << endl;

    cout << "Set includes " << N << " samples with " << no_dims << " values each." << endl;
    cout << "Total number of values is " << (N * no_dims) << "." << endl;
    
    uint8_t* data = new uint8_t[N * no_dims * sizeof(double)];
    input.read((char*)data, N * no_dims * sizeof(double));
    cout << "Loaded data." << endl;

    input.close();
    cout << "Closed input." << endl;

    double* data_d = new double[N * no_dims];
    double max = 0;
    for(int i = 0; i < N * no_dims; i++)
    {
        double d = ((double*)data)[i];
        data_d[i] = d;
        double a = d < 0 ? -d : d;
        if(a > max) max = a;
    }
    double radius = 0.1;
    string viewBox = to_string(-max-radius) + " " + to_string(-max-radius) + " " + to_string(2*max+2*radius) + " " + to_string(2*max+2*radius);

    uint8_t* labels;
    bool useLabels = false;
    if (argc > 2)
    {
        useLabels = true;
        ifstream input_l;
        input_l.open(argv[2], ios::in | ios::binary);
        if (!input_l.is_open()) fail("Could not open labels.");

        uint32_t N_l;
        input_l.read((char*)&N_l, 4);
        cout << "Labels file contains " << to_string(N_l) << " labels." << endl;
        if (N_l < N) fail("Not enough labels for result.");

        labels = new uint8_t[N_l];
        input_l.read((char*)labels, N_l);

        input_l.close();
        cout << "Read labels." << endl;
    }

    auto colors = vector<string>{ "red", "blue", "green", "yellow", "magenta", "brown", "orange", "lime", "lightblue", "pink" };

    ofstream output;
    output.open("result.svg", ios::out | ios::trunc);
    cout << "Opened output." << endl;

    output << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
    output << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"600\" height=\"600\" viewBox=\"" << viewBox << "\">" << endl;
    cout << "Wrote header." << endl;
    
    string color = "black";
    for(int i = 0; i < N; i++)
    {
        double d0 = ((double*)data)[i*2];
        double d1 = ((double*)data)[i*2+1];
        if (useLabels)
            color = labels[i] < colors.size() ? colors[labels[i]] : "black";
        output << "<circle cx=\"" << d0 << "\" cy=\"" << d1 << "\" r=\"0.1\" fill=\"" << color << "\" stroke=\"" << color << "\"/>" << endl;
    }
    output << "</svg>" << endl;
    cout << "Wrote data." << endl;

    output.close();
    cout << "Closed output." << endl;

    return 0;
}