#include <iostream>
#include <fstream>
#include <string>

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

    ofstream output;
    output.open("result.svg", ios::out | ios::trunc);
    cout << "Opened output." << endl;

    output << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
    output << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"600\" height=\"600\" viewBox=\"" << viewBox << "\">" << endl;
    cout << "Wrote header." << endl;
    
    for(int i = 0; i < N * no_dims; i+=2)
    {
        double d0 = ((double*)data)[i];
        double d1 = ((double*)data)[i+1];
        output << "<circle cx=\"" << d0 << "\" cy=\"" << d1 << "\" r=\"0.1\" fill=\"black\" />" << endl;
    }
    output << "</svg>" << endl;
    cout << "Wrote data." << endl;

    output.close();
    cout << "Closed output." << endl;

    return 0;
}