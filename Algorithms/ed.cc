#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <cstring>
#include <unistd.h>
#include "edlib.h"

using namespace std;

void usage(void)
{
    fprintf(stderr, "ed\n");
    fprintf(stderr, "  -x <seq_x_filename> input sequence file.\n");
    fprintf(stderr, "  -y <seq_y_filename> input sequence file.\n");
    fprintf(stderr, "  -o <output_filename> output file\n");
    exit(0);
}


unsigned char* readFile(const char * filename, std::streamsize& size)
{
    ifstream file(filename, ios::binary);

    if (!file)
    {
        cerr << "Failed to open file " << filename << endl;
        size = 0;
        return nullptr;
    }

    vector<unsigned char> temp;

    char ch;

    while (file.get(ch))
    {
        if (ch == ' ' || ch == '\n' || ch == '\r' || ch == '\t')
            continue;
        
        temp.push_back(static_cast<unsigned char>(ch));
    }

    size = temp.size();

    unsigned char * sequence = new unsigned char[size+1];

    for (int i = 0; i < size; i++)
    {
        sequence[i] = temp[i];
    }
    
    sequence[size] = '\0';
    temp.clear();
    
    return sequence;
}

int main(int argc, char** argv)
{
    string x;
    string y;
    string output;
    
    int opt;

    while ((opt = getopt(argc, argv, "x:y:o:")) != -1)
    {
        switch (opt)
        {
            case 'x':
                x = optarg;
                break;

            case 'y':
                y = optarg;
                break;

            case 'o':
                output = optarg;
                break;

            default:
                usage();
        }
    }


    if (output.empty())
    {
        fprintf(stderr, "Error: output file (-o) is required.\n");
        usage();
    }

    auto start = chrono::high_resolution_clock::now();

    streamsize size_x;
    unsigned char * seq_x = readFile(x.c_str(), size_x);

    streamsize size_y;
    unsigned char * seq_y = readFile(y.c_str(), size_y);

    int ed = edlibAlign( (const char*) seq_x, strlen( (char*) seq_x ), (const char*) seq_y, strlen( (char*) seq_y ), edlibDefaultAlignConfig()).editDistance;
    
    ofstream out(output);
    out <<ed<< endl;
   //cout<<ed<<endl;
    
    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout<< "Elapsed time: " << elapsed.count() << " ms\n";

    delete[] seq_x;
    delete[] seq_y;
   
    return 0;
}
