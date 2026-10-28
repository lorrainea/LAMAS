#include <stdio.h>
#include <set>
#include <cstring>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>
#include <climits>
#include <vector>
#include <map>

using namespace std;
void usage(void)
{
    fprintf(stderr, "qgram\n");
    fprintf(stderr, "  -x <seq_x_filename> input sequence file.\n");
    fprintf(stderr, "  -y <seq_y_filename> input sequence file.\n");
    fprintf(stderr, "  -q <q_gram> q.\n");
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

    int q = 0;

    int opt;

    while ((opt = getopt(argc, argv, "x:y:q:o:")) != -1)
    {
        switch (opt)
        {
            case 'x':
                x = optarg;
                break;

            case 'y':
                y = optarg;
                break;

            case 'q':
            {
                stringstream ss(optarg);
                ss >> q;
                break;
            }

            case 'o':
                output = optarg;
                break;

            default:
                usage();
        }
    }

    if (q < 0)
    {
        fprintf(stderr, "Error: q-gram (-q) cannot be less than 0.\n");
        usage();
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

    int seq_x_len = strlen((char*) seq_x);
    int seq_y_len = strlen((char*) seq_y);
    
    
    // frequency maps for q-grams
    map<string, int> freq_x;
    map<string, int> freq_y;

    // extract q-grams from x
    for (int i = 0; i <= seq_x_len - q; i++)
    {
        string qgram((char*)seq_x + i, q);
        freq_x[qgram]++;
    }

    // extract q-grams from y
    for (int i = 0; i <= seq_y_len - q; i++)
    {
        string qgram((char*)seq_y + i, q);
        freq_y[qgram]++;
    }

    // compute q-gram distance
    int distance = 0;

    set<string> all_qgrams;

    for (auto& p : freq_x)
        all_qgrams.insert(p.first);

    for (auto& p : freq_y)
        all_qgrams.insert(p.first);

    for (const auto& qgram : all_qgrams)
    {
        int count_x = freq_x[qgram];
        int count_y = freq_y[qgram];

        distance += abs(count_x - count_y);
    }

    ofstream out(output);



    out << distance << endl;

    out.close();
        
  
    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout<< "Elapsed time: " << elapsed.count() << " ms\n";

    delete[] seq_x;
    delete[] seq_y;

    return 0;
}
