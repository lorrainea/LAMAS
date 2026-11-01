#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>
#include <climits>
#include <cstring>
#include <vector>
#include <unistd.h>
#include <divsufsort64.h>

using namespace std;

void usage(void)
{
    fprintf(stderr, "lamas_linear\n");
    fprintf(stderr, "  -x <seq_x_filename> input sequence X file.\n");
    fprintf(stderr, "  -y <seq_y_filename> input sequence Y file.\n");
    fprintf(stderr, "  -o <output_filename> output file\n");
    exit(EXIT_FAILURE);
}

/* Kasai et al algorithm for O(n)-time LCP construction */
int64_t LCParray(unsigned char *text, int64_t n, int64_t *SA, int64_t *ISA, int64_t *LCP)
{
    int64_t i = 0, j = 0;

    LCP[0] = 0;
    for (i = 0; i < n; i++)
    {
        if (ISA[i] != 0)
        {
            if (i == 0)
            {
                j = 0;
            }
            else
            {
                j = (LCP[ISA[i - 1]] >= 2) ? LCP[ISA[i - 1]] - 1 : 0;
            }
            
            while (text[i + j] == text[SA[ISA[i] - 1] + j])
            {
                j++;
            }
            LCP[ISA[i]] = j;
        }
    }
    return 0;
}

unsigned char *concatenate(unsigned char *seq_x, unsigned char *seq_y, int seq_x_len, int seq_y_len)
{
    unsigned char *seq_xy = new unsigned char[seq_x_len + seq_y_len + 2];
    
    strcpy(reinterpret_cast<char*>(seq_xy), reinterpret_cast<char*>(seq_x));

    int pos = seq_x_len;
    seq_xy[pos] = '~';
    seq_xy[pos + 1] = '\0';
    
    strcat(reinterpret_cast<char*>(seq_xy), reinterpret_cast<char*>(seq_y));
    seq_xy[seq_x_len + seq_y_len + 1] = '\0';
    
    return seq_xy;
}

unsigned char *readFile(const char *filename, std::streamsize &size)
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
    unsigned char *sequence = new unsigned char[size + 1];

    for (int i = 0; i < size; i++)
    {
        sequence[i] = temp[i];
    }
    
    sequence[size] = '\0';
    return sequence;
}


vector<int64_t> compute_matching_statistics(int64_t *SA, int64_t *invSA, int64_t *LCP, int64_t seq_xy_len, int64_t seq_x_len, int64_t seq_y_len)
{
    vector<int64_t> best_len(seq_y_len, 0);
    int64_t start_y = seq_x_len + 1;


    int64_t last_x = -1;
    int64_t cur_min = INT_MAX;
    
    for (int64_t r = 0; r < seq_xy_len; r++)
    {
        int64_t pos = SA[r];
        if (last_x != -1)
        {
            cur_min = min(cur_min, LCP[r]);
        }
        if (pos < seq_x_len)
        {
            last_x = r;
            cur_min = INT_MAX;
        }
        else if (pos >= start_y && pos < start_y + seq_y_len)
        {
            if (last_x != -1)
            {
                int64_t y_i = pos - start_y;
                best_len[y_i] = max(best_len[y_i], cur_min);
            }
        }
    }

    last_x = -1;
    cur_min = INT_MAX;
    
    for (int64_t r = seq_xy_len - 1; r >= 0; r--)
    {
        int64_t pos = SA[r];
        if (last_x != -1 && r + 1 < seq_xy_len)
        {
            cur_min = min(cur_min, LCP[r + 1]);
        }
        if (pos < seq_x_len)
        {
            last_x = r;
            cur_min = INT_MAX;
        }
        else if (pos >= start_y && pos < start_y + seq_y_len)
        {
            if (last_x != -1)
            {
                int64_t y_i = pos - start_y;
                best_len[y_i] = max(best_len[y_i], cur_min);
            }
        }
    }

    return best_len;
}

int main(int argc, char **argv)
{
    string x, y, output;
    int opt;

    while ((opt = getopt(argc, argv, "x:y:o:")) != -1)
    {
        switch (opt)
        {
            case 'x': x = optarg; break;
            case 'y': y = optarg; break;
            case 'o': output = optarg; break;
            default: usage();
        }
    }

    if (x.empty() || y.empty() || output.empty())
    {
        fprintf(stderr, "Error: -x, -y, and -o options are required.\n");
        usage();
    }

    auto start = chrono::high_resolution_clock::now();

    streamsize size_x, size_y;
    unsigned char *seq_x = readFile(x.c_str(), size_x);
    unsigned char *seq_y = readFile(y.c_str(), size_y);

    int seq_x_len = strlen((char *)seq_x);
    int seq_y_len = strlen((char *)seq_y);

    unsigned char *seq_xy = concatenate(seq_x, seq_y, seq_x_len, seq_y_len);
    int64_t seq_xy_len = seq_x_len + seq_y_len + 1;

    int64_t *SA = (int64_t *)malloc(seq_xy_len * sizeof(int64_t));
    int64_t *invSA = (int64_t *)calloc(seq_xy_len, sizeof(int64_t));
    int64_t *LCP = (int64_t *)calloc(seq_xy_len, sizeof(int64_t));

    if (divsufsort64(seq_xy, SA, seq_xy_len) != 0)
    {
        fprintf(stderr, "Error: SA computation failed.\n");
        exit(EXIT_FAILURE);
    }

    for (int64_t i = 0; i < seq_xy_len; i++)
    {
        invSA[SA[i]] = i;
    }

    LCParray(seq_xy, seq_xy_len, SA, invSA, LCP);

    vector<int64_t> ms = compute_matching_statistics(SA, invSA, LCP, seq_xy_len, seq_x_len, seq_y_len);

    ofstream out(output);
    if (!out)
    {
        fprintf(stderr, "Error: Failed to open output file %s.\n", output.c_str());
        return EXIT_FAILURE;
    }

    for (int64_t i = 0; i < seq_y_len; i++)
    {
        out << ms[i] << "\n";
    }
    out.close();

    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Elapsed time: " << elapsed.count() << " ms\n";

    delete[] seq_x;
    delete[] seq_y;
    delete[] seq_xy;
    free(SA);
    free(LCP);
    free(invSA);

    return 0;
}
