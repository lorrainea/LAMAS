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
#include <divsufsort64.h>             
#include <sdsl/bit_vectors.hpp>               
#include <sdsl/rmq_support.hpp>	
#include <sdsl/rmq_support_sparse_table.hpp>

using namespace std;
using namespace sdsl;

void usage(void)
{
    fprintf(stderr, "lamas\n");
    fprintf(stderr, "  -x <seq_x_filename> input sequence file.\n");
    fprintf(stderr, "  -y <seq_y_filename> input sequence file.\n");
    fprintf(stderr, "  -d <d_range> range of d.\n");
    fprintf(stderr, "  -o <output_filename> output file\n");
    exit(0);
}

/* Kasai et al algorithm for O(n)-time LCP construction */
int64_t LCParray (unsigned char * text, int64_t n, int64_t * SA, int64_t * ISA, int64_t * LCP)
{
    int64_t i=0, j=0;

    LCP[0] = 0;
    for ( i = 0; i < n; i++ )
    {
        if ( ISA[i] != 0 )
        {
            if ( i == 0 )
            {
                j = 0;
            }
            else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
            
            while ( text[i+j] == text[SA[ISA[i]-1]+j] )
            {
                j++;
            }
            LCP[ISA[i]] = j;
        }
    }
   
    return 0;
}


int compute_delta(rmq_succinct_sct<>& rmq, int64_t& delta, int64_t& li, int64_t rp, int64_t rs, int64_t* SA, int64_t* invSA, int64_t* LCP, int seq_x_len, int index)
{
    int64_t rank = invSA[index];

    int lp = INT_MIN;
    int ls = INT_MIN;

    if(rp >= 0)
    {
        lp = LCP[rmq(min(rank, rp) + 1, max(rank, rp))];
    }

    if(rs >= 0)
    {
        ls = LCP[rmq(min(rank, rs) + 1, max(rank, rs))];
    }

    if(rp == -1 && rs == -1)
    {
        li = INT_MIN;
        delta = INT_MIN;
    }
    else
    {
        li = max(lp, ls);

        if(lp >= ls)
        {
            delta = SA[rp] - (index - (seq_x_len + 1));
        }
        else
        {
            delta = SA[rs] - (index - (seq_x_len + 1));
        }
    }
    
    return 0;
}

int compute_w(int seq_x_len, int seq_y_len, int d, int64_t * SA, int64_t * invSA, int64_t * LCP, rmq_succinct_sct<>& rmq, const string& output_file)
{
    ofstream out(output_file);

    set<uint64_t> W;

    int left  = max(0, -d);
    int right = min(seq_x_len - 1, d);

    for(int j = left; j <= right; ++j)
        W.insert(invSA[j]);
  
    for(int i = 0; i < seq_y_len; ++i)
    {
        int index = seq_x_len + 1 + i;
        uint64_t rank  = invSA[index];

        int64_t rp = -1;
        int64_t rs = -1;

        auto it = W.upper_bound(rank);

        if(it != W.begin())
        {
            auto pred = prev(it);
            rp = *pred;
        }

        if(it != W.end())
        {
            rs = *it;
        }

        int64_t delta;
        int64_t li;

        compute_delta(rmq, delta, li, rp, rs, SA, invSA, LCP, seq_x_len, index);

        //out<<li<<" "<<delta<<'\n';
        
        out<<li<<'\n';

        int remove_pos = i-d;

        if(remove_pos >= 0 && remove_pos < seq_x_len)
        {
            W.erase(invSA[remove_pos]);
        }

        int add_pos = i + d + 1;

        if(add_pos >= 0 && add_pos < seq_x_len)
        {
            W.insert(invSA[add_pos]);
        }
    }
    
    return 0;
} 

unsigned char * concatenate(unsigned char * seq_x, unsigned char * seq_y, int seq_x_len, int seq_y_len)
{
    unsigned char * seq_xy = new unsigned char[seq_x_len + seq_y_len + 2];
   
    strcpy(reinterpret_cast<char*>(seq_xy), reinterpret_cast<char*>(seq_x));

    int pos = seq_x_len;

    seq_xy[pos] = '~';
    seq_xy[pos+1] = '\0';
    
    strcat(reinterpret_cast<char*>(seq_xy),reinterpret_cast<char*>(seq_y));
    
    seq_xy[seq_x_len + seq_y_len + 1] = '\0';
    
    return seq_xy;
}

unsigned char* readFile(const char * filename, std::streamsize& size)
{
    ifstream file(filename, ios::binary);

    if (!file)
    {
        fprintf(stderr, "Failed to open %s.\n", filename);
        exit( EXIT_FAILURE );
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

    int d = -1;

    int opt;

    while ((opt = getopt(argc, argv, "x:y:d:o:")) != -1)
    {
        switch (opt)
        {
            case 'x':
                x = optarg;
                break;

            case 'y':
                y = optarg;
                break;

            case 'd':
            {
                stringstream ss(optarg);
                ss >> d;
                break;
            }

            case 'o':
                output = optarg;
                break;

            default:
                usage();
        }
    }

    if (d < 0)
    {
        fprintf(stderr, "Error: distance (-d) cannot be less than 0.\n");
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
    
    unsigned char * seq_xy = concatenate( seq_x, seq_y, seq_x_len, seq_y_len);
    
    int seq_xy_len = seq_x_len + seq_y_len + 1;
    
	int64_t * SA;
	int64_t * LCP;
	int64_t * invSA;
	
	SA = ( int64_t * ) malloc( ( seq_xy_len ) * sizeof( int64_t ) );
	invSA = ( int64_t * ) calloc( seq_xy_len , sizeof( int64_t ) );
	LCP = ( int64_t * ) calloc  ( seq_xy_len, sizeof( int64_t ) );

  	if( divsufsort64( seq_xy, SA,  seq_xy_len ) != 0 )
  	{
  		fprintf(stderr, "Error: SA computation failed.\n" );
        exit( EXIT_FAILURE );
  	}
	
	for (int64_t i=0; i <seq_xy_len; i++)
	    invSA[SA[i]] = i;        
	
	LCParray(seq_xy, seq_xy_len, SA, invSA, LCP);
	
	int_vector<> lcp_rmq(seq_xy_len , 0);

	for ( int i = 0; i < seq_xy_len; i++ )
	    lcp_rmq[i] = LCP[i];
	
	rmq_succinct_sct<> rmq(&lcp_rmq);
    ofstream out(output);
    
    int rp = 0;
    int rs = 0;
    int delta = 0;
    int li = 0;
    
    compute_w(seq_x_len, seq_y_len, d, SA, invSA, LCP, rmq, output);
    
    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout<< "Elapsed time: " << elapsed.count() << " ms\n";

    delete[] seq_x;
    delete[] seq_y;
    delete[] seq_xy;
    free( SA );
    free( LCP );
    free( invSA );
   
    return 0;
}
