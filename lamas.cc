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
#include <sdsl/io.hpp>

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
int64_t LCParray ( unsigned char * text, int64_t n, int64_t * SA, int64_t * ISA, int64_t * LCP )
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


int compute_delta(rmq_succinct_sct<> &rmq, int &delta, int &li, int rp, int rs, int64_t * SA, int64_t * invSA, int64_t * LCP, int seq_x_len, int seq_y_len, int seq_xy_len, int i)
{

    int lp = 0;
    int ls = 0;
    
    int index = i-(seq_x_len+1);
        
    if( rp == -1 )
        lp = INT_MIN;
        
    if( rs == -1 )
        ls = INT_MIN;
            
    if( rp > -1 )
    {
        lp = LCP[rmq ( min((int)invSA[i],rp)+1, max(rp,(int)invSA[i]) )];
    }
            
    if( rs > -1 )
    {
        ls = LCP[rmq ( min((int)invSA[i],rs)+1, max(rs,(int)invSA[i]) )];
    }
        
    li = max(lp, ls);
    delta = (lp >= ls) ? SA[rp] - index : SA[rs] - index;
    
    return 0;
}
    
int compute_rp(set<int> &W, int64_t * invSA, int seq_x_len, int seq_xy_len, int i)
{

    int rp = 0;
    
    auto it = W.upper_bound(invSA[i]);
    if (it == W.begin())
    {
        rp = -1;
    }
    else
    {
        --it;
        rp = *it;
    }
    
    return rp;
}


int compute_rs(set<int> &W, int64_t * invSA, int seq_x_len, int seq_xy_len, int i)
{
    
    int rs = 0;
    
    auto it = W.upper_bound(invSA[i]);
    if (it == W.end())
    {
       rs = -1;
    }
    else rs = *it;
    
    return rs;
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

int compute_W(set<int> &W, int seq_x_len, int64_t * invSA, int d, int i)
{

    for(int j=max(0,i-d); j<=min(seq_x_len-1,i+d); j++)
        W.insert(invSA[j]);
       
    return 0;
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
	
	for ( int64_t i = 0; i < seq_xy_len; i ++ )
	{
        invSA [SA[i]] = i;        
	}
    
	LCParray( seq_xy, seq_xy_len, SA, invSA, LCP );
	
	int_vector<> lcp_rmq( seq_xy_len , 0 );

	for ( int i = 0; i < seq_xy_len; i ++ )
	{
		lcp_rmq[i] = LCP[i];
	}
	
	rmq_succinct_sct<> rmq(&lcp_rmq);
    ofstream out(output);
    
    set<int> W;
    int rp = 0;
    int rs = 0;
    int delta = 0;
    int li = 0;
    
    for(int i = 0; i<seq_y_len; i++)
    { 
        int index = seq_x_len+1+i;
        compute_W(W, seq_x_len, invSA, d, i);
        rp = compute_rp(W, invSA, seq_x_len, seq_xy_len, index);
        rs = compute_rs(W, invSA, seq_x_len, seq_xy_len, index);
        compute_delta(rmq, delta, li, rp, rs, SA, invSA, LCP, seq_x_len, seq_y_len, seq_xy_len, index);
        
        out<<li<<endl;
        
        W.clear();
    }
    
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
