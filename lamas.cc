#include <stdio.h>
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


int compute_delta(int * delta, int * li, int * rp, int * rs, int64_t * SA, int64_t * invSA, int64_t * LCP, int seq_x_len, int seq_y_len, int seq_xy_len)
{

    int * lp = new int[seq_y_len];
    int * ls = new int[seq_y_len];
    
    int_vector<> lcp_rmq( seq_xy_len , 0 );

	for ( int i = 0; i < seq_xy_len; i ++ )
	{
		lcp_rmq[i] = LCP[i];
	}
	
	rmq_succinct_sct<> rmq(&lcp_rmq);
    
    for(int i = seq_x_len+1;  i<seq_xy_len; i++)
    {
        int index = i-(seq_x_len+1);
        
        if( rp[index] == -1 )
            lp[index] = INT_MIN;
        
        if( rs[index] == -1 )
            ls[index] = INT_MIN;
            
        if( rp[index] > -1 )
        {
            lp[index] = lcp_rmq[rmq ( min((int)invSA[i],rp[index])+1, max(rp[index],(int)invSA[i]) )];
        }
            
        if( rs[index] > -1 )
        {
            ls[index] = lcp_rmq[rmq ( min((int)invSA[i],rs[index])+1, max(rs[index],(int)invSA[i]) )];
        }
        
        li[index] = max(lp[index], ls[index]);
        delta[index] = (lp[index] >= ls[index]) ? SA[rp[index]] - index: SA[rs[index]] - index;
    }
    
   /* cout<<"ls "<<endl;
    for(int i = 0;  i<seq_y_len; i++)
    {
        cout<<ls[i]<<" ";
    }
    cout<<endl;
    
    cout<<"lp "<<endl;
    for(int i = 0;  i<seq_y_len; i++)
    {
        cout<<lp[i]<<" ";
    }
    cout<<endl;
    
    cout<<"delta "<<endl;
    for(int i = 0;  i<seq_y_len; i++)
    {
        cout<<delta[i]<<" ";
    }
    cout<<endl;*/
    
    delete[] lp;
    delete[] ls;
    
    return 0;
}
    
int compute_rp(set<int> ** W, int * rp, int64_t * invSA, int seq_x_len, int seq_y_len, int seq_xy_len)
{

    for(int i = seq_x_len+1;  i<seq_xy_len; i++)
    {
        int index = i-(seq_x_len+1);
        
        auto it = W[index]->upper_bound(invSA[i]);
        if (it == W[index]->begin())
        {
            rp[index] = -1;
        }
        else
        {
            --it;
            rp[index] = *it;
        }
    }
    
    return 0;
}


int compute_rs(set<int> ** W, int * rs, int64_t * invSA, int seq_x_len, int seq_y_len, int seq_xy_len)
{
    
    for(int i = seq_x_len+1;  i<seq_xy_len; i++)
    {
        int index = i-(seq_x_len+1);
        
        auto it = W[index]->upper_bound(invSA[i]);
        if (it == W[index]->end())
        {
            rs[index] = -1;
        }
        else rs[index] = *it;
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

int compute_W(set<int> ** W, int seq_x_len, int seq_y_len, int64_t * invSA, int d)
{

    for(int i = 0; i<seq_y_len; i++)
    {
        for(int j=max(0,i-d); j<=min(seq_x_len-1,i+d); j++)
        {
            W[i]->insert(invSA[j]);
        }
    }

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

    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);

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
	
    set<int> * W[seq_y_len];
    
    for(int i = 0; i<seq_y_len; i++)
            W[i] = new set<int>;
    
    compute_W(W, seq_x_len, seq_y_len, invSA, d);
    
   /* cout<<" W "<<endl;
    for(int a= 0; a<seq_y_len; a++)
    {
         for (auto x : *W[a])
         {
            cout << x << " ";
         }
         cout<<endl;
    }*/
    
    int * rp = new int[seq_y_len];
    int * rs = new int[seq_y_len];
    
    compute_rp(W, rp, invSA, seq_x_len, seq_y_len, seq_xy_len);
    compute_rs(W, rs, invSA, seq_x_len, seq_y_len, seq_xy_len);

    /* cout<<" rp "<<endl;
    for(int a= 0; a<seq_y_len; a++)
    {
        
         cout<<rp[a]<<endl;
    }
    
    
    cout<<" rs "<<endl;
    for(int a= 0; a<seq_y_len; a++)
    {
        
         cout<<rs[a]<<endl;
    }*/
    
    int * delta = new int[seq_y_len];
    int * li = new int[seq_y_len];
    
    compute_delta(delta, li, rp, rs, SA, invSA, LCP, seq_x_len, seq_y_len, seq_xy_len);
        
    ofstream out(output);

    for(int i = 0; i<seq_y_len; i++)
        out<<li[i]<<" "<<delta[i]<<endl;
        
    cout<< "Elapsed time: " << elapsed.count() << " ms\n";

    delete[] seq_x;
    delete[] seq_y;
    delete[] seq_xy;
    delete[] rp;
    delete[] rs;
    delete[] delta;
    delete[] li;
    free( SA );
    free( LCP );
    free( invSA );
   
    for (int i = 0; i < seq_y_len; i++)
        delete W[i];
   
    return 0;
}
