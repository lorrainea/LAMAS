#include <stdio.h>
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


/*
 * Deletion-only set over the indices [0, m), all initially present, supporting
 *   erase(k)  : mark index k as deleted,
 *   next(k)   : smallest present index k' >= k (or m if none),
 *   prev(k)   : largest present index k' <= k (or -1 if none).
 * Implemented as two union-find forests with path compression (one pointing
 * right, one pointing left). This is the "interval union-find" problem; the
 * simple path-compression version below is near-linear in practice.
 */
struct DeletionOnlySet
{
    int m;
    vector<int> nxt;   /* nxt[k] == k iff k is present; nxt[m] = m is a sentinel      */
    vector<int> prv;   /* stored with an offset of 1: prv[k+1]; prv[0] = 0 is a sentinel (-1) */

    void init(int size)
    {
        m = size;
        nxt.resize(m + 1);
        prv.resize(m + 2);
        for (int k = 0; k <= m; ++k) nxt[k] = k;
        for (int k = 0; k <= m + 1; ++k) prv[k] = k;
    }

    int next(int k)
    {
        int r = k;
        while (nxt[r] != r) r = nxt[r];
        while (nxt[k] != r) { int t = nxt[k]; nxt[k] = r; k = t; }   /* path compression */
        return r;                                                   /* r == m means "none" */
    }

    int prev(int k)
    {
        int s = k + 1;
        int r = s;
        while (prv[r] != r) r = prv[r];
        while (prv[s] != r) { int t = prv[s]; prv[s] = r; s = t; }
        return r - 1;                                               /* -1 means "none" */
    }

    void erase(int k)
    {
        nxt[k] = k + 1;
        prv[k + 1] = k;
    }
};

/*
 * Linear-time LAMAS (Section 3.5 of the paper).
 *
 * Notation: S = X ~ Y; SA/invSA/LCP are those of S; A[k] = invSA[k] for k in [0,|X|) and
 * B[i] = invSA[|X|+1+i] for i in [0,|Y|). For every i we need the predecessor and the
 * successor of B[i] among the values A[max(0,i-d) .. min(|X|-1,i+d)]. Y is cut into blocks
 * of d positions; block b (Y positions [bd,(b+1)d)) only looks at the A positions
 * [(b-1)d, (b+2)d), which we index locally by j in [0,3d). The window of the local
 * position i' is [i', i'+2d] = [i', d) U [d, 2d) U [2d, 2d+i'], i.e. a suffix of the
 * "left" third, the whole "middle" third and a prefix of the "right" third. The middle
 * part is static; the left and right parts are deletion-only predecessor problems.
 * Sorting the elements of every block by rank is done by one global scan of SA.
 */
template<class RMQ>
int compute_lamas(int seq_x_len, int seq_y_len, int d, const int64_t * SA, const int64_t * invSA,
                  const int64_t * LCP, RMQ & rmq, const string & output_file)
{
    ofstream out(output_file);

    if (seq_y_len <= 0) return 0;

    /* LCP of the suffixes of S with (distinct) ranks r1 and r2 */
    auto lcp_ranks = [&](int64_t r1, int64_t r2) -> int64_t
    {
        if (r1 > r2) swap(r1, r2);
        return LCP[rmq(r1 + 1, r2)];
    };

    /* Report the answer for Y position i given the X positions holding the
       predecessor and the successor rank (-1 if undefined).                 */
    auto report = [&](int i, int64_t pred_pos, int64_t succ_pos)
    {
        int64_t rank = invSA[seq_x_len + 1 + i];
        int64_t lp = INT_MIN, ls = INT_MIN;
        int64_t li, delta;

        if (pred_pos >= 0) lp = lcp_ranks(rank, invSA[pred_pos]);
        if (succ_pos >= 0) ls = lcp_ranks(rank, invSA[succ_pos]);

        if (pred_pos < 0 && succ_pos < 0)
        {
            li = INT_MIN;
            delta = INT_MIN;
        }
        else
        {
            li = max(lp, ls);
            delta = (lp >= ls ? pred_pos : succ_pos) - i;
        }

        (void) delta;                 /* delta is computed for the two-column output below */
        //out<<li<<" "<<delta<<'\n';
        out<<li<<'\n';
    };

    /* Any d >= max(|X|,|Y|) yields the same windows as d = max(|X|,|Y|). */
    d = min(d, max(seq_x_len, seq_y_len));

    if (d == 0)
    {
        for (int i = 0; i < seq_y_len; ++i)
            report(i, (i < seq_x_len ? i : -1), -1);
        return 0;
    }

    const int nb = (seq_y_len + d - 1) / d;   /* number of blocks of Y */

    /* ---- Step A: bucket the elements of every block, sorted by rank, by one scan of SA. ---- */
    vector<int64_t> off(nb + 1, 0);
    for (int b = 0; b < nb; ++b)
    {
        int64_t aBase = (int64_t) (b - 1) * d;
        int64_t aL = max<int64_t>(0, aBase), aR = min<int64_t>(seq_x_len, aBase + 3 * (int64_t) d);
        int64_t yL = (int64_t) b * d, yR = min<int64_t>(seq_y_len, yL + d);
        off[b + 1] = off[b] + max<int64_t>(0, aR - aL) + (yR - yL);   /* blocks beyond the end of X have no A elements */
    }
    vector<int> buf(off[nb]);                 /* positions in S, grouped by block, sorted by rank */
    {
        vector<int64_t> fill(off.begin(), off.end() - 1);
        int64_t seq_xy_len = (int64_t) seq_x_len + seq_y_len + 1;
        for (int64_t r = 0; r < seq_xy_len; ++r)
        {
            int64_t p = SA[r];
            if (p < seq_x_len)
            {
                int b0 = (int) (p / d);
                for (int b = max(0, b0 - 1); b <= min(nb - 1, b0 + 1); ++b)
                    buf[fill[b]++] = (int) p;
            }
            else if (p > seq_x_len)
            {
                int b = (int) ((p - seq_x_len - 1) / d);
                buf[fill[b]++] = (int) p;
            }
        }
    }

    /* ---- Step B: solve every block. ---- */
    vector<int> rho(3 * (size_t) d);                     /* rank of A element (local index j) within its third */
    vector<int> cL, cM, cR;                              /* sorted (by rank) X positions of the three thirds   */
    vector<int> piL(d), piM(d), piR(d);                  /* # elements of each third smaller than B[i']        */
    vector<int> predL(d), succL(d), predR(d), succR(d);  /* answers for the left/right thirds (X positions)     */
    DeletionOnlySet ds;
    cL.reserve(d); cM.reserve(d); cR.reserve(d);

    for (int b = 0; b < nb; ++b)
    {
        int64_t aBase = (int64_t) (b - 1) * d;
        int64_t aL = max<int64_t>(0, aBase), aR = min<int64_t>(seq_x_len, aBase + 3 * (int64_t) d);
        int64_t yL = (int64_t) b * d, yR = min<int64_t>(seq_y_len, yL + d);
        int len = (int) (yR - yL);

        /* exists(j): does local A index j correspond to a position of X? */
        auto exists = [&](int j) -> bool { int64_t p = aBase + j; return p >= aL && p < aR; };

        cL.clear(); cM.clear(); cR.clear();

        for (int64_t t = off[b]; t < off[b + 1]; ++t)
        {
            int p = buf[t];
            if (p < seq_x_len)
            {
                int j = (int) (p - aBase);
                if      (j <  d)     { rho[j] = (int) cL.size(); cL.push_back(p); }
                else if (j < 2 * d)  { rho[j] = (int) cM.size(); cM.push_back(p); }
                else                 { rho[j] = (int) cR.size(); cR.push_back(p); }
            }
            else
            {
                int ip = (int) (p - seq_x_len - 1 - yL);
                piL[ip] = (int) cL.size(); piM[ip] = (int) cM.size(); piR[ip] = (int) cR.size();
            }
        }

        /* Right third: prefix [2d, 2d+i'], processed for decreasing i' with deletions. */
        ds.init((int) cR.size());
        for (int j = 2 * d + len; j < 3 * d; ++j)
            if (exists(j)) ds.erase(rho[j]);
        for (int ip = len - 1; ip >= 0; --ip)
        {
            int s = ds.next(piR[ip]);
            succR[ip] = (s < (int) cR.size()) ? cR[s] : -1;
            int q = (piR[ip] > 0) ? ds.prev(piR[ip] - 1) : -1;
            predR[ip] = (q >= 0) ? cR[q] : -1;
            if (exists(2 * d + ip)) ds.erase(rho[2 * d + ip]);
        }

        /* Left third: suffix [i', d), processed for increasing i' with deletions. */
        ds.init((int) cL.size());
        for (int ip = 0; ip < len; ++ip)
        {
            int s = ds.next(piL[ip]);
            succL[ip] = (s < (int) cL.size()) ? cL[s] : -1;
            int q = (piL[ip] > 0) ? ds.prev(piL[ip] - 1) : -1;
            predL[ip] = (q >= 0) ? cL[q] : -1;
            if (exists(ip)) ds.erase(rho[ip]);
        }

        /* Combine the three parts and report. */
        for (int ip = 0; ip < len; ++ip)
        {
            int64_t predM = (piM[ip] > 0) ? cM[piM[ip] - 1] : -1;
            int64_t succM = (piM[ip] < (int) cM.size()) ? cM[piM[ip]] : -1;

            int64_t pred_pos = -1, succ_pos = -1;
            for (int64_t cand : {(int64_t) predL[ip], predM, (int64_t) predR[ip]})
                if (cand >= 0 && (pred_pos < 0 || invSA[cand] > invSA[pred_pos])) pred_pos = cand;
            for (int64_t cand : {(int64_t) succL[ip], succM, (int64_t) succR[ip]})
                if (cand >= 0 && (succ_pos < 0 || invSA[cand] < invSA[succ_pos])) succ_pos = cand;

            report((int) (yL + ip), pred_pos, succ_pos);
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

    compute_lamas(seq_x_len, seq_y_len, d, SA, invSA, LCP, rmq, output);
    
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
