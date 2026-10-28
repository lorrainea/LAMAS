#include <bits/stdc++.h>
#include <unistd.h>

using namespace std;

/*
    ============================================================
    Usage
    ============================================================
*/
void usage(void)
{
    fprintf(stderr, "linear_matching_statistics\n");
    fprintf(stderr, "  -x <seq_x_filename> input sequence S.\n");
    fprintf(stderr, "  -y <seq_y_filename> input sequence T.\n");
    fprintf(stderr, "  -o <output_filename> output file.\n");
    exit(EXIT_FAILURE);
}

/*
    ============================================================
    Read sequence from file
    ============================================================
*/
string readFile(const char* filename)
{
    ifstream file(filename, ios::binary);

    if (!file)
    {
        fprintf(stderr, "Failed to open %s.\n", filename);
        exit(EXIT_FAILURE);
    }

    string sequence;
    char ch;

    while (file.get(ch))
    {
        if (ch == ' ' ||
            ch == '\n' ||
            ch == '\r' ||
            ch == '\t')
        {
            continue;
        }

        sequence.push_back(ch);
    }

    return sequence;
}

/*
    ============================================================
    SA-IS suffix array construction - O(N)
    ============================================================
*/
template <class T>
vector<int> sa_is(const vector<T>& s, int upper)
{
    int n = (int)s.size();

    if (n == 0)
        return {};

    if (n == 1)
        return {0};

    if (n == 2)
    {
        if (s[0] < s[1])
            return {0, 1};
        else
            return {1, 0};
    }

    vector<int> sa(n);
    vector<bool> ls(n);

    ls[n - 1] = true;

    for (int i = n - 2; i >= 0; i--)
    {
        ls[i] =
            (s[i] == s[i + 1])
                ? ls[i + 1]
                : (s[i] < s[i + 1]);
    }

    int alph = upper + 2;
    vector<int> sum_l(alph, 0);
    vector<int> sum_s(alph, 0);

    for (int i = 0; i < n; i++)
    {
        if (!ls[i])
            sum_s[s[i]]++;
        else
            sum_l[s[i] + 1]++;
    }

    for (int i = 0; i < alph; i++)
    {
        sum_s[i] += sum_l[i];

        if (i + 1 < alph)
            sum_l[i + 1] += sum_s[i];
    }

    auto induce = [&](const vector<int>& lms)
    {
        fill(sa.begin(), sa.end(), -1);

        vector<int> buf(alph);

        copy(sum_s.begin(),
             sum_s.end(),
             buf.begin());

        for (auto d : lms)
        {
            if (d == n)
                continue;

            sa[buf[s[d]]++] = d;
        }

        copy(sum_l.begin(),
             sum_l.end(),
             buf.begin());

        sa[buf[s[n - 1]]++] = n - 1;

        for (int i = 0; i < n; i++)
        {
            int v = sa[i];

            if (v >= 1 && !ls[v - 1])
            {
                sa[buf[s[v - 1]]++] = v - 1;
            }
        }

        copy(sum_l.begin(),
             sum_l.end(),
             buf.begin());

        for (int i = n - 1; i >= 0; i--)
        {
            int v = sa[i];

            if (v >= 1 && ls[v - 1])
            {
                sa[--buf[s[v - 1] + 1]] = v - 1;
            }
        }
    };

    vector<int> lms_map(n + 1, -1);
    int m = 0;

    for (int i = 1; i < n; i++)
    {
        if (!ls[i - 1] && ls[i])
        {
            lms_map[i] = m++;
        }
    }

    // FIX: Handle case where no LMS suffixes exist (monotone sequences)
    if (m == 0)
    {
        iota(sa.begin(), sa.end(), 0);
        return sa;
    }

    vector<int> lms;
    lms.reserve(m);

    for (int i = 1; i < n; i++)
    {
        if (!ls[i - 1] && ls[i])
        {
            lms.push_back(i);
        }
    }

    induce(lms);

    vector<int> sorted_lms;
    sorted_lms.reserve(m);

    for (int v : sa)
    {
        if (v >= 0 && lms_map[v] != -1)
        {
            sorted_lms.push_back(v);
        }
    }

    vector<int> rec_s(m);
    int rec_upper = 0;

    rec_s[lms_map[sorted_lms[0]]] = 0;

    for (int i = 1; i < m; i++)
    {
        int l = sorted_lms[i - 1];
        int r = sorted_lms[i];

        int end_l =
            (lms_map[l] + 1 < m)
                ? lms[lms_map[l] + 1]
                : n;

        int end_r =
            (lms_map[r] + 1 < m)
                ? lms[lms_map[r] + 1]
                : n;

        bool same = true;

        if (end_l - l != end_r - r)
        {
            same = false;
        }
        else
        {
            while (l < end_l)
            {
                if (s[l] != s[r])
                {
                    same = false;
                    break;
                }

                l++;
                r++;
            }
        }

        if (!same)
            rec_upper++;

        rec_s[lms_map[sorted_lms[i]]] = rec_upper;
    }

    vector<int> rec_sa =
        sa_is(rec_s, rec_upper);

    for (int i = 0; i < m; i++)
    {
        sorted_lms[i] = lms[rec_sa[i]];
    }

    induce(sorted_lms);

    return sa;
}

/*
    ============================================================
    Suffix array wrapper
    ============================================================
*/
vector<int> suffix_array(const vector<int>& s)
{
    if (s.empty())
        return {};

    int upper =
        *max_element(s.begin(), s.end());

    return sa_is(s, upper);
}

/*
    ============================================================
    Kasai LCP construction - O(N)
    ============================================================
*/
vector<int> build_lcp(
    const vector<int>& text,
    const vector<int>& sa)
{
    int n = text.size();
    vector<int> rank(n);

    for (int i = 0; i < n; i++)
    {
        rank[sa[i]] = i;
    }

    vector<int> lcp(n, 0);
    int h = 0;

    for (int i = 0; i < n; i++)
    {
        int r = rank[i];

        if (r == 0)
            continue;

        int j = sa[r - 1];

        while (i + h < n &&
               j + h < n &&
               text[i + h] == text[j + h])
        {
            h++;
        }

        lcp[r] = h;

        if (h > 0)
            h--;
    }

    return lcp;
}

/*
    ============================================================
    Matching Statistics - O(|S| + |T|)
    ============================================================
*/
vector<int> matching_statistics(
    const string& S,
    const string& T)
{
    int n = S.size();
    int m = T.size();

    vector<int> text;
    text.reserve(n + m + 1);

    for (unsigned char c : S)
    {
        text.push_back((int)c + 1);
    }

    text.push_back(0);

    for (unsigned char c : T)
    {
        text.push_back((int)c + 1);
    }

    vector<int> sa = suffix_array(text);
    vector<int> lcp = build_lcp(text, sa);

    int N = text.size();
    vector<int> ms(n, 0);

    // Left-to-right scan
    int lastT = -1;
    int currentMin = INT_MAX;

    for (int r = 0; r < N; r++)
    {
        int pos = sa[r];

        if (lastT != -1)
        {
            currentMin = min(currentMin, lcp[r]);
        }

        if (pos < n)
        {
            if (lastT != -1)
            {
                ms[pos] = max(ms[pos], currentMin);
            }
        }

        bool isT =
            (pos >= n + 1 &&
             pos < n + 1 + m);

        if (isT)
        {
            lastT = r;
            currentMin = INT_MAX;
        }
    }

    // Right-to-left scan
    int nextT = -1;
    currentMin = INT_MAX;

    for (int r = N - 1; r >= 0; r--)
    {
        int pos = sa[r];

        // Bounds check: ensure r + 1 is valid
        if (nextT != -1 && r + 1 < N)
        {
            currentMin = min(currentMin, lcp[r + 1]);
        }

        if (pos < n)
        {
            if (nextT != -1)
            {
                ms[pos] = max(ms[pos], currentMin);
            }
        }

        bool isT =
            (pos >= n + 1 &&
             pos < n + 1 + m);

        if (isT)
        {
            nextT = r;
            currentMin = INT_MAX;
        }
    }

    return ms;
}

/*
    ============================================================
    MAIN
    ============================================================
*/
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

    if (x.empty())
    {
        fprintf(stderr, "Error: input file (-x) is required.\n");
        usage();
    }

    if (y.empty())
    {
        fprintf(stderr, "Error: input file (-y) is required.\n");
        usage();
    }

    if (output.empty())
    {
        fprintf(stderr, "Error: output file (-o) is required.\n");
        usage();
    }

    auto start = chrono::high_resolution_clock::now();

    string S = readFile(x.c_str());
    string T = readFile(y.c_str());

    if (S.empty())
    {
        fprintf(stderr, "Error: sequence X is empty.\n");
        return EXIT_FAILURE;
    }

    if (T.empty())
    {
        fprintf(stderr, "Error: sequence Y is empty.\n");
        return EXIT_FAILURE;
    }

    vector<int> ms = matching_statistics(S, T);

    ofstream out(output);

cout<<" done "<<endl;
    if (!out)
    {
        fprintf(stderr, "Failed to open output file %s.\n", output.c_str());
        return EXIT_FAILURE;
    }

    for (size_t i = 0; i < ms.size(); i++)
    {
        out << ms[i] << '\n';
    }

    out.close();

    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << "Elapsed time: " << elapsed.count() << " ms\n";
    cout << "Output written to: " << output << '\n';

    return 0;
}
