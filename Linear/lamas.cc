#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
#include <divsufsort.h>
#include <sdsl/rmq_support.hpp>

namespace lamas {

using index_type = saidx_t;
static_assert(sizeof(index_type) == 4, "LAMAS requires 32-bit divsufsort indices");
constexpr index_type no_match = INT_MIN; // Preserve the existing empty-window output.

/* Deletion-only predecessor/successor with sentinels. Each index is erased
 * once; the two forests compress paths. This is the same practical interval
 * union-find scheme as the original implementation, not a new strict-linear
 * worst-case interval-union-find construction. */
struct DeletionOnlySet {
    std::vector<index_type> nxt, prv;
    void init(index_type size) {
        nxt.resize(static_cast<size_t>(size) + 1);
        prv.resize(static_cast<size_t>(size) + 1);
        std::iota(nxt.begin(), nxt.end(), 0);
        std::iota(prv.begin(), prv.end(), 0);
    }
    index_type next(index_type k) {
        index_type r = k;
        while (nxt[r] != r) r = nxt[r];
        while (nxt[k] != r) { const auto t = nxt[k]; nxt[k] = r; k = t; }
        return r;
    }
    index_type prev(index_type k) {
        index_type s = k + 1, r = s;
        while (prv[r] != r) r = prv[r];
        while (prv[s] != r) { const auto t = prv[s]; prv[s] = r; s = t; }
        return r - 1;
    }
    void erase(index_type k) {
        nxt[k] = k + 1;
        prv[k + 1] = k;
    }
};

struct SuffixIndex {
    std::vector<index_type> sa, isa, lcp;
    static constexpr index_type rmq_block = 32;
    std::vector<index_type> block_min;
    sdsl::rmq_succinct_sct<> rmq;

    explicit SuffixIndex(const std::vector<unsigned char>& text) {
        if (text.empty() || text.size() > static_cast<size_t>(INT_MAX))
            throw std::length_error("combined sequence length must be in [1, INT_MAX]");
        const auto n = static_cast<index_type>(text.size());
        sa.resize(n);
        if (divsufsort(text.data(), sa.data(), n) != 0)
            throw std::runtime_error("suffix-array construction failed");
        isa.resize(n);
        for (index_type r = 0; r < n; ++r) isa[sa[r]] = r;
        lcp.resize(n);
        // Kasai: retain h in a register instead of reloading the preceding LCP.
        index_type h = 0;
        for (index_type i = 0; i < n; ++i) {
            const auto r = isa[i];
            if (r == 0) { h = 0; continue; }
            const auto j = sa[r - 1];
            const auto limit = n - std::max(i, j);
            while (h < limit && text[i + h] == text[j + h]) ++h;
            lcp[r] = h;
            if (h) --h;
        }
        // RMQ only over minima of fixed-size blocks. An LCE scans at most two
        // boundary blocks and makes at most one succinct query, still O(1).
        // This shrinks RMQ construction and improves its cache locality.
        block_min.resize(1 + (n - 1) / rmq_block, INT_MAX);
        for (index_type r = 0; r < n; ++r)
            block_min[r / rmq_block] = std::min(block_min[r / rmq_block], lcp[r]);
        rmq = sdsl::rmq_succinct_sct<>(&block_min);
    }

    index_type lce(index_type a, index_type b, index_type best = -1) const {
        if (a > b) std::swap(a, b);
        ++a; // LCP minimum on (a_old,b]. The ranks are always distinct.
        // Both endpoints are upper bounds on the minimum. A candidate that
        // cannot improve the current answer needs no RMQ query.
        index_type bound = std::min(lcp[a], lcp[b]);
        if (bound <= best || a == b) return bound;
        // A fixed-size scan avoids succinct-RMQ overhead on very short ranges.
        if (b - a < 8) {
            for (index_type r = a + 1; r < b; ++r) bound = std::min(bound, lcp[r]);
            return bound;
        }
        const index_type first = a / rmq_block, last = b / rmq_block;
        if (first == last) {
            for (index_type r = a + 1; r < b; ++r) bound = std::min(bound, lcp[r]);
            return bound;
        }
        const index_type left_end = (first + 1) * rmq_block; // <= b, no overflow
        for (index_type r = a + 1; r < left_end; ++r) bound = std::min(bound, lcp[r]);
        for (index_type r = last * rmq_block; r < b; ++r) bound = std::min(bound, lcp[r]);
        if (bound <= best || first + 1 == last) return bound;
        return std::min(bound, block_min[rmq(first + 1, last - 1)]);
    }
};

/* For each Y position emit (i, predecessor rank, successor rank).
 * Missing ranks are -1 and n. The window is inclusive: [i-d,i+d].
 *
 * X and Y are each divided into physical blocks of d positions. One SA scan
 * stores each relevant suffix rank ONCE, in its own sorted block. Query block
 * b reuses X blocks b-1, b, b+1 (left suffix, static middle, right prefix).
 * Three linear merges give insertion points in those lists; the two moving
 * parts use deletion-only sets. All block-local arrays are reused.
 */
template<class Emit>
void neighbors(index_type nx, index_type ny, index_type d,
               const std::vector<index_type>& sa,
               const std::vector<index_type>& isa, Emit emit) {
    const auto n = static_cast<index_type>(sa.size());
    if (ny == 0) return;
    if (d == 0) {
        for (index_type i = 0; i < ny; ++i) {
            const auto r = isa[nx + 1 + i];
            const auto a = i < nx ? isa[i] : -1;
            emit(i, a >= 0 && a < r ? a : -1, a > r ? a : n);
        }
        return;
    }
    d = std::min(d, std::max(nx, ny));
    // X suffixes starting beyond the last possible window are never needed.
    const auto used_x = static_cast<index_type>(std::min<int64_t>(nx, int64_t(ny) + d));
    const auto xb = used_x ? 1 + (used_x - 1) / d : 0;
    const auto yb = 1 + (ny - 1) / d;
    std::vector<index_type> sorted(static_cast<size_t>(used_x) + ny);
    {
        std::vector<index_type> fill(static_cast<size_t>(xb) + yb);
        for (index_type b = 0; b < xb; ++b) fill[b] = b * d;
        for (index_type b = 0; b < yb; ++b) fill[xb + b] = used_x + b * d;
        for (index_type r = 0; r < n; ++r) {
            const auto p = sa[r];
            if (p < used_x) sorted[fill[p / d]++] = r;
            else if (p > nx) sorted[fill[xb + (p - nx - 1) / d]++] = r;
        }
    }
    const auto capacity = std::min(d, ny);
    std::vector<index_type> pred(capacity), succ(capacity), pi(capacity);
    std::vector<index_type> rho(std::min(d, used_x));
    DeletionOnlySet ds;

    for (index_type b = 0; b < yb; ++b) {
        const auto ystart = b * d;
        const auto len = std::min(d, ny - ystart);
        const auto* queries = sorted.data() + used_x + ystart;
        std::fill_n(pred.begin(), len, -1);
        std::fill_n(succ.begin(), len, n);

        // Merge queries with one sorted X block. For a moving part also record
        // each X position's local rank and each Y position's insertion point.
        auto part = [&](int64_t xstart64, int kind) {
            if (xstart64 < 0 || xstart64 >= used_x) return;
            const auto xstart = static_cast<index_type>(xstart64);
            const auto count = std::min(d, used_x - xstart);
            const auto* ranks = sorted.data() + xstart;
            if (kind != 0) {
                for (index_type k = 0; k < count; ++k)
                    rho[sa[ranks[k]] - xstart] = k;
            }
            index_type k = 0;
            for (index_type q = 0; q < len; ++q) {
                const auto rank = queries[q];
                const auto ip = sa[rank] - nx - 1 - ystart;
                while (k < count && ranks[k] < rank) ++k;
                if (kind == 0) {
                    if (k) pred[ip] = ranks[k - 1];
                    if (k < count) succ[ip] = ranks[k];
                } else pi[ip] = k;
            }
            if (kind == 0) return;
            ds.init(count);
            auto query = [&](index_type ip) {
                const auto p = ds.prev(pi[ip] - 1);
                const auto s = ds.next(pi[ip]);
                if (p >= 0) pred[ip] = std::max(pred[ip], ranks[p]);
                if (s < count) succ[ip] = std::min(succ[ip], ranks[s]);
            };
            if (kind < 0) { // left suffix [ip,d)
                for (index_type ip = 0; ip < len; ++ip) {
                    query(ip);
                    if (ip < count) ds.erase(rho[ip]);
                }
            } else { // right prefix [0,ip], including its last position
                for (index_type j = len; j < count; ++j) ds.erase(rho[j]);
                for (index_type ip = len; ip-- > 0;) {
                    query(ip);
                    if (ip < count) ds.erase(rho[ip]);
                }
            }
        };
        part(ystart, 0);
        part(int64_t(ystart) - d, -1);
        part(int64_t(ystart) + d, 1);
        for (index_type ip = 0; ip < len; ++ip)
            emit(ystart + ip, pred[ip], succ[ip]);
    }
}

template<class Emit>
void solve_indexed(index_type nx, index_type ny, index_type d,
                   const SuffixIndex& index, Emit emit) {
    const auto n = static_cast<index_type>(index.sa.size());
    neighbors(nx, ny, d, index.sa, index.isa,
        [&](index_type i, index_type pred, index_type succ) {
            if (pred < 0 && succ == n) { emit(no_match); return; }
            const auto rank = index.isa[nx + 1 + i];
            index_type best = -1;
            if (pred >= 0) best = index.lce(pred, rank);
            if (succ < n) best = std::max(best, index.lce(rank, succ, best));
            emit(best);
        });
}

/* d=0 needs only aligned runs, so bypass suffix preprocessing entirely. */
template<class Emit>
void solve_zero(const std::vector<unsigned char>& text, index_type nx,
                index_type ny, Emit emit) {
    index_type i = 0;
    const auto overlap = std::min(nx, ny);
    while (i < overlap) {
        if (text[i] != text[nx + 1 + i]) { emit(0); ++i; continue; }
        const auto start = i;
        while (i < overlap && text[i] == text[nx + 1 + i]) ++i;
        for (index_type run = i - start; run > 0; --run) emit(run);
    }
    while (i++ < ny) emit(no_match);
}

/* With a fixed small radius, a bounded number of diagonal scans is faster
 * than suffix preprocessing. Along diagonal delta, LCE(i,i+delta) is one plus
 * the next LCE if the symbols match, and zero otherwise. The threshold is a
 * constant, so this branch remains O(nx+ny) time and uses only ny integers. */
constexpr index_type direct_radius = 16;
template<class Emit>
void solve_small(const std::vector<unsigned char>& text, index_type nx,
                 index_type ny, index_type d, Emit emit) {
    if (d < 0 || d > direct_radius)
        throw std::invalid_argument("direct solver radius exceeds fixed threshold");
    std::vector<index_type> answer(ny, no_match);
    for (index_type delta = -d; delta <= d; ++delta) {
        const auto begin = std::max<index_type>(0, -delta);
        const auto end = static_cast<index_type>(std::min<int64_t>(ny, int64_t(nx) - delta));
        index_type length = 0;
        for (index_type i = end; i > begin;) {
            --i;
            length = text[i + delta] == text[nx + 1 + i] ? length + 1 : 0;
            answer[i] = std::max(answer[i], length);
        }
    }
    for (const auto value : answer) emit(value);
}

// Append filtered input directly to the one concatenated text buffer. Newline
// is a safe separator because all input newlines are removed. No strlen/strcpy:
// embedded NUL and '~' are ordinary symbols and cannot truncate/cross X.
void read_append(const std::string& filename, std::vector<unsigned char>& text) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) throw std::runtime_error("failed to open input: " + filename);
    std::array<char, 65536> buffer;
    while (in.read(buffer.data(), buffer.size()) || in.gcount()) {
        for (std::streamsize i = 0; i < in.gcount(); ++i) {
            const auto c = static_cast<unsigned char>(buffer[i]);
            if (c == ' ' || c == '\n' || c == '\r' || c == '\t') continue;
            if (text.size() == static_cast<size_t>(INT_MAX))
                throw std::length_error("combined sequence length exceeds INT_MAX");
            text.push_back(c);
        }
    }
    if (in.bad() || !in.eof()) throw std::runtime_error("failed to read input: " + filename);
}

class Output {
    std::ofstream out;
    std::array<char, 65536> buffer;
    size_t used = 0;
    void flush() {
        out.write(buffer.data(), static_cast<std::streamsize>(used));
        if (!out) throw std::runtime_error("failed to write output");
        used = 0;
    }
public:
    explicit Output(const std::string& filename) : out(filename, std::ios::binary) {
        if (!out) throw std::runtime_error("failed to open output: " + filename);
    }
    void operator()(index_type value) {
        if (buffer.size() - used < 16) flush();
        auto result = std::to_chars(buffer.data() + used, buffer.data() + buffer.size(), value);
        used = static_cast<size_t>(result.ptr - buffer.data());
        buffer[used++] = '\n';
    }
    void finish() {
        flush();
        out.close();
        if (!out) throw std::runtime_error("failed to close output");
    }
};

} // namespace lamas

#ifndef LAMAS_NO_MAIN
int main(int argc, char** argv) {
    try {
        std::string x, y, output;
        int d = -1, opt;
        while ((opt = getopt(argc, argv, "x:y:d:o:h")) != -1) {
            switch (opt) {
                case 'x': x = optarg; break;
                case 'y': y = optarg; break;
                case 'o': output = optarg; break;
                case 'd': {
                    const auto end = optarg + std::strlen(optarg);
                    const auto parsed = std::from_chars(optarg, end, d);
                    if (parsed.ec != std::errc{} || parsed.ptr != end || d < 0)
                        throw std::invalid_argument("-d must be an integer in [0, INT_MAX]");
                    break;
                }
                case 'h':
                    std::cout << "lamas -x <sequence> -y <sequence> -d <radius> -o <output>\n";
                    return 0;
                default: throw std::invalid_argument("invalid command-line option");
            }
        }
        if (x.empty() || y.empty() || output.empty() || d < 0 || optind != argc)
            throw std::invalid_argument("usage: lamas -x <sequence> -y <sequence> -d <radius> -o <output>");
        const auto start = std::chrono::steady_clock::now();
        std::vector<unsigned char> text;
        lamas::read_append(x, text);
        if (text.size() == static_cast<size_t>(INT_MAX))
            throw std::length_error("combined sequence length exceeds INT_MAX");
        const auto nx = static_cast<lamas::index_type>(text.size());
        text.push_back('\n');
        lamas::read_append(y, text);
        const auto ny = static_cast<lamas::index_type>(text.size() - nx - 1);
        lamas::Output out(output);
        auto emit = [&](lamas::index_type value) { out(value); };
        if (ny == 0) { /* Empty output. */ }
        else if (nx == 0) { for (int i = 0; i < ny; ++i) emit(lamas::no_match); }
        else if (d == 0) lamas::solve_zero(text, nx, ny, emit);
        else if (d <= lamas::direct_radius) lamas::solve_small(text, nx, ny, d, emit);
        else {
            lamas::SuffixIndex index(text);
            std::vector<unsigned char>().swap(text); // Release text before block workspace.
            lamas::solve_indexed(nx, ny, d, index, emit);
        }
        out.finish();
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start);
        std::cout << "Elapsed time: " << elapsed.count() << " ms\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
#endif
