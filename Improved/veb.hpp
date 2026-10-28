#ifndef LAMAS_VEB_HPP
#define LAMAS_VEB_HPP
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <vector>

// Eager van Emde Boas tree. Universe is a power of two, rounded up by caller.
// Strict predecessor/successor; -1 means absent. insert/erase require absent/
// present keys respectively. All storage is allocated once and reused.
class VEB {
    unsigned bits_, low_;
    int64_t min_ = -1, max_ = -1;
    std::unique_ptr<VEB> summary_;
    std::vector<std::unique_ptr<VEB>> clusters_;
    int64_t high(int64_t x) const { return x >> low_; }
    int64_t low(int64_t x) const { return x & ((int64_t(1) << low_) - 1); }
    int64_t join(int64_t h, int64_t l) const { return (h << low_) | l; }
public:
    explicit VEB(unsigned bits) : bits_(bits), low_(bits / 2) {
        assert(bits >= 1 && bits <= 34);
        if (bits > 1) {
            summary_.reset(new VEB(bits - low_));
            const size_t count = size_t(1) << (bits - low_);
            clusters_.reserve(count);
            for (size_t i = 0; i < count; ++i)
                clusters_.emplace_back(new VEB(low_));
        }
    }
    int64_t minimum() const { return min_; }
    int64_t maximum() const { return max_; }
    void clear() {
        min_ = max_ = -1;
        if (summary_) {
            summary_->clear();
            for (auto &c : clusters_) c->clear();
        }
    }
    // O(U) bulk construction from a U-byte membership array. Recurse bottom-up,
    // then exclude the minimum from its cluster, as required by the vEB invariant.
    // Per-node O(sqrt(U) + log log U) overhead sums to O(U) over the eager tree.
    void build(const unsigned char *present) {
        min_ = max_ = -1;
        if (bits_ == 1) {
            if (present[0]) min_ = max_ = 0;
            if (present[1]) { if (min_ < 0) min_ = 1; max_ = 1; }
            return;
        }
        std::vector<unsigned char> occupied(clusters_.size(), 0);
        for (size_t h = 0; h < clusters_.size(); ++h) {
            auto &c = *clusters_[h];
            c.build(present + (h << low_));
            if (c.min_ >= 0) {
                if (min_ < 0) min_ = join(h, c.min_);
                max_ = join(h, c.max_);
            }
        }
        if (min_ >= 0) clusters_[high(min_)]->erase(low(min_));
        for (size_t h = 0; h < clusters_.size(); ++h)
            occupied[h] = clusters_[h]->min_ >= 0;
        summary_->build(occupied.data());
    }
    void insert(int64_t x) {
        if (min_ < 0) { min_ = max_ = x; return; }
        if (x < min_) std::swap(x, min_);
        if (bits_ > 1) {
            auto &c = *clusters_[high(x)];
            if (c.min_ < 0) {
                summary_->insert(high(x));
                c.min_ = c.max_ = low(x);
            } else c.insert(low(x));
        }
        if (x > max_) max_ = x;
    }
    void erase(int64_t x) {
        if (min_ == max_) { min_ = max_ = -1; return; }
        if (bits_ == 1) { min_ = max_ = 1 - x; return; }
        if (x == min_) {
            int64_t h = summary_->min_;
            x = min_ = join(h, clusters_[h]->min_);
        }
        auto &c = *clusters_[high(x)];
        c.erase(low(x));
        if (c.min_ < 0) {
            summary_->erase(high(x));
            if (x == max_) {
                int64_t h = summary_->max_;
                max_ = h < 0 ? min_ : join(h, clusters_[h]->max_);
            }
        } else if (x == max_) max_ = join(high(x), c.max_);
    }
    int64_t successor(int64_t x) const {
        if (min_ < 0 || x >= max_) return -1;
        if (x < min_) return min_;
        if (bits_ == 1) return max_;
        const auto &c = *clusters_[high(x)];
        if (c.max_ >= 0 && low(x) < c.max_)
            return join(high(x), c.successor(low(x)));
        int64_t h = summary_->successor(high(x));
        return h < 0 ? -1 : join(h, clusters_[h]->min_);
    }
    int64_t predecessor(int64_t x) const {
        if (min_ < 0 || x <= min_) return -1;
        if (x > max_) return max_;
        if (bits_ == 1) return min_;
        const auto &c = *clusters_[high(x)];
        if (c.min_ >= 0 && low(x) > c.min_)
            return join(high(x), c.predecessor(low(x)));
        int64_t h = summary_->predecessor(high(x));
        return h < 0 ? min_ : join(h, clusters_[h]->max_);
    }
};
#endif
