#ifndef UTILS_H
#define UTILS_H

#include "span.h"

#include <cstdint>
#include <vector>

const char FORWARD_ORIENTATION = '+';
const char REVERSE_ORIENTATION = '-';
const char INVALID_ORIENTATION = 'X';

namespace utils{
    size_t hash_kmers(Span<uint64_t>& kmers, bool flip);

    extern const uint64_t HASH_BASE;
    extern const uint64_t HASH_BASE_INV;
    uint64_t pow_base(size_t e);

    struct RollingHasher{
        const uint64_t* data;
        size_t n;
        size_t k;
        uint64_t Bk1;
        uint64_t hf;
        uint64_t hr;
        size_t l;

        RollingHasher(const uint64_t* kmers, size_t total, size_t window) : data(kmers), n(total), k(window), Bk1(utils::pow_base(window - 1)), hf(0), hr(0), l(0) {
            for (size_t t = 0; t < k; t++) { hf = hf * utils::HASH_BASE + data[t]; }
            for (size_t t = k; t-- > 0;) { hr = hr * utils::HASH_BASE + data[t]; }
        }

        void advance(){
            uint64_t out = data[l];
            uint64_t in  = data[l + k];
            hf = (hf - out * Bk1) * utils::HASH_BASE + in;
            hr = (hr - out) * utils::HASH_BASE_INV + in * Bk1;
            l++;
        }

        size_t hash_at(size_t idx, bool flip){
            while (l < idx) { advance(); }
            return (size_t)(flip ? hr : hf);
        }
    };
}

struct KmerView{
    Span<uint64_t> kmers;
    bool flip;
    // so we don't need to perform the expensive hash computation every time
    // and we know every KmerView will be hashed at least once
    size_t hash;

    KmerView(Span<uint64_t> k, bool f) : kmers(k), flip(f), hash(utils::hash_kmers(k, f)) {}
    KmerView(Span<uint64_t> k, bool f, size_t precomputed) : kmers(k), flip(f), hash(precomputed) {}
    KmerView() : kmers(Span<uint64_t>()), flip(false), hash(0) {}
};

namespace utils{
    char reverse_orientation(char orig);

    struct KmerViewHash{
        size_t operator()(const KmerView& kmers) const {
            return kmers.hash;
        }
    };

    struct KmerViewEqual {
        bool operator()(const KmerView& a, const KmerView& b) const {
            return a.hash == b.hash;
        }
    };
}

#endif
