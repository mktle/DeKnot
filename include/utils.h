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
}

struct KmerView{
    Span<uint64_t> kmers;
    bool flip;
    // so we don't need to perform the expensive hash computation every time
    // and we know every KmerView will be hashed at least once
    size_t hash;

    KmerView(Span<uint64_t> k, bool f) : kmers(k), flip(f), hash(utils::hash_kmers(k, f)) {}
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
