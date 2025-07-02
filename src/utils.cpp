#include "hpc.h"
#include "utils.h"
#include "xxhash.h"

#include <cassert>

size_t utils::hash_kmers(Span<uint64_t>& kmers, bool flip) {
    if (flip){
        std::vector<uint64_t> reversed(kmers.rbegin(), kmers.rend());
        return XXH3_64bits(reversed.data(), reversed.size() * sizeof(uint64_t));
    }
    return XXH3_64bits(kmers.data(), kmers.size() * sizeof(uint64_t));
}

char utils::reverse_orientation(char orig){
    assert(orig == FORWARD_ORIENTATION || orig == REVERSE_ORIENTATION);
    if (orig == FORWARD_ORIENTATION){
        return REVERSE_ORIENTATION;
    }
    return FORWARD_ORIENTATION;
}
