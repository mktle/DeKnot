#include "hpc.h"
#include "utils.h"

#include <cassert>

const uint64_t utils::HASH_BASE = 0x9E3779B97F4A7C15ULL;

static uint64_t mod_inverse_2_64(uint64_t a){
    uint64_t x = a;
    for (int i = 0; i < 5; i++){
        x *= 2 - a * x;
    }
    return x;
}
const uint64_t utils::HASH_BASE_INV = mod_inverse_2_64(utils::HASH_BASE);

uint64_t utils::pow_base(size_t e){
    uint64_t r = 1;
    uint64_t b = utils::HASH_BASE;
    while(e){
        if (e & 1){
            r *= b;
        }
        b *= b;
        e >>= 1;
    }
    return r;
}

size_t utils::hash_kmers(Span<uint64_t>& kmers, bool flip){
    uint64_t h = 0;
    if (!flip){
        for (size_t i = 0; i < kmers.size(); i++){
            h = h * utils::HASH_BASE + kmers[i];
        }
    }else{
        for (size_t i = kmers.size(); i-- > 0;){
            h = h * utils::HASH_BASE + kmers[i];
        }
    }
    return (size_t)h;
}

char utils::reverse_orientation(char orig){
    assert(orig == FORWARD_ORIENTATION || orig == REVERSE_ORIENTATION);
    if (orig == FORWARD_ORIENTATION){
        return REVERSE_ORIENTATION;
    }
    return FORWARD_ORIENTATION;
}
