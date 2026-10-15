#ifndef HASH_H
#define HASH_H

#include "hpc.h"
#include "syncmer.h"

#include <ankerl/unordered_dense.h>
#include <stdint.h>
#include <string>

// to be filled in by initialize_dinucleotide_hashes()
extern ankerl::unordered_dense::set<uint64_t> DINUCLEOTIDE_SUBHASHES;

namespace hash{
    bool hash_entire_kmer(uint64_t& curr_hash, uint64_t& rc_hash, std::string& kmer, int s, bool& is_syncmer);
    bool hash_next(uint64_t& curr_hash, uint64_t& rc_hash, char next_char);
    extern const uint64_t* const BASE_VAL;
    extern const uint64_t* const BASE_RC_VAL;
    inline uint64_t get_smer_hash(int s, int pos, uint64_t hash){
        // shift relevant smer to the least significant bits, then mask
        return (hash >> static_cast<uint64_t>(2 * (32 - pos - s))) & ((uint64_t(1) << static_cast<uint64_t>(2 * s)) - 1);
    }
    inline uint64_t base_to_val(char base){
        return BASE_VAL[(unsigned char)base];
    }
    inline uint64_t base_to_rc_val(char base){
        return BASE_RC_VAL[(unsigned char)base];
    }
    void hash_to_kmer(uint64_t hash, bool flipped, std::string& ret_string);
    void add_to_seq_from_hash(uint64_t hash, bool flipped, int len, std::string& ret_string);
    void prepend_to_seq_from_hash(uint64_t hash, bool flipped, int len, std::string& ret_string);
    uint64_t encode_subkmer(std::string& kmer);
    void initialize_dinucleotide_hashes();
}

#endif
