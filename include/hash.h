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
    uint64_t get_smer_hash(int s, int pos, uint64_t hash);
    uint64_t base_to_val(char base);
    uint64_t base_to_rc_val(char base);
    void hash_to_kmer(uint64_t hash, bool flipped, std::string& ret_string);
    void add_to_seq_from_hash(uint64_t hash, bool flipped, int len, std::string& ret_string);
    void prepend_to_seq_from_hash(uint64_t hash, bool flipped, int len, std::string& ret_string);
    uint64_t encode_subkmer(std::string& kmer);
    void initialize_dinucleotide_hashes();
}

#endif
