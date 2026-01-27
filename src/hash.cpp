#include "hash.h"
#include "syncmer.h"

#include <cassert>
#include <cmath>
#include <cstring>

ankerl::unordered_dense::set<uint64_t> DINUCLEOTIDE_SUBHASHES{};
const char val_to_base[] = "ACGT";

bool hash::hash_entire_kmer(uint64_t& curr_hash, uint64_t& rc_hash, std::string& kmer, int s, bool& is_syncmer){
    assert(kmer.size() == KMER_SIZE);
    curr_hash = 0;
    rc_hash = 0;

    for (size_t i{0}; i < kmer.size(); i++){
        curr_hash += base_to_val(kmer.at(i));
        curr_hash = curr_hash << 2;
        rc_hash = rc_hash >> 2;
        rc_hash += base_to_rc_val(kmer.at(i)) << (2 * KMER_SIZE);
    }

    // check if the forward hash is a syncmer
    is_syncmer = syncmer::hash_is_syncmer((int)kmer.size(), s, curr_hash) || syncmer::hash_is_syncmer((int)kmer.size(), s, rc_hash);

    // shouldn't be equal since we're using an odd k
    assert(rc_hash != curr_hash);
    if (rc_hash < curr_hash){
        return true;
    }
    return false;
}

// returns true if we take the reverse complement as canonical
bool hash::hash_next(uint64_t& curr_hash, uint64_t& rc_hash, char next_char){
    curr_hash += base_to_val(next_char);
    curr_hash = curr_hash << 2;

    rc_hash = rc_hash >> 2;
    rc_hash += base_to_rc_val(next_char) << (2*KMER_SIZE);

    // remove the last two bits / set them to padding
    rc_hash &= ~(uint64_t)3;

    // shouldn't be equal since we're using an odd k
    assert(rc_hash != curr_hash);
    if (rc_hash < curr_hash){
        return true;
    }
    return false;
}

uint64_t hash::get_smer_hash(int s, int pos, uint64_t hash){
    // shift relevant smer to the least significant bits, then mask
    return (hash >> static_cast<uint64_t>(2 * (32 - pos - s))) & ((uint64_t(1) << static_cast<uint64_t>(2 * s)) - 1);
}

uint64_t hash::base_to_val(char base){
    if (base == 'A') return 0;
    if (base == 'C') return 1;
    if (base == 'G') return 2;
    if (base == 'T') return 3;
    assert(0);
    return 4;
}

uint64_t hash::base_to_rc_val(char base){
    if (base == 'A') return 3;
    if (base == 'C') return 2;
    if (base == 'G') return 1;
    if (base == 'T') return 0;
    assert(0);
    return 4;
}

void hash::hash_to_kmer(uint64_t hash, bool flipped, std::string& ret_string){
    ret_string = "";
    if (flipped){
        for (size_t i{0}; i < KMER_SIZE; i++){
            hash = hash >> 2;
            uint64_t rc_val = ~hash;
            rc_val = rc_val & 3;
            ret_string += val_to_base[rc_val];
        }
    }else{
        for (size_t i{0}; i < KMER_SIZE; i++){
            uint64_t curr_val = hash >> ( 2 * (KMER_SIZE - i));
            curr_val &= 3;
            ret_string += val_to_base[curr_val];
        }
    }
}

void hash::add_to_seq_from_hash(uint64_t hash, bool flipped, int len, std::string& ret_string){
    if (flipped){
        for (size_t i{0}; i < KMER_SIZE; i++){
            hash = hash >> 2;
            if (i >= KMER_SIZE - len){
                uint64_t rc_val = ~hash;
                rc_val = rc_val & 3;
                ret_string += val_to_base[rc_val];
            }
        }
    }else{
        for (int i{KMER_SIZE - len}; i < KMER_SIZE; i++){
            uint64_t curr_val = hash >> ( 2 * (KMER_SIZE - i));
            curr_val &= 3;
            ret_string += val_to_base[curr_val];
        }
    }
}

void hash::prepend_to_seq_from_hash(uint64_t hash, bool flipped, int len, std::string& ret_string){
    std::string to_prepend;
    if (flipped){
        for (int i{0}; i < len; i++){
            hash = hash >> 2;
            uint64_t rc_val = ~hash;
            rc_val = rc_val & 3;
            to_prepend += val_to_base[rc_val];
        }
    }else{
        for (int i{0}; i < len; i++){
            uint64_t curr_val = hash >> ( 2 * (KMER_SIZE - i));
            curr_val &= 3;
            to_prepend += val_to_base[curr_val];
        }
    }
    ret_string = to_prepend + ret_string;
}

uint64_t hash::encode_subkmer(std::string& kmer) {
    uint64_t hash = 0;
    for (char c : kmer) {
        hash = (hash << 2) | base_to_val(c);
    }
    return hash;
}

void hash::initialize_dinucleotide_hashes(){
    std::vector<std::string> to_hash = {"ACACAC", "AGAGAG", "ATATAT", "CACACA", "CGCGCG", "CTCTCT", "GAGAGA", "GCGCGC", "GTGTGT", "TATATA", "TCTCTC", "TGTGTG" };
    for (size_t i{0}; i < to_hash.size(); i++){
        DINUCLEOTIDE_SUBHASHES.insert(hash::encode_subkmer(to_hash[i]));
    }

}
