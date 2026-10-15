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

static uint64_t make_val_table(uint64_t a, uint64_t c, uint64_t g, uint64_t t, uint64_t* table){
    for (int i{0}; i < 256; i++) table[i] = 4;
    table[(unsigned char)'A'] = a;
    table[(unsigned char)'C'] = c;
    table[(unsigned char)'G'] = g;
    table[(unsigned char)'T'] = t;
    return 0;
}
static uint64_t base_val_storage[256];
static uint64_t base_rc_val_storage[256];
static const uint64_t base_tables_init = make_val_table(0, 1, 2, 3, base_val_storage) + make_val_table(3, 2, 1, 0, base_rc_val_storage);
const uint64_t* const hash::BASE_VAL = base_val_storage;
const uint64_t* const hash::BASE_RC_VAL = base_rc_val_storage;

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
    char buf[64];
    int n {0};
    if (flipped){
        for (int i{KMER_SIZE - len}; i < KMER_SIZE; i++){
            uint64_t rc_val = (~(hash >> (2 * (i + 1)))) & 3;
            buf[n++] = val_to_base[rc_val];
        }
    }else{
        for (int i{KMER_SIZE - len}; i < KMER_SIZE; i++){
            uint64_t curr_val = (hash >> (2 * (KMER_SIZE - i))) & 3;
            buf[n++] = val_to_base[curr_val];
        }
    }
    ret_string.append(buf, (size_t)n);
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
