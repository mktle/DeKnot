#include "hpc.h"
#include "node.h"
#include "hash.h"
#include "syncmer.h"
#include "utils.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

void syncmer::sketch_sequences(int s, std::vector<std::string>& reads, std::vector<std::string>& names, std::vector<SketchedSeq*>& ret_sketches){
    for (size_t i{0}; i < reads.size(); i++){
        // perform homopolymer compression
        std::string hpc_seq;
        std::vector<int> hp_lengths;
        hpc::get_compressed_sequence(reads[i], hpc_seq, hp_lengths);
        std::string seq_name = names.size() > 0 ? names[i] : "read_" + std::to_string(i);
        SketchedSeq* curr_seq = new SketchedSeq(seq_name, true);
        std::vector<int> indices;

        syncmer::get_all_syncmers(hpc_seq, s, curr_seq->kmers, curr_seq->flipped, curr_seq->unique_seq_len, indices);
        ret_sketches.push_back(curr_seq);
    }
}

SketchedSeq* syncmer::split_sequence(SketchedSeq* init_seq, size_t idx_first, size_t idx_second){
    assert(idx_first > 0);
    assert(idx_second < init_seq->kmers.size() - 2);
    SketchedSeq* new_seq = new SketchedSeq(init_seq->name + "2", init_seq->is_read);
    init_seq->name = init_seq->name + "1";
    new_seq->assignment = init_seq->assignment;
    new_seq->num_tried = init_seq->num_tried;

    // move latter part into new sequence object
    new_seq->kmers.assign(init_seq->kmers.begin() + (long)idx_second, init_seq->kmers.end());
    new_seq->flipped.assign(init_seq->flipped.begin() + (long)idx_second, init_seq->flipped.end());
    new_seq->unique_seq_len.assign(init_seq->unique_seq_len.begin() + (long)idx_second, init_seq->unique_seq_len.end());

    // shrink original to only keep first part
    init_seq->kmers.resize(idx_first + 1);
    init_seq->flipped.resize(idx_first + 1);
    init_seq->unique_seq_len.resize(idx_first + 1);
    return new_seq;
}

void syncmer::get_all_syncmers(std::string& sequence, int s, std::vector<uint64_t>& kmer_list, std::vector<uint8_t>& flipped_list, std::vector<int>& unique_seq_len, std::vector<int>& indices){
    if (sequence.size() < (size_t)KMER_SIZE) return;
    const int k {KMER_SIZE};
    uint64_t curr_hash {0};
    uint64_t rc_hash {0};
    constexpr size_t RING = 64;
    constexpr size_t RMASK = RING - 1;
    uint64_t fv[RING];
    int fp[RING];
    uint64_t rv[RING];
    int rp[RING];
    size_t f_head {0};
    size_t f_tail {0};
    size_t r_head {0};
    size_t r_tail {0};
    int prev_pos {0};
    const size_t n = sequence.size();
    const char* seq = sequence.data();
    for (int pos{0}; pos < (int)n; pos++){
        curr_hash += hash::base_to_val(seq[(size_t)pos]);
        curr_hash = curr_hash << 2;
        rc_hash = rc_hash >> 2;
        rc_hash += hash::base_to_rc_val(seq[(size_t)pos]) << (2 * KMER_SIZE);
        rc_hash &= ~(uint64_t)3;
        if (pos + 1 < s) continue;
        int q = pos + 1 - s;
        uint64_t fval = hash::get_smer_hash(s, k - s, curr_hash);
        uint64_t rval = hash::get_smer_hash(s, 0, rc_hash);
        while (f_tail > f_head && fv[(f_tail - 1) & RMASK] > fval) f_tail--;
        fv[f_tail & RMASK] = fval;
        fp[f_tail & RMASK] = q;
        f_tail++;
        while (r_tail > r_head && rv[(r_tail - 1) & RMASK] >= rval) r_tail--;
        rv[r_tail & RMASK] = rval;
        rp[r_tail & RMASK] = q;
        r_tail++;
        if (pos + 1 < k) continue;
        int p = pos + 1 - k;
        while (fp[f_head & RMASK] < p) f_head++;
        while (rp[r_head & RMASK] < p) r_head++;
        bool fwd_sync = (fp[f_head & RMASK] == p) || (fp[f_head & RMASK] == p + k - s);
        bool rc_sync = (rp[r_head & RMASK] == p) || (rp[r_head & RMASK] == p + k - s);
        if (fwd_sync || rc_sync){
            bool took_reverse = rc_hash < curr_hash;
            kmer_list.push_back(took_reverse ? rc_hash : curr_hash);
            flipped_list.push_back(took_reverse);
            unique_seq_len.push_back(pos + 1 - prev_pos);
            prev_pos = pos + 1;
            indices.push_back(p);
        }
    }
}
bool syncmer::next_syncmer(int k, int s, std::string& seq, int& pos_next, uint64_t&hash, uint64_t& rc_hash){
    assert(pos_next > 0);
    if (pos_next == seq.size()){
        pos_next = -1;
        return false;
    }
    bool took_reverse = hash::hash_next(hash, rc_hash, seq.at((size_t)pos_next));
    pos_next++;
    while (!(hash_is_syncmer(k, s, hash) || hash_is_syncmer(k, s, rc_hash))){
        if (pos_next == seq.size()){
            pos_next = -1;
            return false;
        }
        took_reverse = hash::hash_next(hash, rc_hash, seq.at((size_t)pos_next));
        pos_next++;
    }
    return took_reverse;
}

// TODO: can roll for consecutive syncmers so we don't have to keep checking every single position
bool syncmer::hash_is_syncmer(int k, int s, uint64_t hash){
    int smallest_pos {0};
    uint64_t smallest_hash {hash::get_smer_hash(s, 0, hash)};
    for(int i{1}; i < k - s + 1; i++){
        uint64_t next_hash {hash::get_smer_hash(s, i, hash)};
        if (next_hash < smallest_hash){
            smallest_pos = i;
            smallest_hash = next_hash;
        }
    }
    return (smallest_pos == 0) || (smallest_pos == k - s);
}

// returns whether the canonical k-min-mer is the reverse of the current representation
