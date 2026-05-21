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
    if(init_seq->name == "m84039_230117_233243_s1/56365526/ccsM") std::cout << "splitting" << std::endl;//mkl
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
    uint64_t curr_hash {0};
    uint64_t rc_hash {0};
    bool first_is_syncmer {false};
    std::string start_kmer = sequence.substr(0, KMER_SIZE);
    bool took_reverse = hash::hash_entire_kmer(curr_hash, rc_hash, start_kmer, s, first_is_syncmer);
    int next_pos {KMER_SIZE};
    if (!first_is_syncmer){
        took_reverse = syncmer::next_syncmer(KMER_SIZE, s, sequence, next_pos, curr_hash, rc_hash);
    }
    int prev_pos {0};
    //  index will be set to -1 when the hash function gets to the end of the sequence
    while (next_pos != -1){
        uint64_t hash = (took_reverse) ? rc_hash : curr_hash;
        kmer_list.push_back(hash);
        flipped_list.push_back(took_reverse);
        unique_seq_len.push_back(next_pos - prev_pos);
        prev_pos = next_pos;
        indices.push_back(next_pos - KMER_SIZE);
        took_reverse = syncmer::next_syncmer(KMER_SIZE, s, sequence, next_pos, curr_hash, rc_hash);
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
char syncmer::canonical_syncmer(Span<uint64_t>& to_check){
    assert(to_check.size() > 1);
    size_t start {0};
    size_t end {to_check.size() - 1};
    // break ties by checking next kmers
    while (to_check[start] == to_check[end] && (start < end - 1)){
        start++;
        end--;
    }
    if (to_check[start] >= to_check[end]){
        return REVERSE_ORIENTATION;
    }
    return FORWARD_ORIENTATION;
}
