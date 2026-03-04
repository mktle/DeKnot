#ifndef SYNCMER_H
#define SYNCMER_H

#include "utils.h"

#include <ankerl/unordered_dense.h>
#include <string>
#include <vector>

const int KMER_SIZE = 31;

struct SketchedSeq{
    // for debugging purposes
    std::string name;
    // as opposed to a unitig
    bool is_read;
    std::vector<uint64_t> kmers;
    std::vector<uint8_t> flipped;
    // track original position on HPC so we can later reconstruct sequence
    std::vector<int> unique_seq_len;
    // track how many times this read has had dropped portions that are not fixed
    int num_tried;
    // for phasing
    int assignment;
    SketchedSeq(std::string n, bool r) : name(n), is_read(r), num_tried(0), assignment(-1) {}
};


namespace syncmer{
    void sketch_sequences(int s, std::vector<std::string>& reads, std::vector<std::string>& names, std::vector<SketchedSeq*>& ret_sketches);
    SketchedSeq* split_sequence(SketchedSeq* init_seq, size_t idx_first, size_t idx_second);
    void get_all_syncmers(std::string& sequence, int s, std::vector<uint64_t>& kmer_list, std::vector<uint8_t>& flipped_list, std::vector<int>& unique_seq_len, std::vector<int>& indices);
    bool next_syncmer(int k, int s, std::string& seq, int& pos_next, uint64_t&hash, uint64_t& rc_hash);
    bool hash_is_syncmer(int k, int s, uint64_t hash);
    char canonical_syncmer(Span<uint64_t>& to_check);
}

#endif
