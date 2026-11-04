#ifndef HPC_H
#define HPC_H

#include "argument_parser.h"
#include "syncmer.h"

#include <array>
#include <string>
#include <vector>

struct HPCCounts{
    int* coverages = NULL;
    int* sums = NULL;
    int len;

    HPCCounts(int size){
        len = size;
        if (size > 0){
            coverages = new int[(size_t)size]{};
            sums = new int[(size_t)size]{};
        }
    }
};

struct SeqIdx{
    // we are aligning to a list of sequences; this idx is the index into the list
    int seq_num;
    int seq_pos;
    bool reversed;
    // if we select all the syncmers from the sequence, this is the position in that ordered list
    int syncmer_idx;

    SeqIdx() : seq_num(0), seq_pos(0), reversed(false), syncmer_idx(0) {}
    SeqIdx(int n, int p, bool r, int s) : seq_num(n), seq_pos(p), reversed(r), syncmer_idx(s) {}
};

struct ChainHit{
    int utg_num;
    int start;
    int len;
    int read_start;
    bool forward;
    ChainHit(int u, int s, int l, int r, bool f) : utg_num(u), start(s), len(l), read_start(r), forward(f) {}
};

namespace hpc{
    void get_compressed_sequence(std::string& orig, std::string& ret, std::vector<int>& hp_lengths);
    void align_single_read(std::vector<uint64_t>& kmer_list, ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>>& kmer_index, std::vector<SketchedSeq>& unitig_sketches, std::vector<int>& hp_lengths, std::vector<HPCCounts>& hpc_tallies, std::vector<int>& kmer_indices, std::vector<uint8_t>& flipped_list);
    void collect_all_chains(std::vector<uint64_t>& kmer_list, ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>>& kmer_index, std::vector<SketchedSeq>& unitig_sketches, std::vector<int>& kmer_indices, std::vector<uint8_t>& flipped_list, std::vector<ChainHit>& hits);
    size_t check_chain_direction(bool forward, std::vector<SketchedSeq>& unitig_sketches, SeqIdx& first_hit, std::vector<uint64_t>& kmer_list, size_t start_kmer_idx, std::vector<uint8_t>& flipped_list);
    void add_to_hpc_tally(std::vector<int>& hp_lengths, HPCCounts& hpc_tally, int utg_pos, int read_pos, int len, bool forward);
}

#endif
