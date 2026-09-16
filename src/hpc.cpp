#include "hpc.h"
#include "interval_tree.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

struct AlignedSegment{
    int length;
    int utg_num;
    int utg_pos;
    int read_pos;
    bool forward;

    AlignedSegment(int l, int u_n, int u_p, int r_p, bool f) : length(l), utg_num(u_n), utg_pos(u_p), read_pos(r_p), forward(f) {}
};

bool compareAlignedSegment(const AlignedSegment& a, const AlignedSegment& b){
    return a.length > b.length;
}

void hpc::get_compressed_sequence(std::string& orig, std::string& ret, std::vector<int>& hp_lengths){
    ret.reserve(orig.length());

    int block_start {0};
    ret += orig[0];
    int i{1};
    for (; i < orig.length(); i++){
        if (orig[(size_t)i] != orig[(size_t)block_start]) {
            ret += orig[(size_t)i];
            hp_lengths.push_back(i - block_start);
            block_start = i;
        }
    }
    // need to add the length of the last block
    hp_lengths.push_back(i - block_start);
    assert(hp_lengths.size() == ret.length());
}

void hpc::align_single_read(std::vector<uint64_t>& kmer_list, ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>>& kmer_index, std::vector<SketchedSeq>& unitig_sketches, std::vector<int>& hp_lengths, std::vector<HPCCounts>& hpc_tallies, std::vector<int>& kmer_indices, std::vector<uint8_t>& flipped_list){
    /*********************************************************************************************
     * Find all segments with exact match to a unitig, seeking the unitig with the longest match *
     *********************************************************************************************/
    std::vector<AlignedSegment> alignments;
    if (kmer_list.size() < 2) return;
    // check segments for hits
    size_t i{0};
    while (i < kmer_list.size() - 1){
        uint64_t curr_kmer = kmer_list[i];
        // check the matches to the utg sequences
        if (kmer_index[curr_kmer].size() > 0){
            AlignedSegment longest_segment(-1, -1, -1, -1, false);
            size_t longest_match_len {1};
            // check each hit to see if there is a large segment that matches
            // we can assume this kmer is the start of the match or else it would have already been added with a previous segment
            for (size_t j{0}; j < kmer_index[curr_kmer].size(); j++){
                SeqIdx hit = kmer_index[curr_kmer][j];

                // if in a repetitive region, may have both a hit both in front and behind of kmer
                // so need to check both directions...
                size_t num_forward = check_chain_direction(true, unitig_sketches, hit, kmer_list, i, flipped_list);
                size_t num_backward = check_chain_direction(false,  unitig_sketches, hit, kmer_list, i, flipped_list);

                bool use_forward = num_forward > num_backward;
                size_t match_len = use_forward ? num_forward : num_backward;

                if (match_len > longest_match_len){
                    int length = kmer_indices[(size_t)i + match_len - 1] + KMER_SIZE - kmer_indices[(size_t)i];
                    assert(length > 0);
                    longest_match_len = match_len;
                    longest_segment = AlignedSegment(length, hit.seq_num, hit.seq_pos, kmer_indices[i], use_forward);
                }
            }
            if (longest_match_len > 30){
                alignments.push_back(std::move(longest_segment));
                // skip to next non-overlapping kmer
                size_t last_kmer = i + longest_match_len - 1;
                size_t check_overlap = last_kmer + 1;
                while(check_overlap < kmer_indices.size() && kmer_indices[check_overlap] < kmer_indices[last_kmer] + KMER_SIZE){
                    check_overlap++;
                    longest_match_len++;
                }
            }
            i += longest_match_len;
        }else{
            i++;
        }
    }
    /*************************************************************************
     * Only add matches to tally for the unitig with the longest exact match *
     *************************************************************************/
    IntervalNode* query_intervals = NULL;
    if (alignments.size() > 0){
        // sort by aligned length
        std::sort(alignments.begin(), alignments.end(), compareAlignedSegment);
        for (size_t i{0}; i < alignments.size(); i++){
            if (alignments[i].length < 150) break;
            int read_start = alignments[i].read_pos;
            int read_end = read_start + alignments[i].length - 1;
            if (query_intervals == NULL || !IntervalTree::overlap_search(query_intervals, read_start, read_end)){
                hpc::add_to_hpc_tally(hp_lengths, hpc_tallies[(size_t)alignments[i].utg_num], alignments[i].utg_pos, alignments[i].read_pos, alignments[i].length, alignments[i].forward);
                IntervalTree::insert(query_intervals, read_start, read_end);
            }
        }
    }
    IntervalTree::destroy_tree(query_intervals);
}

// return number of matching kmers
size_t hpc::check_chain_direction(bool forward, std::vector<SketchedSeq>& unitig_sketches, SeqIdx& first_hit, std::vector<uint64_t>& kmer_list, size_t start_kmer_idx, std::vector<uint8_t>& flipped_list){
    bool init_dir_matches {flipped_list[start_kmer_idx + 1] == first_hit.reversed};
    size_t num_hits {1};
    // advancing forward on the read
    while(start_kmer_idx + num_hits < kmer_list.size()){
        // next kmer on this read
        uint64_t next_kmer = kmer_list[start_kmer_idx + num_hits];
        // calculate which kmer idx we need depending on direction we're checking
        int next_utg_idx = forward ? first_hit.syncmer_idx + (int)num_hits : first_hit.syncmer_idx - (int)num_hits;
        if (next_utg_idx <= -1 || next_utg_idx >= (int)unitig_sketches[(size_t)first_hit.seq_num].kmers.size()){
            // we reached the end of the unitig sequence
            break;
        }
        uint64_t next_utg_kmer = unitig_sketches[(size_t)first_hit.seq_num].kmers[(size_t)next_utg_idx];
        bool next_dir_matches {flipped_list[start_kmer_idx + num_hits] == unitig_sketches[(size_t)first_hit.seq_num].flipped[(size_t)next_utg_idx]};

        if (!(next_kmer == next_utg_kmer && init_dir_matches == next_dir_matches)){
            break;
        }
        num_hits++;
    }
    return num_hits;
}

void hpc::add_to_hpc_tally(std::vector<int>& hp_lengths, HPCCounts& hpc_tally, int utg_pos, int read_pos, int len, bool forward){
    if (!forward){
        utg_pos += KMER_SIZE - 1;
    }
    for (int i{0}; i < len; i++){
        int utg_idx = forward ? utg_pos + i : utg_pos - i;
        hpc_tally.coverages[utg_idx]++;
        hpc_tally.sums[utg_idx] += hp_lengths[(size_t)(read_pos + i)];
    }
}
