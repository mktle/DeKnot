#include "hpc.h"
#include "interval_tree.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

// one placement of read segment on a unitig
struct SegmentTarget{
    int utg_num;
    int utg_pos;
    bool forward;

    SegmentTarget(int u_n, int u_p, bool f) : utg_num(u_n), utg_pos(u_p), forward(f) {}
};

struct AlignedSegment{
    int length;
    int read_pos;
    std::vector<SegmentTarget> targets;

    AlignedSegment(int l, int r_p) : length(l), read_pos(r_p), targets() {}
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

static void find_alignments(std::vector<uint64_t>& kmer_list, ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>>& kmer_index, std::vector<SketchedSeq>& unitig_sketches, std::vector<int>& kmer_indices, std::vector<uint8_t>& flipped_list, std::vector<AlignedSegment>& alignments){
    /*********************************************************************************************
     * Find all segments with exact match to a unitig, seeking the unitig with the longest match *
     *********************************************************************************************/
    if (kmer_list.size() < 2) return;
    // check segments for hits
    size_t i{0};
    while (i < kmer_list.size() - 1){
        uint64_t curr_kmer = kmer_list[i];
        // check the matches to the utg sequences
        if (kmer_index[curr_kmer].size() > 0){
            AlignedSegment longest_segment(-1, kmer_indices[i]);
            size_t longest_match_len {1};
            // check each hit to see if there is a large segment that matches
            // we can assume this kmer is the start of the match or else it would have already been added with a previous segment
            for (size_t j{0}; j < kmer_index[curr_kmer].size(); j++){
                SeqIdx hit = kmer_index[curr_kmer][j];

                // if in a repetitive region, may have both a hit both in front and behind of kmer
                // so need to check both directions...
                size_t num_forward = hpc::check_chain_direction(true, unitig_sketches, hit, kmer_list, i, flipped_list);
                size_t num_backward = hpc::check_chain_direction(false,  unitig_sketches, hit, kmer_list, i, flipped_list);

                bool use_forward = num_forward > num_backward;
                size_t match_len = use_forward ? num_forward : num_backward;

                if (match_len > longest_match_len){
                    int length = kmer_indices[(size_t)i + match_len - 1] + KMER_SIZE - kmer_indices[(size_t)i];
                    assert(length > 0);
                    longest_match_len = match_len;
                    longest_segment.length = length;
                    longest_segment.targets.clear();
                    longest_segment.targets.emplace_back(hit.seq_num, hit.seq_pos, use_forward);
                }else if (match_len == longest_match_len && !longest_segment.targets.empty()){
                    // read matches equally well at multiple places, so store all
                    longest_segment.targets.emplace_back(hit.seq_num, hit.seq_pos, use_forward);
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
}

void hpc::align_single_read(std::vector<uint64_t>& kmer_list, ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>>& kmer_index, std::vector<SketchedSeq>& unitig_sketches, std::vector<int>& hp_lengths, std::vector<HPCCounts>& hpc_tallies, std::vector<int>& kmer_indices, std::vector<uint8_t>& flipped_list){
    std::vector<AlignedSegment> alignments;
    find_alignments(kmer_list, kmer_index, unitig_sketches, kmer_indices, flipped_list, alignments);
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
                for (const SegmentTarget& target : alignments[i].targets){
                    hpc::add_to_hpc_tally(hp_lengths, hpc_tallies[(size_t)target.utg_num], target.utg_pos, alignments[i].read_pos, alignments[i].length, target.forward);
                }
                IntervalTree::insert(query_intervals, read_start, read_end);
            }
        }
    }
    IntervalTree::destroy_tree(query_intervals);
}

void hpc::collect_all_chains(std::vector<uint64_t>& kmer_list, ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>>& kmer_index, std::vector<SketchedSeq>& unitig_sketches, std::vector<int>& kmer_indices, std::vector<uint8_t>& flipped_list, std::vector<ChainHit>& hits){
    if (kmer_list.size() < 2) return;
    std::vector<ChainHit> raw;
    ankerl::unordered_dense::map<int, size_t> skip_until;
    for (size_t i{0}; i + 1 < kmer_list.size(); i++){
        auto it = kmer_index.find(kmer_list[i]);
        if (it == kmer_index.end()) continue;
        for (size_t j{0}; j < it->second.size(); j++){
            SeqIdx hit = it->second[j];
            auto sk = skip_until.find(hit.seq_num);
            if (sk != skip_until.end() && i < sk->second) continue;
            size_t num_forward = hpc::check_chain_direction(true, unitig_sketches, hit, kmer_list, i, flipped_list);
            size_t num_backward = hpc::check_chain_direction(false, unitig_sketches, hit, kmer_list, i, flipped_list);
            bool use_forward = num_forward > num_backward;
            size_t match_len = use_forward ? num_forward : num_backward;
            if (match_len <= 30) continue;
            int len = kmer_indices[i + match_len - 1] + KMER_SIZE - kmer_indices[i];
            if (len < 150) continue;
            int start = use_forward ? hit.seq_pos : hit.seq_pos + KMER_SIZE - len;
            raw.emplace_back(hit.seq_num, start, len, kmer_indices[i], use_forward);
            size_t until = i + match_len;
            if (sk == skip_until.end()) skip_until[hit.seq_num] = until;
            else if (until > sk->second) sk->second = until;
        }
    }
    const int MAX_GAP = 100;
    const int MAX_DIAG_DIFF = 20;
    std::sort(raw.begin(), raw.end(), [](const ChainHit& x, const ChainHit& y){
        if (x.utg_num != y.utg_num) return x.utg_num < y.utg_num;
        if (x.forward != y.forward) return x.forward;
        return x.start < y.start;
    });
    for (size_t i{0}; i < raw.size(); i++){
        ChainHit cur = raw[i];
        while (i + 1 < raw.size() && raw[i + 1].utg_num == cur.utg_num && raw[i + 1].forward == cur.forward){
            const ChainHit& nxt = raw[i + 1];
            int gap_t = nxt.start - (cur.start + cur.len);
            int gap_r = cur.forward ? nxt.read_start - (cur.read_start + cur.len) : cur.read_start - (nxt.read_start + nxt.len);
            if (gap_t < -MAX_DIAG_DIFF || gap_t > MAX_GAP || gap_r < -MAX_DIAG_DIFF || gap_r > MAX_GAP || std::abs(gap_t - gap_r) > MAX_DIAG_DIFF) break;
            int new_end = nxt.start + nxt.len;
            if (!cur.forward) cur.read_start = nxt.read_start;
            cur.len = new_end - cur.start;
            i++;
        }
        hits.push_back(cur);
    }
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
