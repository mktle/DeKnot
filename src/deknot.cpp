#include "deknot.h"
#include "graph.h"
#include "hash.h"
#include "hpc.h"
#include "io.h"
#include "phase.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

const double COLLAPSE_FRAC = 0.25;

const int DEBUG_CHECKPOINTS[] = {99999999};
int DEBUG_IDX = 0;

typedef ankerl::unordered_dense::set<size_t> IdxSet;
static int prepare_iteration(int k_prime, int max_k_prime, std::vector<SketchedSeq*>& input_sequences, Graph* prev_graph, IdxSet& remove_set);
static void build_graph(std::vector<SketchedSeq*>& input_sequences, const IdxSet& base_remove, Graph* prev_graph, Graph* dbg, IdxSet& remove_set);
static void trim_walk_ends(dk_opt_t* opt, std::vector<std::string>& reads, std::vector<std::string>& walks);
static void retire_sequences(std::vector<SketchedSeq*>& input_sequences, const IdxSet& remove_set, std::vector<SketchedSeq*>& prev_seqs);
static bool graph_collapsed(Graph* curr_graph, Graph* prev_graph);

void dk_opt_init(dk_opt_t* opt){
    memset(opt, 0, sizeof(dk_opt_t));
    opt->s = 9;
    opt->init_k_prime = 2;
    opt->max_iterations = 2000;
    opt->compact_thresh = 4;
    opt->return_hpc_utg = 0;
    opt->out_dir = NULL;
    opt->verbose = 0;
    opt->check_strand = 0;
}

int dk_set_opt(const char* preset, dk_opt_t* opt){
    //TODO: establish presets
    assert(preset != NULL);
    assert(opt != NULL);
    return 0;
}

static void syncmers_without_counterpart(Unitig* a, Unitig* b, int s, std::vector<uint64_t>& neutral){
    std::string seq_a, seq_b;
    a->get_hpc_sequence(seq_a, false, false);
    b->get_hpc_sequence(seq_b, false, false);
    if (seq_b.size() > seq_a.size()) std::swap(seq_a, seq_b); // seq_a = longer
    std::vector<uint64_t> ka, kb; std::vector<uint8_t> fa, fb; std::vector<int> ua, ub, ia, ib;
    syncmer::get_all_syncmers(seq_a, s, ka, fa, ua, ia);
    syncmer::get_all_syncmers(seq_b, s, kb, fb, ub, ib);
    ankerl::unordered_dense::map<uint64_t, long> pos_in_a;
    for (size_t i{0}; i < ka.size(); i++){
        auto it = pos_in_a.find(ka[i]);
        if (it == pos_in_a.end()) pos_in_a[ka[i]] = (long)i; else it->second = -1;
    }
    auto median_pos = [&](size_t from, size_t to) -> long {
        std::vector<long> hits;
        for (size_t i{from}; i < to; i++){
            auto it = pos_in_a.find(kb[i]);
            if (it != pos_in_a.end() && it->second >= 0) hits.push_back(it->second);
        }
        if (hits.size() < 30) return -1;
        std::sort(hits.begin(), hits.end());
        return hits[hits.size() / 2];
    };
    size_t n_end = std::min<size_t>(300, kb.size() / 4);
    long p_first = (kb.size() >= 4) ? median_pos(0, n_end) : -1;
    long p_last = (kb.size() >= 4) ? median_pos(kb.size() - n_end, kb.size()) : -1;
    long lo, hi;
    if (p_first >= 0 && p_last >= 0){
        lo = std::min(p_first, p_last); hi = std::max(p_first, p_last);
    }else{
        lo = -1; hi = -1;
        ankerl::unordered_dense::set<uint64_t> in_b(kb.begin(), kb.end());
        for (size_t i{0}; i < ka.size(); i++){
            if (in_b.contains(ka[i])){ if (lo < 0) lo = (long)i; hi = (long)i; }
        }
        if (lo < 0) return;
    }
    long margin = 200;
    for (size_t i{0}; i < ka.size(); i++){
        if ((long)i < lo - margin || (long)i > hi + margin) neutral.push_back(ka[i]);
    }
}

static long env_long(const char* name, long dflt){
    const char* v = getenv(name);
    return (v == NULL) ? dflt : atol(v);
}

static void trim_walk_ends(dk_opt_t* opt, std::vector<std::string>& reads, std::vector<std::string>& walks){
    if (env_long("DK_END_TRIM", 1) == 0) return;
    const long W = env_long("DK_END_TRIM_W", 300);
    const long MIN = env_long("DK_END_TRIM_MIN", 2);
    const long MAXLEN = env_long("DK_END_TRIM_MAXLEN", 20000);
    if (walks.empty()) return;

    ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>> kmer_index;
    std::vector<SketchedSeq> walk_sketches;
    for (size_t i{0}; i < walks.size(); i++){
        SketchedSeq curr_seq = SketchedSeq("walk" + std::to_string(i), false);
        if (walks[i].size() > 0){
            std::vector<int> indices;
            syncmer::get_all_syncmers(walks[i], opt->s, curr_seq.kmers, curr_seq.flipped, curr_seq.unique_seq_len, indices);
            for (size_t j{0}; j < curr_seq.kmers.size(); j++){
                kmer_index[curr_seq.kmers[j]].emplace_back((int)i, indices[j], curr_seq.flipped[j], (int)j);
            }
        }
        walk_sketches.emplace_back(curr_seq);
    }

    std::vector<std::vector<int>> front_diff, back_diff;
    for (size_t i{0}; i < walks.size(); i++){
        front_diff.emplace_back(walks[i].size() + 1, 0);
        back_diff.emplace_back(walks[i].size() + 1, 0);
    }
    for (size_t i{0}; i < reads.size(); i++){
        std::string hpc_seq;
        std::vector<int> hp_lengths;
        hpc::get_compressed_sequence(reads[i], hpc_seq, hp_lengths);
        std::vector<uint64_t> kmer_list;
        std::vector<uint8_t> flipped_list;
        std::vector<int> unique_seq_len;
        std::vector<int> indices;
        syncmer::get_all_syncmers(hpc_seq, opt->s, kmer_list, flipped_list, unique_seq_len, indices);
        std::vector<ChainHit> hits;
        hpc::collect_all_chains(kmer_list, kmer_index, walk_sketches, indices, flipped_list, hits);
        for (const ChainHit& h : hits){
            long L = (long)walks[(size_t)h.utg_num].size();
            long a = std::max<long>(0, h.start);
            long b = std::min<long>(L, (long)h.start + h.len);
            if (b <= a) continue;
            long fa = (a == 0) ? a : std::min(b, a + W);
            front_diff[(size_t)h.utg_num][(size_t)fa]++;
            front_diff[(size_t)h.utg_num][(size_t)b]--;
            long bb = (b == L) ? b : std::max(a, b - W);
            back_diff[(size_t)h.utg_num][(size_t)a]++;
            back_diff[(size_t)h.utg_num][(size_t)bb]--;
        }
    }

    for (size_t i{0}; i < walks.size(); i++){
        long L = (long)walks[i].size();
        if (L < 2 * W + 2) continue;
        std::vector<int> front(walks[i].size(), 0), back(walks[i].size(), 0);
        int cf {0}, cb {0};
        for (long p{0}; p < L; p++){
            cf += front_diff[i][(size_t)p]; front[(size_t)p] = cf;
            cb += back_diff[i][(size_t)p]; back[(size_t)p] = cb;
        }
        long mid = L / 2;
        if (getenv("DK_END_TRIM_DUMP") != NULL){
            for (long p{0}; p < L; p += 25) std::cout << "[trim_dump] walk " << i << " pos " << p << " front " << front[(size_t)p] << " back " << back[(size_t)p] << std::endl;
        }
        long cut_front {0};
        for (long p{mid}; p >= 0; p--){
            if (front[(size_t)p] < MIN){
                cut_front = std::max<long>(0, p + 1 - W);
                break;
            }
        }
        long cut_back {L};
        for (long p{mid}; p < L; p++){
            if (back[(size_t)p] < MIN){
                cut_back = std::min<long>(L, p + W);
                break;
            }
        }
        if (cut_front > MAXLEN){
            if (opt->verbose) std::cout << "[trim_walk_ends] walk " << i << ": unsupported from " << cut_front << " bases inside the start, not trimming" << std::endl;
            cut_front = 0;
        }
        if (L - cut_back > MAXLEN){
            if (opt->verbose) std::cout << "[trim_walk_ends] walk " << i << ": unsupported from " << (L - cut_back) << " bases inside the end, not trimming" << std::endl;
            cut_back = L;
        }
        if (cut_front > 0 || cut_back < L){
            if (opt->verbose) std::cout << "[trim_walk_ends] walk " << i << " (" << L << " bp): trimming " << cut_front << " from the start and " << (L - cut_back) << " from the end" << std::endl;
            walks[i] = walks[i].substr((size_t)cut_front, (size_t)(cut_back - cut_front));
        }
    }
}

static dk_graph_t* dk_assemble_core(dk_opt_t* opt, char** reads_arr, char** read_names, int n_reads, dk_walks_t** walks_out){
    if (n_reads <= 0 || reads_arr == NULL) return NULL;
    /************************************************************
     * Initialization (convert to C++ structures, sketch reads) *
     ************************************************************/
    static const bool dinuc_initialized = (hash::initialize_dinucleotide_hashes(), true);
    (void)dinuc_initialized;
    std::vector<std::string> reads;
    std::vector<std::string> names;
    for (int i{0}; i < n_reads; i++){
        if (reads_arr[i] == NULL) return NULL;
        reads.emplace_back(reads_arr[i]);
        if (read_names != NULL && read_names[i] != NULL){
            names.emplace_back(read_names[i]);
        }else{
            names.emplace_back("read_" + std::to_string(i));
        }
    }
    std::vector<SketchedSeq*> input_seqs;
    syncmer::sketch_sequences(opt->s, reads, names, input_seqs);
    std::vector<SketchedSeq*> prev_seqs;
    std::vector<SketchedSeq*> backing_prev_graph;
    // to be used for phasing
    std::vector<SketchedSeq*> corrected_seqs;

    std::vector<size_t> read_kmer_counts;
    for (auto* seq : input_seqs){
        if (seq->is_read) read_kmer_counts.push_back(seq->kmers.size());
    }
    std::sort(read_kmer_counts.begin(), read_kmer_counts.end());
    int max_k_prime = (int)(read_kmer_counts[read_kmer_counts.size() / 2] / 2);

    /*************************************************************************
     * Loop until unitigs are all disconnected or max iterations are reached *
     *************************************************************************/
    Graph* prev_graph {NULL};
    int init_k_prime = opt->init_k_prime;
    int k_prime = init_k_prime;
    Graph* curr_graph = new Graph(k_prime);
    while(k_prime <= init_k_prime + opt->max_iterations && k_prime <= max_k_prime && (prev_graph == NULL || prev_graph->any_connected_unitigs())){
        if(opt->verbose){
            std::cout << "\n*****\nbeginning construction:" << std::endl;
        }
        IdxSet base_remove;
        int target_k_prime = prepare_iteration(k_prime, max_k_prime, input_seqs, prev_graph, base_remove);
        int k_lo = (prev_graph == NULL) ? target_k_prime : prev_graph->get_k_prime();
        bool keep_prev {false};
        IdxSet remove_set;
        while (true){
            curr_graph->set_k_prime(target_k_prime);
            build_graph(input_seqs, base_remove, prev_graph, curr_graph, remove_set);
            k_prime = curr_graph->get_k_prime();
            if (opt->verbose){
                std::cout << "k' = " << k_prime << std::endl;
                std::cout << input_seqs.size() - remove_set.size() << " sequences" << std::endl;
                std::cout << curr_graph->num_nodes() << " initial nodes" << std::endl;
                std::cout << curr_graph->num_unitigs() << " initial unitigs" << std::endl;
            }
            bool is_checkpoint = (opt->out_dir != NULL) && check_for_checkpoint(k_prime);
            if ((opt->out_dir != NULL) && (k_prime == init_k_prime || k_prime == init_k_prime + opt->max_iterations || is_checkpoint)){
                std::string out_unitigs {std::string(opt->out_dir) + "/initial_k" + std::to_string(k_prime) + ".gfa"};
                //std::string out_nodes {std::string(opt->out_dir) + "/initial_k" + std::to_string(k_prime) + "_nodes.gfa"};
                curr_graph->write_hpc_unitigs_to_gfa(out_unitigs, true);
                //curr_graph->write_nodes_to_gfa(out_nodes, true);
                std::string out_reads {std::string(opt->out_dir) + "/k" + std::to_string(k_prime) + "_reads.fa"};
                std::vector<SketchedSeq*> kept_seqs;
                for (size_t i{0}; i < input_seqs.size(); i++){
                    if (!remove_set.count(i)) kept_seqs.push_back(input_seqs[i]);
                }
                io::write_sequences(kept_seqs, out_reads);
            }
            bool check_dinuc {k_prime - init_k_prime > 150};
            curr_graph->trim_low_freq_unitigs(opt->compact_thresh, check_dinuc, (bool)opt->check_strand, opt->verbose);
            curr_graph->remove_floaters(opt->compact_thresh);
            if ((opt->out_dir != NULL) && (k_prime == init_k_prime || k_prime == init_k_prime + opt->max_iterations || is_checkpoint)){
                std::string out_hpc {std::string(opt->out_dir) + "/compacted_" + std::to_string(k_prime) + ".hpc.gfa"};
                curr_graph->write_hpc_unitigs_to_gfa(out_hpc, true);
            }
            if (opt->verbose){
                std::cout << curr_graph->num_nodes() << " nodes after compaction" << std::endl;
                std::cout << curr_graph->num_unitigs() << " unitigs after compaction" << std::endl;
            }
            if (prev_graph != NULL && graph_collapsed(curr_graph, prev_graph)){
                // k' step overshot what the reads support
                int k_hi = k_prime;
                if (k_hi - k_lo <= 1){
                    if (opt->verbose) std::cout << "graph collapsed at k' = " << k_hi << ", keeping k' = " << k_lo << std::endl;
                    keep_prev = true;
                    break;
                }
                target_k_prime = k_lo + (k_hi - k_lo) / 2;
                if (opt->verbose) std::cout << "graph collapsed at k' = " << k_hi << ", retrying at k' = " << target_k_prime << std::endl;
                delete curr_graph;
                curr_graph = new Graph(target_k_prime);
                continue;
            }
            break;
        }
        if (keep_prev) break;
        if (k_prime == 3){
            // deep copy the sequences so they aren't lost later
            // TODO: maybe just perform correction again but with the final unitigs?
            for (size_t i{0}; i < input_seqs.size(); i++){
                if (!remove_set.count(i) && input_seqs[i]->is_read) corrected_seqs.push_back(new SketchedSeq(*input_seqs[i]));
            }
        }
        curr_graph->trim_tips((float)opt->compact_thresh);
        retire_sequences(input_seqs, remove_set, prev_seqs);
        curr_graph->add_unitigs_to_input(input_seqs, opt->s);
        if (prev_graph != NULL){
            delete prev_graph;
            for (SketchedSeq* retired : backing_prev_graph){
                delete retired;
            }
        }
        backing_prev_graph = std::move(prev_seqs);
        k_prime = curr_graph->get_k_prime() + 1;
        prev_graph = curr_graph;
        curr_graph = new Graph(k_prime);
    }
    if (prev_graph == NULL){
        delete curr_graph;
        return NULL;
    }
    prev_graph->trim_tips((float)opt->compact_thresh);
    if(prev_graph->get_largest_component() == -1){
        prev_graph->identify_connected_components();
    }
    
    std::vector<std::string> phased_walks;

    std::unordered_map<std::string, int> name_to_assignment;
    std::vector<uint64_t> neutral_syncmers;
    phase::set_sketch_s(opt->s);
    if(prev_graph->num_unitigs() > 2){
        phase::phase_component(prev_graph->get_largest_component(), phased_walks, prev_graph, corrected_seqs, opt->verbose);
        if (prev_graph->num_components() > 1){
            for (size_t i{1}; i <= prev_graph->num_components(); i++){
                if (i != prev_graph->get_largest_component()){
                    phase::phase_component(i, phased_walks, prev_graph, corrected_seqs, opt->verbose);
                }
            }
        }
        for (SketchedSeq* seq : corrected_seqs){
            if (seq->assignment == HAP1 || seq->assignment == HAP2) name_to_assignment[seq->orig_name] = seq->assignment;
        }
    }else if (prev_graph->num_unitigs() == 2){
        std::vector<Unitig*> pair;
        const std::vector<Unitig*> utgs = prev_graph->view_unitigs();
        for (size_t i{0}; i < utgs.size(); i++) if (utgs[i] != NULL) pair.push_back(utgs[i]);
        if (pair.size() == 2){
            pair[0]->set_assignment(HAP1);
            pair[1]->set_assignment(HAP2);
            syncmers_without_counterpart(pair[0], pair[1], opt->s, neutral_syncmers);
        }
    }
    {
        std::vector<SketchedSeq*> all_reads;
        syncmer::sketch_sequences(opt->s, reads, names, all_reads);
        for (SketchedSeq* r : all_reads){
            auto it = name_to_assignment.find(r->orig_name);
            if (it != name_to_assignment.end()) r->assignment = it->second;
        }
        if (getenv("DK_NO_SYNCMER_VOTE") == NULL) phase::assign_reads_by_syncmers(prev_graph, all_reads, neutral_syncmers, opt->verbose);
        for (SketchedSeq* r : all_reads){
            if (r->assignment == HAP1 || r->assignment == HAP2) name_to_assignment[r->orig_name] = r->assignment;
            delete r;
        }
    }
    dk_graph_t* ret_graph = prev_graph->convert_to_gfa_representation();
    char** expanded_sequences;
    std::vector<std::string> final_unitigs;
    if (prev_graph->num_unitigs() <= 2){
        for (int i{0}; i < ret_graph->n_unitigs; i++){
            phased_walks.emplace_back(ret_graph->unitigs[i]);
        }
    }
    trim_walk_ends(opt, reads, phased_walks);
    {
        std::vector<char*> c_args;
        for (auto& s : phased_walks) {
            c_args.push_back(s.data());
        }
        c_args.push_back(nullptr);
        expanded_sequences = dk_expand_sequences(opt, reads_arr, n_reads, c_args.data(), (int)phased_walks.size());
    }
    for (size_t i{0}; i < phased_walks.size(); i++){
        final_unitigs.emplace_back(expanded_sequences[i]);
    }

    //if (opt->verbose) prev_graph->print_unitigs(true, false);
    if (opt->out_dir != NULL){
        //std::string out_reads {std::string(opt->out_dir) + "/final_reads.fa"};
        //io::write_sequences(input_seqs, out_reads);

        std::string out_hpc {std::string(opt->out_dir) + "/final_hpc.gfa"};
        prev_graph->write_hpc_unitigs_to_gfa(out_hpc, true);

        std::string out_broken {std::string(opt->out_dir) + "/broken_hpc.gfa"};
        curr_graph->write_hpc_unitigs_to_gfa(out_broken, true);

        std::string out_hpc_walks {std::string(opt->out_dir) + "/final_walks.hpc.fa"};
        io::write_walks(phased_walks, out_hpc_walks);

        std::string out_walks {std::string(opt->out_dir) + "/final_walks.fa"};
        io::write_walks(final_unitigs, out_walks);
    }

    if (opt->return_hpc_utg == 0){
        char** expanded_sequences = dk_expand_sequences(opt, reads_arr, n_reads, ret_graph->unitigs, ret_graph->n_unitigs);
        //free HPC versions of the unitigs
        for (int i{0}; i < ret_graph->n_unitigs; i++){
            delete[] ret_graph->unitigs[i];
        }
        ret_graph->unitigs = expanded_sequences;
    }

    if (walks_out != NULL){
        dk_walks_t* walks = new dk_walks_t();
        walks->n_walks = (int)final_unitigs.size();
        walks->walks = new char*[final_unitigs.size()];
        for (size_t i{0}; i < final_unitigs.size(); i++){
            walks->walks[i] = new char[final_unitigs[i].size() + 1];
            std::strcpy(walks->walks[i], final_unitigs[i].c_str());
        }
        walks->read_assignments = new int[(size_t)n_reads]();
        std::unordered_map<std::string, int> name_to_idx;
        for (int i{0}; i < n_reads; i++){
            name_to_idx[names[(size_t)i]] = i;
        }
        for (const auto& kv : name_to_assignment){
            auto it = name_to_idx.find(kv.first);
            if (it != name_to_idx.end() && kv.second > 0){
                walks->read_assignments[it->second] = kv.second;
            }
        }
        *walks_out = walks;
    }
    return ret_graph;
}

dk_graph_t* dk_assemble_reads(dk_opt_t* opt, char** reads_arr, char** read_names, int n_reads){
    return dk_assemble_core(opt, reads_arr, read_names, n_reads, NULL);
}

dk_walks_t* dk_assemble_walks(dk_opt_t* opt, char** reads_arr, char** read_names, int n_reads){
    dk_walks_t* walks = NULL;
    dk_graph_t* ret_graph = dk_assemble_core(opt, reads_arr, read_names, n_reads, &walks);
    if (ret_graph != NULL) dk_graph_free(ret_graph);
    return walks;
}

void dk_walks_free(dk_walks_t* to_free){
    for (int i{0}; i < to_free->n_walks; i++){
        delete[] to_free->walks[i];
    }
    delete[] to_free->walks;
    delete[] to_free->read_assignments;
    delete to_free;
}

char** dk_expand_sequences(dk_opt_t* opt, char** reads_arr, int n_reads, char** hpc_seqs_arr, int n_unitigs){
    if (n_unitigs <= 0){
        return new char*[1]();
    }
    /*****************************
     * Convert to C++ structures *
     *****************************/
    std::vector<std::string> hpc_seqs(hpc_seqs_arr, hpc_seqs_arr + n_unitigs);
    std::vector<std::string> reads(reads_arr, reads_arr + n_reads);

    /***************************************
     * Initialize HPC block data structure *
     ***************************************/
    std::vector<HPCCounts> hpc_tallies;
    for (size_t i{0}; i < hpc_seqs.size(); i++){
        hpc_tallies.emplace_back(hpc_seqs[i].size());
    }


    /**************************
     * Index unitig sequences *
     **************************/
    // add the selected syncmers within the unitig sequences to the index
    ankerl::unordered_dense::map<uint64_t, std::vector<SeqIdx>> kmer_index;
    // also keep the syncmer informtion for exact matching later
    std::vector<SketchedSeq> unitig_sketches;
    for (size_t i{0}; i < hpc_seqs.size(); i++){
        if (hpc_seqs[i].size() > 0){
            std::string curr_utg_name = "utg" + std::to_string(i);
            SketchedSeq curr_seq = SketchedSeq(curr_utg_name, false);
            std::vector<int> indices;
            syncmer::get_all_syncmers(hpc_seqs[i], opt->s, curr_seq.kmers, curr_seq.flipped, curr_seq.unique_seq_len, indices);
            for (size_t j{0}; j < curr_seq.kmers.size(); j++){
                kmer_index[curr_seq.kmers[j]].emplace_back(i, indices[j], curr_seq.flipped[j], j);
            }
            unitig_sketches.emplace_back(curr_seq);
        }else{
            unitig_sketches.emplace_back("", false);
        }
    }

    /***************
     * Align reads *
     ***************/
    for (size_t i{0}; i < reads.size(); i++){
        // first perform homopolymer compression on the reads
        std::string hpc_seq;
        std::vector<int> hp_lengths;
        hpc::get_compressed_sequence(reads[i], hpc_seq, hp_lengths);

        // select syncmers
        std::vector<uint64_t> kmer_list;
        std::vector<uint8_t> flipped_list;
        std::vector<int> unique_seq_len;
        std::vector<int> indices;
        syncmer::get_all_syncmers(hpc_seq, opt->s, kmer_list, flipped_list, unique_seq_len, indices);

        // align to the unitigs and record compression lengths
        hpc::align_single_read(kmer_list, kmer_index, unitig_sketches, hp_lengths, hpc_tallies, indices, flipped_list);
    }

    /********************
     * Expand sequences *
     ********************/
    std::vector<std::string> ret_seqs;
    for (size_t i{0}; i < hpc_seqs.size(); i++){
        std::string final_seq {""};
        int end {-1};
        for (size_t j{0}; j < hpc_tallies[i].len; j++){
            // some positions may not have coverage due to being at the end of a unitig
            if (hpc_tallies[i].coverages[j] < 3){
                if (j == 0){
                    if (opt->verbose) std::cout << "[dk_expand_sequences] trimming start of unitig " << i << std::endl;
                    continue;
                }else if (final_seq.size() == 0){
                    // still skipping over beginning of unitig
                    continue;
                }else{
                    // check if this is a temporary dip or actually the end
                    if (end == -1){
                        end = (int)hpc_tallies[i].len - 1;
                        while(end > 0 && hpc_tallies[i].coverages[end] < 3){
                            end--;
                        }
                    }
                    if (j <= end){
                        int num_to_add = (hpc_tallies[i].coverages[j] == 0) ? 1 : (int)std::round((float)hpc_tallies[i].sums[j] / (float)hpc_tallies[i].coverages[j]);
                        for (int k{0}; k < num_to_add; k++){
                            final_seq += hpc_seqs[i][j];
                        }
                    }else{
                        if (opt->verbose) std::cout << "[dk_hpc::expand_sequences] truncating unitig " << i << " at HPC location " << j << " (out of " << hpc_tallies[i].len << ")" << std::endl;
                        break;
                    }
                }
            }else{
                int num_to_add = (int)std::round((float)hpc_tallies[i].sums[j] / (float)hpc_tallies[i].coverages[j]);
                for (int k{0}; k < num_to_add; k++){
                    final_seq += hpc_seqs[i][j];
                }
            }
        }
        ret_seqs.push_back(std::move(final_seq));
    }
    // convert to C-style array
    char** ret_arr = new char*[ret_seqs.size()];
    for (size_t i{0}; i < ret_seqs.size(); i++){
        ret_arr[i] = new char[ret_seqs[i].size() + 1];
        std::strcpy(ret_arr[i], ret_seqs[i].c_str());
    }
    return ret_arr;
}

void dk_write_gfa(const char* file_path, dk_graph_t* graph){
    std::ofstream out_file(static_cast<const char*>(file_path));
    for (int i{0}; i < graph->n_unitigs; i++){
        out_file << "S\tutg" << i << '\t' << graph->unitigs[i] << '\n';
    }
    for (int i{0}; i < graph->n_edges; i++){
        out_file << "L\tutg" << graph->edges[i].from_utg << '\t' << graph->edges[i].from_dir << "\tutg" << graph->edges[i].to_utg << '\t' << graph->edges[i].to_dir << "\t0M\n";
    }
    out_file.close();
}

void dk_graph_free(dk_graph_t* to_free){
    for (int i{0}; i < to_free->n_unitigs; i++){
        delete[] to_free->unitigs[i];
    }
    delete[] to_free->unitigs;
    delete[] to_free->edges;
    delete[] to_free;
}

static int prepare_iteration(int k_prime, int max_k_prime, std::vector<SketchedSeq*>& input_sequences, Graph* prev_graph, IdxSet& remove_set){
    std::vector<SketchedSeq*> corrected_seqs;
    for (size_t seq_idx {0}; seq_idx < input_sequences.size(); seq_idx++){
        if (input_sequences[seq_idx]->kmers.size() < k_prime){
            remove_set.insert(seq_idx);
        }else if (prev_graph != NULL && input_sequences[seq_idx]->is_read){
            SketchedSeq* corrected_seq = NULL;
            SketchedSeq* split_seq = NULL;
            bool contained = prev_graph->check_sequence(input_sequences[seq_idx], (size_t)k_prime - 1, corrected_seq, split_seq);
            if (contained || corrected_seq != NULL){
                remove_set.insert(seq_idx);
                if (corrected_seq != NULL){
                    corrected_seqs.push_back(corrected_seq);
                }
            }
            if (split_seq != NULL){
                corrected_seqs.push_back(split_seq);
            }
        }
    }
    // if we're skipping iterations, update k' value for current graph
    if(prev_graph != NULL){
        k_prime += std::max(prev_graph->min_connect_size() - 2, 0);
        if (k_prime > max_k_prime) k_prime = max_k_prime;
    }
    for (size_t i{0}; i < corrected_seqs.size(); i++){
        input_sequences.emplace_back(corrected_seqs[i]);
    }
    return k_prime;
}

static void build_graph(std::vector<SketchedSeq*>& input_sequences, const IdxSet& base_remove, Graph* prev_graph, Graph* dbg, IdxSet& remove_set){
    remove_set = base_remove;
    if (prev_graph != NULL){
        dbg->reserve_nodes((size_t)prev_graph->num_nodes() + (size_t)prev_graph->num_nodes() / 8);
    }
    for (size_t seq_idx {0}; seq_idx < input_sequences.size(); seq_idx++){
        if((remove_set.find(seq_idx) == remove_set.end())  && (input_sequences[seq_idx]->kmers.size() < dbg->get_k_prime() || !dbg->add_sequence(*input_sequences[seq_idx], prev_graph))){
            remove_set.insert(seq_idx);
        }
    }
}

static void retire_sequences(std::vector<SketchedSeq*>& input_sequences, const IdxSet& remove_set, std::vector<SketchedSeq*>& prev_seqs){
    prev_seqs.clear();
    std::vector<SketchedSeq*> kept;
    kept.reserve(input_sequences.size() - remove_set.size());
    for (size_t i{0}; i < input_sequences.size(); i++){
        if (remove_set.count(i)){
            prev_seqs.push_back(input_sequences[i]);
        }else{
            kept.push_back(input_sequences[i]);
        }
    }
    input_sequences = std::move(kept);
}

static bool graph_collapsed(Graph* curr_graph, Graph* prev_graph){
    if (curr_graph->num_unitigs() == 0) return true;
    static const double collapse_frac = getenv("DK_COLLAPSE_FRAC") ? atof(getenv("DK_COLLAPSE_FRAC")) : COLLAPSE_FRAC;
    if (collapse_frac <= 0) return false;
    return (double)curr_graph->total_unitig_length() < collapse_frac * (double)prev_graph->total_unitig_length();
}

bool check_for_checkpoint(int curr_k_prime){
    const int num_checkpoints = (int)(sizeof(DEBUG_CHECKPOINTS) / sizeof(DEBUG_CHECKPOINTS[0]));
    int new_idx = DEBUG_IDX;
    while(new_idx < num_checkpoints && curr_k_prime >= DEBUG_CHECKPOINTS[new_idx]){
        new_idx++;
    }
    if (new_idx != DEBUG_IDX){
        DEBUG_IDX = new_idx;
        return true;
    }
    return false;
}
