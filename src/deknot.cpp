#include "deknot.h"
#include "graph.h"
#include "hash.h"
#include "io.h"
#include "phase.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <mutex>
#include <string>
#include <unordered_map>

const int DEBUG_CHECKPOINTS[] = {2, 3, 50, 100, 110, 120, 130, 140, 150, 200, 210, 220, 230, 238, 240, 250, 300, 350, 400, 450, 466, 468, 469, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950};
//const int DEBUG_CHECKPOINTS[] = {2};
int DEBUG_IDX = 0;

bool construct_graph(int k_prime, int max_k_prime, std::vector<SketchedSeq*>& prev_seqs, std::vector<SketchedSeq*>& input_sequences, Graph* prev_graph, Graph* dbg, bool verbose);

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

static dk_graph_t* dk_assemble_core(dk_opt_t* opt, char** reads_arr, char** read_names, int n_reads, dk_walks_t** walks_out){
    if (n_reads <= 0 || reads_arr == NULL) return NULL;
    /************************************************************
     * Initialization (convert to C++ structures, sketch reads) *
     ************************************************************/
    static std::once_flag dinuc_init_flag;
    std::call_once(dinuc_init_flag, hash::initialize_dinucleotide_hashes);
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
        if (!construct_graph(k_prime, max_k_prime, prev_seqs, input_seqs, prev_graph, curr_graph, opt->verbose)){
            return NULL;
        }
        k_prime = curr_graph->get_k_prime();
        if (k_prime == 3){
            // deep copy the sequences so they aren't lost later
            // TODO: maybe just perform correction again but with the final unitigs?
            for (auto* seq : input_seqs){
                if (seq->is_read) corrected_seqs.push_back(new SketchedSeq(*seq));

            }
        }
        if (opt->verbose){
            std::cout << "k' = " << k_prime << std::endl;
            std::cout << input_seqs.size() << " sequences" << std::endl;
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
            io::write_sequences(input_seqs, out_reads);
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
        curr_graph->trim_tips();
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
    prev_graph->trim_tips();
    if(prev_graph->get_largest_component() == -1){
        prev_graph->identify_connected_components();
    }
    
    std::vector<std::string> phased_walks;

    if(prev_graph->num_unitigs() > 2){
        phase::phase_component(prev_graph->get_largest_component(), phased_walks, prev_graph, corrected_seqs, opt->verbose);
        if (prev_graph->num_components() > 1){
            for (size_t i{1}; i <= prev_graph->num_components(); i++){
                if (i != prev_graph->get_largest_component()){
                    phase::phase_component(i, phased_walks, prev_graph, corrected_seqs, opt->verbose);
                }
            }
        }
    }
    dk_graph_t* ret_graph = prev_graph->convert_to_gfa_representation();
    char** expanded_sequences;
    std::vector<std::string> final_unitigs;
    if (prev_graph->num_unitigs() <= 2){
        expanded_sequences = dk_expand_sequences(opt, reads_arr, n_reads, ret_graph->unitigs, ret_graph->n_unitigs);
        for (int i{0}; i < ret_graph->n_unitigs; i++){
            phased_walks.emplace_back(ret_graph->unitigs[i]);
        }
    }else{
        std::vector<char*> c_args;
        for (auto& s : phased_walks) {
            c_args.push_back(s.data());
        }
        c_args.push_back(nullptr);
        expanded_sequences = dk_expand_sequences(opt, reads_arr, n_reads, c_args.data(), phased_walks.size());
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
        for (SketchedSeq* seq : corrected_seqs){
            auto it = name_to_idx.find(seq->name);
            if (it != name_to_idx.end() && seq->assignment > 0){
                walks->read_assignments[it->second] = seq->assignment;
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

bool construct_graph(int k_prime, int max_k_prime, std::vector<SketchedSeq*>& prev_seqs, std::vector<SketchedSeq*>& input_sequences, Graph* prev_graph, Graph* dbg, bool verbose){
    ankerl::unordered_dense::set<size_t> remove_set;
    std::vector<SketchedSeq*> corrected_seqs;
    // check the sequences to see if they should be removed (and also calculate jump length)
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
    dbg->set_k_prime(k_prime);

    for (size_t i{0}; i < corrected_seqs.size(); i++){
        input_sequences.emplace_back(corrected_seqs[i]);
    }
    // actually add the sequences now
    for (size_t seq_idx {0}; seq_idx < input_sequences.size(); seq_idx++){
        if((remove_set.find(seq_idx) == remove_set.end())  && (input_sequences[seq_idx]->kmers.size() < dbg->get_k_prime() || !dbg->add_sequence(*input_sequences[seq_idx], prev_graph))){
            remove_set.insert(seq_idx);
        }
    }
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
    return true;
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
