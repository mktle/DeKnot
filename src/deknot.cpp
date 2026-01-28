#include "deknot.h"
#include "graph.h"
#include "hash.h"
#include "io.h"
#include "phase.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

const int DEBUG_CHECKPOINTS[] = {3, 24};

bool construct_graph(int k_prime, std::vector<SketchedSeq*>& prev_seqs, std::vector<SketchedSeq*>& input_sequences, Graph* prev_graph, Graph* dbg, bool verbose);

void dk_opt_init(dk_opt_t* opt){
    memset(opt, 0, sizeof(dk_opt_t));
    opt->s = 9;
    opt->init_k_prime = 350;
    opt->max_iterations = 1000;//mkl
    opt->compact_thresh = 4;
    opt->return_hpc_utg = 0;
    opt->out_dir = NULL;
    opt->verbose = 0;
}

int dk_set_opt(const char* preset, dk_opt_t* opt){
    //TODO: establish presets
    assert(preset != NULL);
    assert(opt != NULL);
    return 0;
}

dk_graph_t* dk_assemble_reads(dk_opt_t* opt, char** reads_arr, char** read_names, int n_reads){
    /************************************************************
     * Initialization (convert to C++ structures, sketch reads) *
     ************************************************************/
    hash::initialize_dinucleotide_hashes();
    std::vector<std::string> reads(reads_arr, reads_arr + n_reads);
    std::vector<std::string> names;
    if (read_names != NULL){
        for (int i{0}; i < n_reads; i++){
            names.emplace_back(read_names[i]);
        }
    }
    std::vector<SketchedSeq*> input_seqs;
    syncmer::sketch_sequences(opt->s, reads, names, input_seqs);
    std::vector<SketchedSeq*> prev_seqs;

    /*************************************************************************
     * Loop until unitigs are all disconnected or max iterations are reached *
     *************************************************************************/
    Graph* prev_graph {NULL};
    int init_k_prime = opt->init_k_prime;
    int k_prime = init_k_prime;
    Graph* curr_graph = new Graph(k_prime);
    while(k_prime <= init_k_prime + opt->max_iterations && (prev_graph == NULL || prev_graph->any_connected_unitigs())){
        if(opt->verbose){
            std::cout << "\n*****\nbeginning construction:" << std::endl;
        }
        if (!construct_graph(k_prime, prev_seqs, input_seqs, prev_graph, curr_graph, opt->verbose)){
            return NULL;
        }
        k_prime = curr_graph->get_k_prime();
        if (opt->verbose){
            std::cout << "k' = " << k_prime << std::endl;
            std::cout << input_seqs.size() << " sequences" << std::endl;
            std::cout << curr_graph->num_nodes() << " initial nodes" << std::endl;
            std::cout << curr_graph->num_unitigs() << " initial unitigs" << std::endl;
        }
        if ((opt->out_dir != NULL) && (k_prime == init_k_prime || k_prime == init_k_prime + opt->max_iterations || (k_prime - init_k_prime) % 100 == 0 || std::find(std::begin(DEBUG_CHECKPOINTS), std::end(DEBUG_CHECKPOINTS), k_prime) != std::end(DEBUG_CHECKPOINTS))){
            std::string out_unitigs {std::string(opt->out_dir) + "/initial_k" + std::to_string(k_prime) + ".gfa"};
            std::string out_nodes {std::string(opt->out_dir) + "/initial_k" + std::to_string(k_prime) + "_nodes.gfa"};
            curr_graph->write_hpc_unitigs_to_gfa(out_unitigs, true);
            curr_graph->write_nodes_to_gfa(out_nodes, true);
            std::string out_reads {std::string(opt->out_dir) + "/k" + std::to_string(k_prime) + "_reads.fa"};
            io::write_sequences(input_seqs, out_reads);
        }
        int times {0};
        bool check_dinuc {k_prime - init_k_prime > 20};
        bool graph_changed = curr_graph->trim_low_freq_unitigs(opt->compact_thresh, check_dinuc);
        while(graph_changed && times < 10){
            curr_graph->trim_tips();
            graph_changed = curr_graph->trim_low_freq_unitigs(opt->compact_thresh, check_dinuc);
            times++;
        }
        curr_graph->remove_floaters(opt->compact_thresh);
        if ((opt->out_dir != NULL) && (k_prime == init_k_prime || k_prime == init_k_prime + opt->max_iterations || (k_prime - init_k_prime) % 100 == 0 || std::find(std::begin(DEBUG_CHECKPOINTS), std::end(DEBUG_CHECKPOINTS), k_prime) != std::end(DEBUG_CHECKPOINTS))){
            std::string out_unitigs_compacted {std::string(opt->out_dir) + "/compacted_k" + std::to_string(k_prime) + ".gfa"};
            std::string out_hpc {std::string(opt->out_dir) + "/compacted_" + std::to_string(k_prime) + ".hpc.gfa"};
            curr_graph->write_hpc_unitigs_to_gfa(out_hpc, true);
            dk_graph_t* curr_gfa = curr_graph->convert_to_gfa_representation();
            char** curr_expanded = dk_expand_sequences(opt, reads_arr, n_reads, curr_gfa->unitigs, curr_gfa->n_unitigs);
            //free HPC versions of the unitigs
            for (int i{0}; i < curr_gfa->n_unitigs; i++){
                delete[] curr_gfa->unitigs[i];
            }
            curr_gfa->unitigs = curr_expanded;
            dk_write_gfa(out_unitigs_compacted.c_str(), curr_gfa);
            dk_graph_free(curr_gfa);
        }
        if (opt->verbose){
            std::cout << curr_graph->num_nodes() << " nodes after compaction" << std::endl;
            std::cout << curr_graph->num_unitigs() << " unitigs after compaction" << std::endl;
            //curr_graph->print_unitigs(true, false);
        }
        curr_graph->add_unitigs_to_input(input_seqs, opt->s);
        if (prev_graph != NULL){
            delete prev_graph;
        }
        k_prime = curr_graph->get_k_prime() + 1;
        prev_graph = curr_graph;
        curr_graph = new Graph(k_prime);
    }
    prev_graph->print_unitigs(true, false);//mkl
    
    // this will be filled by next function call
    std::vector<Unitig*> utg_list;
    phase::initialize_unitigs_for_tagging(prev_graph, input_seqs, utg_list);
    
    while(utg_list.size() > 0){
        phase::choose_initial_bubble(utg_list);
        int i{0};
        while(phase::reassign_partitions(prev_graph, input_seqs)){
            std::cout << "After phasing round " << i << std::endl;//mkl
            i++;
            prev_graph->print_unitigs(true, false);//mkl
        }
        phase::remove_tagged_unitigs(utg_list);
    }

    std::cout << "After:" << std::endl;//mkl
    prev_graph->print_unitigs(true, false);//mkl

    std::string out_reads {std::string(opt->out_dir) + "/final_reads.fa"};//mkl
    io::write_sequences(input_seqs, out_reads);//mkl

    if (opt->out_dir != NULL){
        std::string out_hpc {std::string(opt->out_dir) + "/final_hpc.gfa"};
        prev_graph->write_hpc_unitigs_to_gfa(out_hpc, true);
    }
    std::vector<std::string> phased_walks;//mkl
    phase::get_haplotype_sequences(prev_graph, phased_walks);//mkl
    std::string out_walks {std::string(opt->out_dir) + "/final_walks.fa"};//mkl
    io::write_walks(phased_walks, out_walks);//mkl

    dk_graph_t* ret_graph = prev_graph->convert_to_gfa_representation();
    if (opt->return_hpc_utg == 0){
        char** expanded_sequences = dk_expand_sequences(opt, reads_arr, n_reads, ret_graph->unitigs, ret_graph->n_unitigs);
        //free HPC versions of the unitigs
        for (int i{0}; i < ret_graph->n_unitigs; i++){
            delete[] ret_graph->unitigs[i];
        }
        ret_graph->unitigs = expanded_sequences;
    }
    return ret_graph;
}

char** dk_expand_sequences(dk_opt_t* opt, char** reads_arr, int n_reads, char** hpc_seqs_arr, int n_unitigs){
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
        for (size_t j{0}; j < hpc_tallies[i].len; j++){
            // some positions may not have coverage due to being at the end of a unitig
            if (hpc_tallies[i].coverages[j] < 3){
                if (j == 0){
                    std::cout << "[hpc::expand_sequences] trimming start of unitig " << i << std::endl;
                    continue;
                }else if (final_seq.size() == 0){
                    // still skipping over beginning of unitig
                    continue;
                }else{
                    std::cout << "[hpc::expand_sequences] truncating unitig " << i << " at HPC location " << j << " (out of " << hpc_tallies[i].len << ")" << std::endl;
                    break;
                }
            }
            int num_to_add = (int)std::round((float)hpc_tallies[i].sums[j] / (float)hpc_tallies[i].coverages[j]);
            for (int k{0}; k < num_to_add; k++){
                final_seq += hpc_seqs[i][j];
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

bool construct_graph(int k_prime, std::vector<SketchedSeq*>& prev_seqs, std::vector<SketchedSeq*>& input_sequences, Graph* prev_graph, Graph* dbg, bool verbose){
    std::vector<size_t> remove_list;
    std::vector<SketchedSeq*> corrected_seqs;
    // check the sequences to see if they should be removed (and also calculate jump length)
    int lowest_span {-1};
    for (size_t seq_idx {0}; seq_idx < input_sequences.size(); seq_idx++){
        int curr_span {-1};
        if (input_sequences[seq_idx]->kmers.size() < k_prime){
            remove_list.push_back(seq_idx);
        }else if (prev_graph != NULL && input_sequences[seq_idx]->is_read){
            SketchedSeq* corrected_seq = NULL;
            bool contained = prev_graph->check_sequence(input_sequences[seq_idx], (size_t)k_prime - 1, curr_span, corrected_seq);
            if (contained || corrected_seq != NULL){
                remove_list.push_back(seq_idx);
                if (corrected_seq != NULL){
                    corrected_seqs.push_back(corrected_seq);
                }
            }
        }
        if (curr_span != -1){
            // valid sequence
            if (lowest_span == -1 || curr_span < lowest_span) lowest_span = curr_span;
        }
    }

    // if we're skipping iterations, update k' value for current graph
    if (lowest_span != -1){
        k_prime += lowest_span - 1;
        dbg->set_k_prime(k_prime);
    }

    if (verbose){
        std::cout << "Lowest span: " << lowest_span << std::endl;
    }
    for (size_t i{0}; i < corrected_seqs.size(); i++){
        input_sequences.emplace_back(corrected_seqs[i]);
    }
    // actually add the sequences now
    for (size_t seq_idx {0}; seq_idx < input_sequences.size(); seq_idx++){
        //TODO: switch from linear search for index
        if((std::find(remove_list.begin(), remove_list.end(), seq_idx) == remove_list.end())  && (input_sequences[seq_idx]->kmers.size() < dbg->get_k_prime() || !dbg->add_sequence(*input_sequences[seq_idx], prev_graph))){
            remove_list.push_back(seq_idx);
        }
    }
    std::sort(remove_list.begin(), remove_list.end());

    // free the sequences that were removed from the previous graph
    for (size_t i{0}; i < prev_seqs.size(); i++){
        delete prev_seqs[i];
    }
    prev_seqs.clear();

    // remove discarded sequences from the next graph's input
    // cannot delete them just yet because the next graph may need to access their nodes
    for(int i{(int)remove_list.size() - 1}; i >= 0; i--){
        prev_seqs.push_back(input_sequences[remove_list[(size_t)i]]);
        input_sequences.erase(input_sequences.begin() + (long)remove_list[(size_t)i]);
    }
    return true;
}
