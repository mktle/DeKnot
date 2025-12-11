#ifndef GRAPH_H
#define GRAPH_H

#include "deknot.h"
#include "node.h"
#include "span.h"
#include "syncmer.h"
#include "unitig.h"
#include "utils.h"
#include "xxhash.h"

#include <ankerl/unordered_dense.h>
#include <vector>

// corresponds to a part of an input sequence that is not found in the cleaned graph
// i.e., this part was probably dropped because of an error
// so we will fix it
struct DroppedBlock {
    Node* end1;
    Node* end2;
    size_t start_idx;
    size_t end_idx;
    bool across_utgs;

    DroppedBlock(Node* e1, Node* e2, size_t s_idx, size_t e_idx, bool a) : end1(e1), end2(e2), start_idx(s_idx), end_idx(e_idx), across_utgs(a) {}
};

class Graph {
    public:
        Graph(int set_k_prime);
        ~Graph();
        int get_k_prime();
        void set_k_prime(int set_val);
        int num_nodes();
        int num_unitigs();
        bool check_sequence(SketchedSeq* orig_seq, size_t check_k_prime, int& min_span, SketchedSeq*& corrected_seq);
        void correct_read_error_within_unitig(DroppedBlock& to_correct, SketchedSeq*& corrected_seq);
        void correct_read_error_across_unitigs(DroppedBlock& dropped_block, SketchedSeq*& corrected_seq);
        bool add_sequence(SketchedSeq& input_seq, Graph* prev_graph);
        bool add_connection(KmerView& first, char first_orient, KmerView& second, char second_orient);
        void handle_self_connection(Node* node);
        void check_endpoints(Node* first_node, Node* second_node);

        Unitig* split_unitig(Unitig* to_split, Node* breakpoint);
        void merge_unitigs(Node* first_node, Node* second_node);
        int next_node_id();
        int next_unitig_id();

        void increment_k_min_mer_count(int increment, KmerView& k_syncmer, Span<uint8_t>& kmer_flipped, Span<int>& unique_seq_idx);
        void increment_unitig(int increment, Unitig* to_increment);
        int k_plus_one_mer_count(Span<uint64_t>& k_syncmer);
        void destroy_unitig(Unitig* to_remove);
        void destroy_node(Node* to_remove);
        void add_node_to_list(Node* new_node);
        void add_unitig_to_list(Unitig* new_utg);
        void create_unitig_for_node(Node* node);

        void trim_tips();
        void remove_floaters(int threshold);
        bool trim_low_freq_unitigs(int threshold, bool check_dinuc);
        bool trim_low_freq_end(Node* node_end, float thresh, std::vector<Unitig*>& sorted_list, bool check_dinuc);
        bool drop_unitig(std::vector<Unitig*>& sorted_list, Node* curr_node, Unitig* neigh_utg);
        bool dinuc_topology(Node* curr_end, Unitig* dropped_utg);
        bool disconnect_unitig_end(Node* remove_node, Unitig* keep_utg, std::vector<Unitig*>& sorted_list);
        void add_node_to_lookup_table(Node* to_add);
        Node* get_node_from_lookup_table(KmerView& k_syncmer);

        bool any_connected_unitigs();

        void write_nodes_to_gfa(std::string& file_name, bool verbose);
        void write_hpc_unitigs_to_gfa(std::string& file_name, bool verbose);
        void free_unitig_if_empty(Unitig* to_check);
        void print_nodes(bool verbose);
        void print_unitigs(bool verbose, bool print_seq);
        void get_all_hpc_seqs(std::vector<std::string>& hpc_seqs);
        void add_unitigs_to_input(std::vector<SketchedSeq*>& input_seqs, int s);

        void find_paths(std::vector<Unitig*> path1, std::vector<Unitig*> path2);
        void path_to_sequence(std::vector<Unitig*> path, std::string& ret_string);

        dk_graph_t* convert_to_gfa_representation();

    private: std::vector<Node*> node_list;
        int k_prime;
        std::vector<int> node_free_list;
        std::vector<Unitig*> unitig_list;
        std::vector<int> unitig_free_list;

        ankerl::unordered_dense::map<KmerView, Node*, utils::KmerViewHash, utils::KmerViewEqual> node_table;
        // inspect and trim branches that have these kmers at the branch site
        // e.g., dinucleotide repeats
        ankerl::unordered_dense::set<uint64_t> drop_sites;
};

#endif
