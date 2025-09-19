#ifndef NODE_H
#define NODE_H

#include "hpc.h"
#include "span.h"
#include "syncmer.h"
#include "unitig.h"

#include <functional>
#include <numeric>
#include <vector>

class Node;
class Unitig;

struct Edge{
    Node* neighbor;
    char self_orientation;
    char neighbor_orientation;

    Edge(Node* neigh, char self_o, char n_o): neighbor(neigh), self_orientation(self_o), neighbor_orientation(n_o) {}
};

class Node{
    public:
        Node(KmerView& set_ksyncmer, Span<uint8_t>& set_flipped, Span<int>& set_unique_len, int set_node, Unitig* set_unitig);
        void set_count(int n);
        void add_neighbor(Node* new_neighbor, char curr_orient, char neigh_orient);
        void remove_neighbor(Node* remove_neighbor);
        int num_edges_in_dir(char check_dir);
        void assign_unitig(Unitig* new_ptr);
        bool is_assigned_to_unitig();
        bool is_neighbors_with(Node* to_check);
        uint64_t unique_kmer_from_neigh(Node* neigh, bool& flipped, int& unique_seq_len);
        int num_neighbors();
        int num_kmers();
        int get_node_id();
        int get_count();
        uint64_t ith_kmer(int i);
        bool ith_kmer_flipped(int i);
        int ith_unique_seq_len(int i);
        Unitig* get_unitig();
        KmerView get_kmer_view();
        uint64_t get_front_hash();
        uint64_t get_back_hash();
        bool neighbors_diff_utg();
        int longest_neighbor_utg_length();
        const std::vector<Edge>& view_neighbors();
        bool is_tip();
        void get_hpc_sequence(std::string& ret_str, bool reverse);
        bool end_is_dinucleotide();
        void print_kmers();
        void print(bool verbose);

    private:
        KmerView ksyncmer;
        // must use uint8_t instead of bool because the bool vectors are handled differently in memory
        Span<uint8_t> flip_kmer;
        Span<int> unique_seq_len;

        std::vector<Edge> neighbors;

        // frequency of this node
        int count;

        int node_idx;
        Unitig* utg_ptr;

        // track how many edges are on each side of the node
        // used for unitig merging/splitting code
        int num_forward_edges;
        int num_reverse_edges;
};

#endif
