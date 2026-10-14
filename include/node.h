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

    Edge() : neighbor(NULL), self_orientation(0), neighbor_orientation(0) {}
    Edge(Node* neigh, char self_o, char n_o): neighbor(neigh), self_orientation(self_o), neighbor_orientation(n_o) {}
};

class EdgeList{
    public:
        static constexpr unsigned INLINE_CAP = 3;
        using iterator = Edge*;
        using const_iterator = const Edge*;

        EdgeList() : n(0), cap(INLINE_CAP), heap(NULL) {}
        ~EdgeList() { delete[] this->heap; }
        EdgeList(const EdgeList&) = delete;
        EdgeList& operator=(const EdgeList&) = delete;

        size_t size() const { return this->n; }
        bool empty() const { return this->n == 0; }
        Edge& operator[](size_t i) { return this->data()[i]; }
        const Edge& operator[](size_t i) const { return this->data()[i]; }
        iterator begin() { return this->data(); }
        iterator end() { return this->data() + this->n; }
        const_iterator begin() const { return this->data(); }
        const_iterator end() const { return this->data() + this->n; }
        void clear() { this->n = 0; }
        void emplace_back(Node* neigh, char self_o, char n_o){
            if (this->n == this->cap) this->grow();
            this->data()[this->n] = Edge(neigh, self_o, n_o);
            this->n++;
        }
        iterator erase(iterator pos){
            for (Edge* p = pos; p + 1 < this->end(); p++) *p = *(p + 1);
            this->n--;
            return pos;
        }
        std::vector<Edge> to_vector() const { return std::vector<Edge>(this->begin(), this->end()); }

    private:
        Edge inl[INLINE_CAP];
        unsigned n;
        unsigned cap;
        Edge* heap;

        Edge* data() { return this->heap != NULL ? this->heap : this->inl; }
        const Edge* data() const { return this->heap != NULL ? this->heap : this->inl; }
        void grow(){
            unsigned new_cap = this->cap * 2;
            Edge* grown = new Edge[new_cap];
            for (unsigned i{0}; i < this->n; i++) grown[i] = this->data()[i];
            delete[] this->heap;
            this->heap = grown;
            this->cap = new_cap;
        }
};

class Node{
    public:
        Node(KmerView& set_ksyncmer, Span<uint8_t>& set_flipped, Span<int>& set_unique_len, int set_node, Unitig* set_unitig);
        void reset(KmerView& set_ksyncmer, Span<uint8_t>& set_flipped, Span<int>& set_unique_len, int set_node, Unitig* set_unitig);
        void set_count(int n);
        void add_neighbor(Node* new_neighbor, char curr_orient, char neigh_orient);
        void remove_neighbor(Node* remove_neighbor);
        int num_edges_in_dir(char check_dir);
        void assign_unitig(Unitig* new_ptr);
        bool is_assigned_to_unitig() { return this->utg_ptr != NULL; }
        bool is_neighbors_with(Node* to_check);
        uint64_t unique_kmer_from_neigh(Node* neigh, bool& flipped, int& unique_seq_len);
        int num_neighbors() { return (int)this->neighbors.size(); }
        int num_kmers() { return (int)this->ksyncmer.kmers.size(); }
        int get_node_id() { return this->node_idx; }
        int get_count() { return this->count; }
        void add_strand_obs(bool flipped, int inc);
        int get_strand_count(bool flipped);
        int minority_strand_count();
        uint64_t ith_kmer(int i);
        bool ith_kmer_flipped(int i);
        int ith_unique_seq_len(int i);
        Unitig* get_unitig() { return this->utg_ptr; }
        KmerView get_kmer_view() { return this->ksyncmer; }
        uint64_t get_front_hash();
        uint64_t get_back_hash();
        bool neighbors_diff_utg();
        Node* get_diff_unitig_neighbor();
        Node* get_parallel_neigh(Unitig* in_utg);
        int longest_neighbor_utg_length();
        const EdgeList& view_neighbors() { return this->neighbors; }
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

        EdgeList neighbors;

        // frequency of this node
        int count;
        // strand counts for ONT error correction
        int strand_counts[2] = {0, 0};

        int node_idx;
        Unitig* utg_ptr;

        // track how many edges are on each side of the node
        // used for unitig merging/splitting code
        int num_forward_edges;
        int num_reverse_edges;
};

#endif
