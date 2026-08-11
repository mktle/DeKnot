#ifndef UNITIG_H
#define UNITIG_H

#include "hash.h"
#include "node.h"

#include <vector>

class Node;

class Unitig{
    public:
        Unitig(std::vector<Node*>& set_nodes, int index);

        bool neighbors_unitig();
        bool node_is_endpoint(Node* to_check);
        int num_nodes();
        Node* get_front();
        Node* get_back();
        Node* opposite_end(Node* in_end);
        int get_id();

        const std::vector<Node*>& view_nodes();
        Node* node_at(size_t idx);

        void add_node_front(Node* to_add);
        void add_node_back(Node* to_add);
        bool contains_node(Node* check);

        void drain_nodes(std::vector<Node*>& recipient, Node* breakpoint);
        void clear_nodes();
        Node* pop_front();
        Node* pop_back();

        float get_median_freq();
        float get_base_freq();
        int min_minority_strand_count(int& support);
        bool self_connected();
        bool connects_bubbles();
        void update_median_freq();
        void set_base_freq(float set_freq);
        bool is_floating();
        Node* get_tip_inner_end(Node*& adj_node);
        char get_in_dir(Node* node_to_check);
        bool is_neighbors_with(Unitig* check_utg);
            
        int get_assignment();
        void set_assignment(int id);
        void set_hap_count(int id, int count);
        int get_hap_count(int id);
        float get_hap_frac(int id);
        void set_hap_updated(bool set);
        bool get_hap_updated();

        void set_parallel_utg();
        Unitig* get_parallel_utg();
        void add_read(SketchedSeq*);
        void assign_reads(bool lock);
        void tally_reads();
        void lock_haplotype();
        bool is_locked();
        void set_times_walked(int val);
        int get_times_walked();
        void set_walk_limit(int val);
        int get_walk_limit();
        const std::vector<SketchedSeq*>& view_aligned_seqs();

        void list_nodes();
        void get_hpc_sequence(std::string& ret_seq, bool reverse, bool append);
        void print_kmers();
        void print(bool verbose);

    private:
        std::vector<Node*> nodes;
        int index;
        float median_freq;
        float base_freq;
        // min minority strand count (for ONT errors)
        int min_minority_count;
        int min_minority_supp;
        /***************
         * For phasing *
         ***************/
        std::vector<SketchedSeq*> aligned_seqs;
        // haplotype assignment
        int assignment;
        // counts of reads assigned to them
        int hap1_count;
        int hap2_count;
        bool hap_updated;
        // to prevent flipping during haplotype reassignment if this is a seed unitig
        bool locked;
        // e.g., for DUP nodes
        int times_walked;
        int walk_limit;
        Unitig* parallel_utg;
};

bool compare_unitigs(Unitig* a, Unitig* b);

#endif
