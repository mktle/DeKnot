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
        int get_id();

        const std::vector<Node*>& view_nodes();
        Node* node_at(size_t idx);

        void add_node_front(Node* to_add);
        void add_node_back(Node* to_add);

        void drain_nodes(std::vector<Node*>& recipient, Node* breakpoint);
        void clear_nodes();
        Node* pop_front();
        Node* pop_back();

        float get_median_freq();
        float get_base_freq();
        void update_median_freq();
        void set_base_freq(float set_freq);
        bool is_floating();
        Node* get_tip_inner_end(Node*& adj_node);
        char get_in_dir(Node* node_to_check);

        void list_nodes();
        void get_hpc_sequence(std::string& ret_seq);
        void print_kmers();
        void print(bool verbose);

    private:
        std::vector<Node*> nodes;
        int index;
        float median_freq;
        float base_freq;
};

bool compare_unitigs(Unitig* a, Unitig* b);

#endif
