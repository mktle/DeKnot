#ifndef UNITIG_H
#define UNITIG_H

#include "hash.h"
#include "node.h"

#include <algorithm>
#include <iterator>
#include <vector>

class Node;

// node list of unitig
class NodeBuf{
    public:
        using iterator = Node**;
        using const_iterator = Node* const*;

        NodeBuf() : head(0) {}
        NodeBuf& operator=(const std::vector<Node*>& v){
            this->buf = v;
            this->head = 0;
            return *this;
        }
        size_t size() const { return this->buf.size() - this->head; }
        bool empty() const { return this->size() == 0; }
        Node*& operator[](size_t i) { return this->buf[this->head + i]; }
        Node* const& operator[](size_t i) const { return this->buf[this->head + i]; }
        Node*& at(size_t i) { return this->buf.at(this->head + i); }
        Node*& front() { return this->buf[this->head]; }
        Node*& back() { return this->buf.back(); }
        Node* front() const { return this->buf[this->head]; }
        Node* back() const { return this->buf.back(); }
        iterator begin() { return this->buf.data() + this->head; }
        iterator end() { return this->buf.data() + this->buf.size(); }
        const_iterator begin() const { return this->buf.data() + this->head; }
        const_iterator end() const { return this->buf.data() + this->buf.size(); }
        std::reverse_iterator<const_iterator> rbegin() const { return std::reverse_iterator<const_iterator>(this->end()); }
        std::reverse_iterator<const_iterator> rend() const { return std::reverse_iterator<const_iterator>(this->begin()); }
        void clear(){
            this->buf.clear();
            this->head = 0;
        }
        void push_back(Node* x) { this->buf.push_back(x); }
        void pop_back() { this->buf.pop_back(); }
        template <typename It>
        void insert(iterator pos, It first, It last){
            size_t m = (size_t)std::distance(first, last);
            if (m == 0) return;
            if (pos == this->begin()){
                if (this->head < m){
                    // make room at the front
                    size_t n = this->size();
                    size_t new_head = std::max(m, n) + 16;
                    std::vector<Node*> grown(new_head + n);
                    std::copy(this->begin(), this->end(), grown.begin() + (long)new_head);
                    this->buf.swap(grown);
                    this->head = new_head;
                }
                this->head -= m;
                std::copy(first, last, this->buf.begin() + (long)this->head);
            }else if (pos == this->end()){
                this->buf.insert(this->buf.end(), first, last);
            }else{
                this->buf.insert(this->buf.begin() + (pos - this->buf.data()), first, last);
            }
        }
        void insert(iterator pos, Node* x){
            Node* one[1] = {x};
            this->insert(pos, one, one + 1);
        }
        void erase(iterator first, iterator last){
            if (first == this->begin()){
                this->head += (size_t)(last - first);
                if (this->head == this->buf.size()) this->clear();
            }else{
                this->buf.erase(this->buf.begin() + (first - this->buf.data()), this->buf.begin() + (last - this->buf.data()));
            }
        }
        void erase(iterator pos) { this->erase(pos, pos + 1); }

    private:
        std::vector<Node*> buf;
        size_t head;
};

class Unitig{
    public:
        Unitig(std::vector<Node*>& set_nodes, int index);
        void reset(std::vector<Node*>& set_nodes, int index);
        void reset(Node* single_node, int index);

        bool neighbors_unitig();
        bool node_is_endpoint(Node* to_check);
        int num_nodes() { return (int)this->nodes.size(); }
        Node* get_front() { return this->nodes.front(); }
        Node* get_back() { return this->nodes.back(); }
        Node* opposite_end(Node* in_end);
        int get_id() { return this->index; }

        const NodeBuf& view_nodes() { return this->nodes; }
        Node* node_at(size_t idx);

        void add_node_front(Node* to_add);
        void add_node_back(Node* to_add);
        void prepend_nodes(const NodeBuf& to_add, bool reversed);
        void append_nodes(const NodeBuf& to_add, bool reversed);
        bool contains_node(Node* check);

        void drain_nodes(std::vector<Node*>& recipient, Node* breakpoint);
        void clear_nodes();
        Node* pop_front();
        Node* pop_back();

        float get_median_freq() { return this->median_freq + this->base_freq; }
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
            
        int get_assignment() { return this->assignment; }
        void set_assignment(int id);
        void set_hap_count(int id, int count);
        int get_hap_count(int id);
        float get_hap_frac(int id);
        void set_hap_updated(bool set);
        bool get_hap_updated();
        void set_component(int id);
        int get_component() { return this->component; }

        void set_parallel_utg();
        Unitig* get_parallel_utg();
        void add_read(SketchedSeq*);
        void assign_reads(bool lock);
        void tally_reads();
        void lock_haplotype();
        bool is_locked() { return this->locked; }
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
        NodeBuf nodes;
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
        // connected component assignment
        int component;
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
        void init_fields(int index);
};

bool compare_unitigs(Unitig* a, Unitig* b);

#endif
