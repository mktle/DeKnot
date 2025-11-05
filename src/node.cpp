#include "node.h"

#include "hash.h"
#include "utils.h"

#include <cassert>
#include <iostream>

//15 length substrings of "ACACACACACACACACACACACACACACACA", "ATATATATATATATATATATATATATATATA", "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC", their reverse complements, etc.
const ankerl::unordered_dense::set<uint64_t> DINUCLEOTIDE_SUBHASHES = {71582788ULL, 1002159035ULL, 143165576ULL,214748364ULL, 286331153ULL, 501079517ULL, 572662306ULL, 787410670ULL, 858993459ULL, 930576247ULL};

Node::Node(KmerView& set_ksyncmer, Span<uint8_t>& set_flipped, Span<int>& set_unique_len, int set_node, Unitig* set_unitig) : ksyncmer(set_ksyncmer), flip_kmer(set_flipped), unique_seq_len(set_unique_len) {
    this->node_idx = set_node;
    this->utg_ptr = set_unitig;
    this->count = 1;
    this->num_forward_edges = 0;
    this->num_reverse_edges = 0;
}

void Node::set_count(int n){
    this->count = n;
}

void Node::add_neighbor(Node* new_neighbor, char curr_orient, char neigh_orient) {
    this->neighbors.emplace_back(new_neighbor, curr_orient, neigh_orient);
    if (curr_orient == FORWARD_ORIENTATION){
        this->num_forward_edges++;
    }else{
        assert(curr_orient == REVERSE_ORIENTATION);
        this->num_reverse_edges++;
    }
}

void Node::remove_neighbor(Node* remove_neighbor){
    for (auto it = this->neighbors.begin(); it != this->neighbors.end(); it++) {
        if (it->neighbor == remove_neighbor) {
            if (it->self_orientation == FORWARD_ORIENTATION){
                this->num_forward_edges--;
                assert(this->num_forward_edges >= 0);
            }else{
                assert(it->self_orientation == REVERSE_ORIENTATION);
                this->num_reverse_edges--;
            }
            it = this->neighbors.erase(it);
            return;
        } 
    }
    assert(0);
}

int Node::num_edges_in_dir(char check_dir){
    if (check_dir == FORWARD_ORIENTATION){
        return this->num_forward_edges;
    }else if(check_dir == REVERSE_ORIENTATION){
        return this->num_reverse_edges;
    }else{
        assert(check_dir == INVALID_ORIENTATION);
        return 0;
    }
}

void Node::assign_unitig(Unitig* new_ptr) {
    this->utg_ptr = new_ptr;
}

bool Node::is_assigned_to_unitig() {
    return this->utg_ptr != NULL;
}

bool Node::is_neighbors_with(Node* to_check){
    for (int i{0}; i < this->neighbors.size(); i++){
        if (this->neighbors[(size_t)i].neighbor == to_check){
            return true;
        }
    }
    return false;
}

uint64_t Node::unique_kmer_from_neigh(Node* neigh, bool& flipped, int& unique_seq_len){
    for (int i{0}; i < this->neighbors.size(); i++){
        if (this->neighbors[(size_t)i].neighbor == neigh){
            if (this->neighbors[(size_t)i].self_orientation == FORWARD_ORIENTATION){
                flipped = this->ith_kmer_flipped(this->num_kmers() - 1);
                unique_seq_len = this->ith_unique_seq_len(this->num_kmers() - 2);
                return this->get_back_hash();
            }
            flipped = !this->ith_kmer_flipped(0);
            unique_seq_len = this->ith_unique_seq_len(0);
            return this->get_front_hash();
        }
    }
    assert(0);
}

int Node::num_kmers(){
    return (int)this->ksyncmer.kmers.size();
}

int Node::num_neighbors(){
    return (int)this->neighbors.size();
}

int Node::get_node_id(){
    return this->node_idx;
}

Unitig* Node::get_unitig(){
    return this->utg_ptr;
}

int Node::get_count(){
    return this->count;
}

uint64_t Node::ith_kmer(int i){
    if (this->ksyncmer.flip){
        return this->ksyncmer.kmers[this->ksyncmer.kmers.size() - 1 - (size_t)i];
    }
    return this->ksyncmer.kmers[(size_t)i];
}

bool Node::ith_kmer_flipped(int i){
    assert(i >= 0 && i < this->ksyncmer.kmers.size());
    if (this->ksyncmer.flip){
        return !this->flip_kmer[this->flip_kmer.size() - 1 - (size_t)i];
    }
    return this->flip_kmer[(size_t)i];
}

int Node::ith_unique_seq_len(int i){
    if (this->ksyncmer.flip){
        return this->unique_seq_len[this->unique_seq_len.size() - 1 - (size_t)i];
    }
    return this->unique_seq_len[(size_t)i];
}

uint64_t Node::get_front_hash(){
    if (this->ksyncmer.flip){
        return *(this->ksyncmer.kmers.end() - 1);
    }
    return *(this->ksyncmer.kmers.begin());
}

uint64_t Node::get_back_hash(){
    if (this->ksyncmer.flip){
        return *(this->ksyncmer.kmers.begin());
    }
    return *(this->ksyncmer.kmers.end() - 1);
}

bool Node::neighbors_diff_utg(){
    for (size_t i{0}; i < neighbors.size(); i++){
        if (neighbors[i].neighbor->get_unitig() != this->get_unitig()){
            return true;
        }
    }
    return false;
}

int Node::longest_neighbor_utg_length(){
    int longest_length {0};
    for (size_t i{0}; i < neighbors.size(); i++){
        Unitig* neigh_utg = neighbors[i].neighbor->get_unitig();
        if (neigh_utg != this->get_unitig() && neigh_utg->num_nodes() > longest_length){
            longest_length = neigh_utg->num_nodes();
        }
    }
    return longest_length;
}

KmerView Node::get_kmer_view(){
    return this->ksyncmer;
}

const std::vector<Edge>& Node::view_neighbors(){
    return this->neighbors;
}

bool Node::is_tip(){
    int num_forward {0};
    int num_reverse {0};
    for (size_t i{0}; i < this->neighbors.size(); i++){
        if (neighbors[i].self_orientation == FORWARD_ORIENTATION){
            num_forward++;
        }else{
            assert(neighbors[i].self_orientation == REVERSE_ORIENTATION);
            num_reverse++;
        }
    }
    return (num_forward == 0 || num_reverse == 0);
}

void Node::get_hpc_sequence(std::string& ret_str, bool reverse){
    if (reverse){
        int i = (int)this->ksyncmer.kmers.size() - 1;
        hash::hash_to_kmer(this->ith_kmer(i), !this->ith_kmer_flipped(i), ret_str);
        i--;
        for (; i >= 0; i--){
            hash::add_to_seq_from_hash(this->ith_kmer(i), !this->ith_kmer_flipped(i), this->ith_unique_seq_len(i), ret_str);
        }
    }else{
        hash::hash_to_kmer(this->ith_kmer(0), this->ith_kmer_flipped(0), ret_str);
        for (int i{1}; i < this->ksyncmer.kmers.size(); i++){
            hash::add_to_seq_from_hash(this->ith_kmer(i), this->ith_kmer_flipped(i), this->ith_unique_seq_len(i - 1), ret_str);
        }
    }
}

// returns true if a dinucleotide hash is in the first or last kmer of this node
// checking for substrings of size 15
bool Node::end_is_dinucleotide(){
    int last_idx = this->num_kmers() - 1;
    for (int side : {0, last_idx}){
        uint64_t hash = this->ksyncmer.kmers[(size_t)side];
        uint64_t least_sig = (hash >> 2) & ((1ULL << 30) - 1);
        uint64_t most_sig = hash >> (64 - 30); 
        for (uint64_t subhash : {least_sig, most_sig}){
            if (DINUCLEOTIDE_SUBHASHES.contains(subhash)){
                return true;
            }
        }
    }
    return false;
}

void Node::print_kmers(){
    if (this->ksyncmer.flip){
        for (int i{(int)this->ksyncmer.kmers.size() - 1}; i >= 0 ; i--){
            std::cout << this->ksyncmer.kmers[(size_t)i] << std::endl;
        }
    }else{
        for (size_t i{0}; i < this->ksyncmer.kmers.size(); i++){
            std::cout << "(" << i << ")" << this->ksyncmer.kmers[(size_t)i] << std::endl;
        }
    }
}

void Node::print(bool verbose){
    std::cout << "n" << this->get_node_id() << " (utg" << this->get_unitig()->get_id() << ")" << std::endl;
    if (verbose){
        std::string seq {""};
        this->get_hpc_sequence(seq, false);
        std::cout << seq << std::endl;
    }
    for (size_t i{0}; i < this->neighbors.size(); i++){
        std::cout << "\tn" << this->neighbors[i].neighbor->get_node_id() << " (utg" << this->neighbors[i].neighbor->get_unitig()->get_id() << ") " << this->neighbors[i].self_orientation << " " << this->neighbors[i].neighbor_orientation << std::endl;
    }
}
