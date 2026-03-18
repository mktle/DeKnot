#include "phase.h"
#include "unitig.h"
#include "utils.h"

#include <cassert>
#include <iostream>

Unitig::Unitig(std::vector<Node*>& set_nodes, int index){
    // deliberately copy
    this->nodes = set_nodes;
    this->index = index;
    this->base_freq = 0;
    this->update_median_freq();
    this->hap1_count = 0;
    this->hap2_count = 0;
}

bool Unitig::neighbors_unitig(){
    // check if either endpoint has a neighboring node that belongs to a different unitig
    if (this->nodes.front()->neighbors_diff_utg()){
        return true;
    }
    if (this->nodes.front() != this->nodes.back() && this->nodes.back()->neighbors_diff_utg()){
        return true;
    }
    return false;
}

bool Unitig::node_is_endpoint(Node* to_check){
    return (this->nodes.front()->get_node_id() == to_check->get_node_id() || this->nodes.back()->get_node_id() == to_check->get_node_id());
}

int Unitig::num_nodes(){
    return (int)this->nodes.size();
}

Node* Unitig::get_front(){
    return this->nodes.front();
}

Node* Unitig::get_back(){
    return this->nodes.back();
}

Node* Unitig::opposite_end(Node* in_end){
    if (in_end == this->nodes.back()){
        return this->nodes.front();
    }
    assert(in_end == this->nodes.front());
    return this->nodes.back();
}

int Unitig::get_id(){
    return this->index;
}

const std::vector<Node*>& Unitig::view_nodes(){
    return this->nodes;
}

Node* Unitig::node_at(size_t idx){
    return this->nodes.at(idx);
}

void Unitig::add_node_front(Node* to_add){
    this->nodes.insert(this->nodes.begin(), to_add);
}

void Unitig::add_node_back(Node* to_add){
    this->nodes.push_back(to_add);
}

bool Unitig::contains_node(Node* check){
    for (size_t i{0}; i < this->nodes.size(); i++){
        if (this->nodes[i]->get_node_id() == check->get_node_id()){
            return true;
        }
    }
    return false;
}

void Unitig::clear_nodes(){
    this->nodes.clear();
}

void Unitig::drain_nodes(std::vector<Node*>& recipient, Node* breakpoint){
    assert(this->nodes.front() != breakpoint);
    while (this->nodes.front() != breakpoint){
        recipient.push_back(this->nodes.front());
        this->nodes.erase(this->nodes.begin());
    }
    // figure out whether the breakpoint should be popped into the new unitig (depends on direction of the edge that was just added)
    int num_forward {0};
    char prev_dir = 'X';

    const std::vector<Edge> neighs = breakpoint->view_neighbors();
    for (size_t i{0}; i < neighs.size(); i++){
        if (neighs[i].neighbor == recipient.back()){
            prev_dir = neighs[i].self_orientation;
        }
        if (neighs[i].self_orientation == FORWARD_ORIENTATION){
            num_forward++;
        }
    }
    assert(prev_dir == FORWARD_ORIENTATION || prev_dir == REVERSE_ORIENTATION);
    if ((prev_dir == FORWARD_ORIENTATION && num_forward == 1) || (prev_dir == REVERSE_ORIENTATION && num_forward == 2)){
        recipient.push_back(this->nodes.front());
        this->nodes.erase(this->nodes.begin());
    }
}

Node* Unitig::pop_front(){
    Node* to_return = this->nodes.front();
    this->nodes.erase(this->nodes.begin());
    return to_return;
}

Node* Unitig::pop_back(){
    Node* to_return = this->nodes.back();
    this->nodes.pop_back();
    return to_return;
}

float Unitig::get_base_freq(){
    return this->base_freq;
}

float Unitig::get_median_freq(){
    return this->median_freq + this->base_freq;
}

void Unitig::update_median_freq(){
    std::vector<int> counts;
    for (size_t i{0}; i < this->num_nodes(); i++){
        int val {this->nodes[i]->get_count()};
        auto pos = std::lower_bound(counts.begin(), counts.end(), val);
        counts.insert(pos, val);
    }
    assert(counts.size() > 0);
    if (counts.size() % 2 == 0){
        this->median_freq = (float)(counts[counts.size() / 2] + counts[counts.size() / 2 - 1]) / (float)2.0;
    }else{
        this->median_freq = (float)counts[counts.size() / (float)2.0];
    }
}

void Unitig::set_base_freq(float set_freq){
    this->base_freq = set_freq;
}

// i.e., no neighbors
bool Unitig::is_floating(){
    return this->get_front()->is_tip() && this->get_back()->is_tip();
}

// a tip only has one neighbor, and that neighbor has another neighbor in the same direction
// return NULL if this is not a tip
Node* Unitig::get_tip_inner_end(Node*& adj_node){
    if (this->get_front()->is_tip()){
        // front doesn't have any other neighbors, so check the back's neighbor
        if (this->get_back()->num_neighbors() == 2 || (this->get_front() == this->get_back() && this->get_back()->num_neighbors() == 1)){
            const std::vector<Edge>& neighbors = this->get_back()->view_neighbors();
            for (size_t i{0}; i < neighbors.size(); i++){
                // get the neighbor on the different unitig
                if (neighbors[i].neighbor->get_unitig() != this){
                    assert(neighbors[i].neighbor->num_edges_in_dir(utils::reverse_orientation(neighbors[i].neighbor_orientation)) > 0);
                    if (neighbors[i].neighbor->num_edges_in_dir(utils::reverse_orientation(neighbors[i].neighbor_orientation)) > 1){
                        adj_node = neighbors[i].neighbor;
                        return this->get_back();
                    }
                }
            }
        }
    }else if (this->get_back()->is_tip()){
        // back doesn't have any other neighbors, so check the front's neighbor
        if (this->get_front()->num_neighbors() == 2){
            const std::vector<Edge>& neighbors = this->get_front()->view_neighbors();
            for (size_t i{0}; i < neighbors.size(); i++){
                // get the neighbor on the different unitig
                if (neighbors[i].neighbor->get_unitig() != this){
                    assert(neighbors[i].neighbor->num_edges_in_dir(utils::reverse_orientation(neighbors[i].neighbor_orientation)) > 0);
                    if (neighbors[i].neighbor->num_edges_in_dir(utils::reverse_orientation(neighbors[i].neighbor_orientation)) > 1){
                        adj_node = neighbors[i].neighbor;
                        return this->get_front();
                    }
                }
            }
        }
    }
    return NULL;
}

// get the direction of the endpoint that connects to the next node on this unitig
char Unitig::get_in_dir(Node* to_check){
    if (this->num_nodes() == 1){
        return INVALID_ORIENTATION;
    }
    assert(node_is_endpoint(to_check));
    const std::vector<Edge> neighs = to_check->view_neighbors();
    for (size_t i{0}; i < neighs.size(); i++){
        if (neighs[i].neighbor->get_unitig() == to_check->get_unitig()){
            return neighs[i].self_orientation;
        }
    }
    assert(0);
    return 'X';
}

bool Unitig::is_neighbors_with(Unitig* check_utg){
    return (this->get_front()->is_neighbors_with(check_utg->get_front())
            || this->get_front()->is_neighbors_with(check_utg->get_back())
            || this->get_back()->is_neighbors_with(check_utg->get_front())
            || this->get_back()->is_neighbors_with(check_utg->get_back())
            );
}

int Unitig::get_assignment(){
        return this->assignment;
}

void Unitig::set_assignment(int id){
        this->assignment = id;
}

int Unitig::get_hap_count(int id){
    assert(id == HAP1 || id == HAP2);
    if (id == HAP1) return this->hap1_count;
    return this->hap2_count;
}

Unitig* Unitig::get_parallel_utg(){
    std::vector<Unitig*> parallel_utgs;
    const std::vector<Edge> front_neighs = this->get_front()->view_neighbors();
    const std::vector<Edge> back_neighs = this->get_back()->view_neighbors();
    for (size_t i{0}; i < front_neighs.size(); i++){
        Node* neigh = front_neighs[i].neighbor;
        // get the direction of connection to unitig's front
        const std::vector<Edge> fnns = neigh->view_neighbors();
        char check_dir = INVALID_ORIENTATION;
        for (size_t j{0}; j < fnns.size(); j++){
            if (fnns[j].neighbor == this->get_front()){
                check_dir = fnns[j].self_orientation;
                break;
            }
        }
        assert(check_dir != INVALID_ORIENTATION);
        for (size_t j{0}; j < fnns.size(); j++){
            Node* parallel_node = fnns[j].neighbor;
            if (parallel_node != this->get_front() && fnns[j].self_orientation == check_dir){
                Unitig* parallel_utg = parallel_node->get_unitig();
                bool connected {false};
                for (size_t k{0}; k < back_neighs.size(); k++){
                    if (parallel_utg->opposite_end(parallel_node)->is_neighbors_with(back_neighs[k].neighbor)){
                        connected = true;
                        break;
                    }
                }
                if (connected) parallel_utgs.push_back(parallel_utg);
            }
        }
    }
    // only return if it is a simple bubble (i.e., exactly one parallel unitig)
    if (parallel_utgs.size() == 1){
        return parallel_utgs[0];
    }
    return NULL;
}

void Unitig::add_read(SketchedSeq* read){
    aligned_seqs.push_back(read);
}

void Unitig::assign_reads(){
    for (size_t i{0}; i < this->aligned_seqs.size(); i++){
        this->aligned_seqs[i]->assignment = this->get_assignment();
    }
}
void Unitig::tally_reads(){
    this->hap1_count = 0;
    this->hap2_count = 0;
    for (size_t i{0}; i < this->aligned_seqs.size(); i++){
        if (this->aligned_seqs[i]->assignment == HAP1) this->hap1_count++;
        if (this->aligned_seqs[i]->assignment == HAP2) this->hap2_count++;
    }
}

bool Unitig::reassign_haplotype(){
    int c1 = this->get_hap_count(HAP1);
    int c2 = this->get_hap_count(HAP2);
    int old_assign = this->get_assignment();
    if (c1 > 0 && c2 > 0){
        this->set_assignment(HOM);
    }else if (c1 > 0){
        this->set_assignment(HAP1);
    }else if (c2 > 0){
        this->set_assignment(HAP2);
    }else{
        this->set_assignment(UNASSIGNED);
    }
    return old_assign != this->get_assignment();
}

void Unitig::list_nodes(){
    for (size_t i{0}; i < this->nodes.size(); i++){
        std::cout << "(" << this->nodes[i]->get_node_id() << ") " << this->nodes[i]<< std::endl;
    }
}

void Unitig::get_hpc_sequence(std::string& ret_seq, bool reverse, bool append){
    if (this->nodes.size() == 1){
        assert(!append);
        this->get_front()->get_hpc_sequence(ret_seq, reverse);
    }else{
        Node* curr_node = reverse ? this->get_back() : this->get_front();
        int k_prime = curr_node->num_kmers();
        char curr_direction = INVALID_ORIENTATION;
        char next_direction = INVALID_ORIENTATION;
        // figure out the direction to read the initial sequence
        const std::vector<Edge>& neighbors = curr_node->view_neighbors();
        // add the initial sequence
        if (!append){
            for (size_t i{0}; i < neighbors.size(); i++){
                if (neighbors[i].neighbor->get_unitig()->get_id() == curr_node->get_unitig()->get_id()){
                    curr_direction = neighbors[i].self_orientation;
                    next_direction = neighbors[i].neighbor_orientation;
                    break;
                }
            }
            assert(curr_direction != INVALID_ORIENTATION);
            assert(next_direction != INVALID_ORIENTATION);
            curr_node->get_hpc_sequence(ret_seq, curr_direction == REVERSE_ORIENTATION);
        }else{
            // appending, so add just the last kmer from the first node, rather than the whole sequence
            for (size_t i{0}; i < neighbors.size(); i++){
                if (neighbors[i].neighbor->get_unitig()->get_id() == curr_node->get_unitig()->get_id()){
                    curr_direction = neighbors[i].self_orientation;
                    next_direction = neighbors[i].neighbor_orientation;
                }
            }
            assert(curr_direction != INVALID_ORIENTATION);
            assert(next_direction != INVALID_ORIENTATION);
            if (curr_direction == FORWARD_ORIENTATION){
                hash::add_to_seq_from_hash(curr_node->get_back_hash(), curr_node->ith_kmer_flipped(k_prime - 1), curr_node->ith_unique_seq_len(k_prime - 2), ret_seq);
            }else{
                assert(curr_direction == REVERSE_ORIENTATION);
                hash::add_to_seq_from_hash(curr_node->get_front_hash(), !curr_node->ith_kmer_flipped(0), curr_node->ith_unique_seq_len(0), ret_seq);
            }
        }
        // now go through the rest of the nodes in the unitig and add the last kmer of that node
        if(!reverse){
            for (size_t i{1}; i < this->nodes.size(); i++){
                curr_node = this->nodes[i];
                // add the kmer that extends past the current node
                if (next_direction == FORWARD_ORIENTATION){
                    hash::add_to_seq_from_hash(curr_node->get_back_hash(), curr_node->ith_kmer_flipped(k_prime - 1), curr_node->ith_unique_seq_len(k_prime - 2), ret_seq);
                }else{
                    assert(next_direction == REVERSE_ORIENTATION);
                    hash::add_to_seq_from_hash(curr_node->get_front_hash(), !curr_node->ith_kmer_flipped(0), curr_node->ith_unique_seq_len(0), ret_seq);
                }

                if (i < this->nodes.size() - 1){
                    char save_direction = next_direction;
                    next_direction = INVALID_ORIENTATION;
                    // figure out which node is the next one
                    const std::vector<Edge>& curr_neighbors = curr_node->view_neighbors();
                    for (size_t j{0}; j < curr_neighbors.size(); j++){
                        Node* neigh = curr_neighbors[j].neighbor;
                        if (neigh != this->nodes[i - 1] && neigh->get_unitig()->get_id() == curr_node->get_unitig()->get_id()){
                            if (curr_neighbors[j].self_orientation != save_direction){
                                next_direction = utils::reverse_orientation(curr_neighbors[j].neighbor_orientation);
                            }else{
                                next_direction = curr_neighbors[j].neighbor_orientation;
                            }
                        }
                    }
                    assert(next_direction != INVALID_ORIENTATION);
                }
            }
        }else{
            int i = (int)this->nodes.size() - 2;
            for (; i >= 0; i--){
                curr_node = this->nodes[(size_t)i];
                // add the kmer that extends past the current node
                if (next_direction == FORWARD_ORIENTATION){
                    hash::add_to_seq_from_hash(curr_node->get_back_hash(), curr_node->ith_kmer_flipped(k_prime - 1), curr_node->ith_unique_seq_len(k_prime - 2), ret_seq);
                }else{
                    assert(next_direction == REVERSE_ORIENTATION);
                    hash::add_to_seq_from_hash(curr_node->get_front_hash(), !curr_node->ith_kmer_flipped(0), curr_node->ith_unique_seq_len(0), ret_seq);
                }

                if (i > 0){
                    char save_direction = next_direction;
                    next_direction = INVALID_ORIENTATION;
                    // figure out which node is the next one
                    const std::vector<Edge>& curr_neighbors = curr_node->view_neighbors();
                    for (size_t j{0}; j < curr_neighbors.size(); j++){
                        Node* neigh = curr_neighbors[j].neighbor;
                        if (neigh != this->nodes[(size_t)i + 1] && neigh->get_unitig()->get_id() == curr_node->get_unitig()->get_id()){
                            if (curr_neighbors[j].self_orientation != save_direction){
                                next_direction = utils::reverse_orientation(curr_neighbors[j].neighbor_orientation);
                            }else{
                                next_direction = curr_neighbors[j].neighbor_orientation;
                            }
                        }
                    }
                    assert(next_direction != INVALID_ORIENTATION);
                }
            }
        }
    }
}

// for debugging
// prints the unique kmers of each node
void Unitig::print_kmers(){
    std::cout << "kmers for utg" << this->get_id() << std::endl;
    Node* curr_node = this->get_front();
    curr_node->print_kmers();
    if (this->get_front() == this->get_back()) return;
    int k_prime = curr_node->num_kmers();
    char next_direction = INVALID_ORIENTATION;
    // figure out the direction to read the initial sequence
    const std::vector<Edge>& neighbors = curr_node->view_neighbors();
    for (size_t i{0}; i < neighbors.size(); i++){
        if (neighbors[i].neighbor->get_unitig()->get_id() == curr_node->get_unitig()->get_id()){
            next_direction = neighbors[i].neighbor_orientation;
            break;
        }
    }
    assert(next_direction != INVALID_ORIENTATION);
    size_t offset {(size_t)k_prime - 1};

    // now go through the rest of the nodes in the unitig and add the last kmer of that node
    for (size_t i{1}; i < this->nodes.size(); i++){
        curr_node = this->nodes[i];
        // add the kmer that extends past the current node
        if (next_direction == FORWARD_ORIENTATION){
            std::string res;
            hash::hash_to_kmer(curr_node->ith_kmer(k_prime - 1), curr_node->ith_kmer_flipped(k_prime - 1), res);
            std::cout << "(" << i + offset << ") " << curr_node->ith_kmer(k_prime - 1) << " " << res << " " << curr_node->ith_unique_seq_len(k_prime - 2) << std::endl;
        }else{
            assert(next_direction == REVERSE_ORIENTATION);
            std::string res;
            hash::hash_to_kmer(curr_node->ith_kmer(0), curr_node->ith_kmer_flipped(0), res);
            std::cout << "(" << i + offset << ") " << curr_node->ith_kmer(0) << " " << res << " " << curr_node->ith_unique_seq_len(0) << std::endl;
        }

        if (i < this->nodes.size() - 1){
            char save_direction = next_direction;
            next_direction = INVALID_ORIENTATION;
            // figure out which node is the next one
            const std::vector<Edge>& curr_neighbors = curr_node->view_neighbors();
            for (size_t j{0}; j < curr_neighbors.size(); j++){
                Node* neigh = curr_neighbors[j].neighbor;
                if (neigh != this->nodes[i - 1] && neigh->get_unitig()->get_id() == curr_node->get_unitig()->get_id()){
                    if (curr_neighbors[j].self_orientation != save_direction){
                        next_direction = utils::reverse_orientation(curr_neighbors[j].neighbor_orientation);
                    }else{
                        next_direction = curr_neighbors[j].neighbor_orientation;
                    }
                }
            }
            assert(next_direction != INVALID_ORIENTATION);
        }
    }
    std::cout << "----------" << std::endl;
}


void Unitig::print(bool verbose){
    this->update_median_freq();
    std::cout << "[" << this->get_id() << "][" << this->get_assignment() << "]";
    std::cout << " n" << this->get_front()->get_node_id() << "\tn" << this->get_back()->get_node_id() << "\tmedian freq " << this->get_median_freq() << '\n';
    if (verbose){
        this->list_nodes();
    }
}

// return true if a > b
bool compare_unitigs(Unitig* a, Unitig* b){
    // first comparison: length of unitig
    if (a->num_nodes() > b->num_nodes()){
        return true;
    }else if (b->num_nodes() > a->num_nodes()){
        return false;
    }

    // second comparison: median frequency of unitig
    if (a->get_median_freq() > b->get_median_freq()){
        return true;
    }else if (b->get_median_freq() > a->get_median_freq()){
        return false;
    }

    // third comparison: arbitrary; utg ID
    if (a->get_id() < b->get_id()){
        return true;
    }
    return false;
}
