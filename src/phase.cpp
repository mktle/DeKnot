#include "phase.h"

#include <cassert>
#include <iostream>
#include <unordered_set>
#include <cstdlib>
#include <cmath>
#include <algorithm>

void phase::initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs, int comp_num){
    graph->update_median_coverage();
    if (!graph->aligned_reads()) graph->align_reads_to_unitigs(reads);
    const std::vector<Unitig*> unitig_list = graph->view_unitigs();
    for (size_t i{0}; i < unitig_list.size(); i++){
        Unitig* utg = unitig_list[i];
        if (utg != NULL && utg->get_component() == comp_num){
            utg->set_assignment(UNASSIGNED);
            utg->set_parallel_utg();
            sorted_unitigs.push_back(utg);
        }
    }
    std::sort(sorted_unitigs.begin(), sorted_unitigs.end(), compare_unitigs);
}

bool phase::seed_initial_haplotypes(std::vector<Unitig*>& sorted_unitigs){
    /***************************************
     * Choose bubble to start phasing from *
     ***************************************/
    Unitig* initial_utg = NULL;
    Unitig* parallel_utg = NULL;
    for (size_t i{0}; i < sorted_unitigs.size(); i++){
        parallel_utg = sorted_unitigs[i]->get_parallel_utg();
        if (parallel_utg != NULL){
            initial_utg = sorted_unitigs[i];
            break;
        }
    }
    if (initial_utg == NULL) return false;
    initial_utg->set_assignment(HAP1);
    initial_utg->assign_reads(true);
    initial_utg->lock_haplotype();
    initial_utg->tally_reads();
    parallel_utg->set_assignment(HAP2);
    parallel_utg->assign_reads(true);
    parallel_utg->lock_haplotype();
    parallel_utg->tally_reads();
    return true;
}

// return true if any assignments changed
bool phase::reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& utg_list){
    std::vector<int> before(utg_list.size());
    for (size_t i{0}; i < utg_list.size(); i++){
        utg_list[i]->set_hap_updated(false);
        before[i] = utg_list[i]->get_assignment();
    }

    /**************************
     * Reassign unitig labels *
     **************************/
    bool changed {false};
    for (size_t i{0}; i < utg_list.size(); i++){
        graph->reassign_haplotype(utg_list[i]);
    }
    for (size_t i{0}; i < utg_list.size(); i++){
        if (utg_list[i]->get_assignment() != before[i]){
            changed = true;
        }
    }

    /************************
     * Reassign read labels *
     ************************/
    for (size_t i{0}; i < reads.size(); i++){
        if (reads[i]->is_read && !reads[i]->locked && reads[i]->kmers.size() >= (size_t)graph->get_k_prime()){
            int num1 {0};
            int num2 {0};
            size_t check_thresh = reads[i]->kmers.size() - (size_t)graph->get_k_prime() + 1;
            utils::RollingHasher rh(reads[i]->kmers.data(), reads[i]->kmers.size(), (size_t)graph->get_k_prime());
            for (size_t j{0}; j < check_thresh; j++){
                Span<uint64_t> curr_kmers(reads[i]->kmers.data() + j, (size_t)graph->get_k_prime());
                char curr_orient = syncmer::canonical_syncmer(curr_kmers);
                KmerView curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION, rh.hash_at(j, curr_orient == REVERSE_ORIENTATION)};
                Node* curr_node = graph->get_node_from_lookup_table(curr_canon);
                //if (curr_node != NULL && curr_node->get_unitig()->get_parallel_utg() != NULL){
                if (curr_node != NULL){
                    if (curr_node->get_unitig()->get_assignment() == HAP1) num1++;
                    if (curr_node->get_unitig()->get_assignment() == HAP2) num2++;
                }
            }
            if (num1 > num2 && (reads[i]->assignment == HAP2 || reads[i]->assignment == UNASSIGNED)){
                reads[i]->assignment = HAP1;
                changed =  true;
            }else if (num2 > num1 && (reads[i]->assignment == HAP1 || reads[i]->assignment == UNASSIGNED)){
                reads[i]->assignment = HAP2;
                changed = true;
            }
        }
    }
    return changed;
}

void phase::remove_tagged_unitigs(std::vector<Unitig*>& utg_list){
    auto assigned = [](Unitig* utg) { return utg->get_assignment() != UNASSIGNED; };
    auto new_end = std::remove_if(utg_list.begin(), utg_list.end(), assigned);
    utg_list.erase(new_end, utg_list.end());
}

void phase::identify_endpoints(Graph* graph, std::vector<std::vector<Unitig*>>& hap_endpoints, std::vector<Unitig*>& hom_endpoints, int comp_num){
    hap_endpoints.resize(2);
    auto& hap1_endpoints = hap_endpoints[0];
    auto& hap2_endpoints = hap_endpoints[1];

    const std::vector<Unitig*> utgs = graph->view_unitigs();
    for (size_t i{0}; i < utgs.size(); i++){
        if(utgs[i] != NULL && utgs[i]->get_component() == comp_num){
            int front_neighs, back_neighs;
            if (utgs[i]->get_front() == utgs[i]->get_back()){
                int n_neighs {utgs[i]->get_front()->num_neighbors()};
                front_neighs = (n_neighs == 1) ? 0 : n_neighs;
                back_neighs = n_neighs;
            }else{
                // subtract one for the neighbor on the same unitig
                front_neighs = utgs[i]->get_front()->num_neighbors() - 1;
                back_neighs = utgs[i]->get_back()->num_neighbors() - 1;
            }

            if ((front_neighs == 0 && back_neighs > 0) || (back_neighs == 0 && front_neighs > 0)){
                if (utgs[i]->get_assignment() == HOM){
                    hom_endpoints.push_back(utgs[i]);
                }else if (utgs[i]->get_assignment() == HAP1){
                    hap1_endpoints.push_back(utgs[i]);
                }else if (utgs[i]->get_assignment() == HAP2){
                    hap2_endpoints.push_back(utgs[i]);
                }
            }
        }
    }
    std::sort(hap1_endpoints.begin(), hap1_endpoints.end(), compare_unitigs);
    std::sort(hap2_endpoints.begin(), hap2_endpoints.end(), compare_unitigs);
    std::sort(hom_endpoints.begin(), hom_endpoints.end(), compare_unitigs);
}

bool phase::get_haplotype_sequences(Graph* graph, std::vector<std::string>& ret_seqs, int comp_num, std::vector<Unitig*>* walked){
    std::vector<std::vector<Unitig*>> hap_endpoints;
    std::vector<Unitig*> hom_endpoints;
    phase::identify_endpoints(graph, hap_endpoints, hom_endpoints, comp_num);
    if(hap_endpoints[0].size() == 0 && hap_endpoints[1].size() == 0 && hom_endpoints.size() == 0){
        return false;
    }

    /**********************
     * Get walk sequences *
     **********************/
    Node* start_node = NULL;
    graph->update_median_coverage();
     for (size_t i = 0; i < 2; i++) {
        graph->reset_haplotypes();
        auto& endpoints = hap_endpoints[i];
        auto& other_endpoints = hap_endpoints[1 - i];

        int hap_num = (i == 0) ? HAP1 : HAP2;
        // if none in this hap, both ends must be homozygous
        if (endpoints.size() > 0){
            start_node = endpoints[0]->get_front()->neighbors_diff_utg() ? endpoints[0]->get_back() : endpoints[0]->get_front();
        }else if (hom_endpoints.size() > 0){
            start_node = hom_endpoints[0]->get_front()->neighbors_diff_utg() ? hom_endpoints[0]->get_back() : hom_endpoints[0]->get_front();
        }else{
            // navigate to first
            Unitig* start_utg = other_endpoints[other_endpoints.size() - 1];
            Node* first_node = start_utg->get_front()->neighbors_diff_utg() ? start_utg->get_front() : start_utg->get_back();
            const EdgeList& neighs = first_node->view_neighbors();
            start_node = NULL;
            for (size_t j{0}; j < neighs.size(); j++){
                Node* neigh = neighs[j].neighbor;
                Unitig* neigh_utg = neigh->get_unitig();
                if (neigh_utg != start_utg && (neigh_utg->get_assignment() == hap_num || (neigh_utg->get_assignment() == HOM && start_node == NULL))){
                    start_node = neigh;
                    if (neigh_utg->get_assignment() == hap_num) break;
                }
            }
            if (start_node == NULL){
                start_node = first_node;
            }
        }
        {
            std::vector<Node*> cands;
            for (Unitig* u : endpoints) cands.push_back(u->get_front()->neighbors_diff_utg() ? u->get_back() : u->get_front());
            for (Unitig* u : hom_endpoints) cands.push_back(u->get_front()->neighbors_diff_utg() ? u->get_back() : u->get_front());
            cands.push_back(start_node);
            long best_score {-1};
            Node* best_start = start_node;
            for (Node* cand : cands){
                graph->reset_haplotypes();
                Walk trial(NULL, 0);
                phase::build_walk(trial, cand, hap_num);
                long hap_nodes {0}, total_nodes {0};
                for (Unitig* u : trial.utgs){
                    total_nodes += u->num_nodes();
                    if (u->get_assignment() == hap_num) hap_nodes += u->num_nodes();
                }
                long score = hap_nodes * 1000000L + total_nodes;
                if (score > best_score){ best_score = score; best_start = cand; }
            }
            graph->reset_haplotypes();
            start_node = best_start;
        }
        Walk best_walk(NULL, 0);
        phase::build_walk(best_walk, start_node, hap_num);
        if (getenv("DK_THREAD_DUMP")){
            std::cout << "WALKPATH hap" << hap_num << ":";
            for (Unitig* u : best_walk.utgs) std::cout << " utg" << u->get_id();
            std::cout << std::endl;
        }
        if (walked != NULL) walked->insert(walked->end(), best_walk.utgs.begin(), best_walk.utgs.end());

        std::string walk_seq;
        assert(best_walk.utgs[0] == start_node->get_unitig());
        auto single_node_reverse = [](Node* node, Node* other, bool entering) -> bool {
            const EdgeList& edges = entering ? other->view_neighbors() : node->view_neighbors();
            Node* target = entering ? node : other;
            for (size_t m{0}; m < edges.size(); m++){
                if (edges[m].neighbor == target){
                    return (entering ? edges[m].neighbor_orientation : edges[m].self_orientation) == REVERSE_ORIENTATION;
                }
            }
            return false;
        };
        Unitig* first_utg = start_node->get_unitig();
        bool first_reverse = (start_node == first_utg->get_back());
        if (first_utg->get_front() == first_utg->get_back() && best_walk.utgs.size() > 1){
            Unitig* second = best_walk.utgs[1];
            Node* second_entry = second->get_front()->is_neighbors_with(start_node) ? second->get_front() : second->get_back();
            first_reverse = single_node_reverse(start_node, second_entry, false);
        }
        first_utg->get_hpc_sequence(walk_seq, first_reverse, false);

        Node* exit_node = first_utg->opposite_end(start_node);
        for (size_t j{1}; j < best_walk.utgs.size(); j++){
            Unitig* next_utg = best_walk.utgs[j];
            Node* next_node = (next_utg->get_front()->is_neighbors_with(exit_node)) ? next_utg->get_front() : next_utg->get_back();
            bool next_reverse = (next_utg->get_front() == next_utg->get_back()) ? single_node_reverse(next_node, exit_node, true) : (next_node == next_utg->get_back());
            next_utg->get_hpc_sequence(walk_seq, next_reverse, true);
            exit_node = next_utg->opposite_end(next_node);
        }
        ret_seqs.emplace_back(walk_seq);
    }
    return true;
}


namespace {
    const size_t THREAD_W = 16;
    const int THREAD_MIN_RUN = 3;
    const size_t THREAD_MAX_SUFFIX = 8;
    const int THREAD_MIN_READS = 2;
    const double THREAD_MIN_FRAC = 0.2;
    const int THREAD_SIG_DEPTH = 4;
    const double THREAD_REVISIT_FRAC = 0.1;
    const double THREAD_DETECT_FRAC = 0.3;
    struct ThreadData{
        bool active {false};
        std::vector<std::vector<int>> paths;
        ankerl::unordered_dense::set<int> placeable;
    };
    thread_local ThreadData g_threads;
    thread_local int g_sketch_s {0};

    uint64_t mix64(uint64_t x){
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    uint64_t window_hash(const uint64_t* k, size_t w){
        uint64_t f {0}, r {0};
        for (size_t i{0}; i < w; i++){
            f = mix64(f ^ k[i]);
            r = mix64(r ^ k[w - 1 - i]);
        }
        return f < r ? f : r;
    }
}

void phase::set_sketch_s(int s){ g_sketch_s = s; }

static void build_threads(Graph* graph, std::vector<SketchedSeq*>& reads, int comp_num, bool verbose){
    g_threads.paths.clear();
    g_threads.active = false;
    if (g_sketch_s <= 0 || getenv("DK_NO_THREAD") != NULL) return;
    ankerl::unordered_dense::map<uint64_t, int> owner;
    const std::vector<Unitig*> utgs = graph->view_unitigs();
    for (Unitig* u : utgs){
        if (u == NULL || u->get_component() != comp_num) continue;
        std::string seq;
        u->get_hpc_sequence(seq, false, false);
        std::vector<uint64_t> kmers; std::vector<uint8_t> flipped; std::vector<int> uniq, idx;
        syncmer::get_all_syncmers(seq, g_sketch_s, kmers, flipped, uniq, idx);
        if (kmers.size() < THREAD_W) continue;
        for (size_t j{0}; j + THREAD_W <= kmers.size(); j++){
            uint64_t h = window_hash(kmers.data() + j, THREAD_W);
            auto it = owner.find(h);
            if (it == owner.end()) owner.emplace(h, u->get_id());
            else it->second = -1;
        }
    }
    g_threads.placeable.clear();
    {
        ankerl::unordered_dense::map<int, int> n_unique;
        for (auto& kv : owner) if (kv.second >= 0) n_unique[kv.second]++;
        for (auto& kv : n_unique) if (kv.second >= THREAD_MIN_RUN) g_threads.placeable.insert(kv.first);
    }
    size_t n_multi {0};
    for (SketchedSeq* r : reads){
        if (!r->is_read || r->kmers.size() < THREAD_W) continue;
        std::vector<std::pair<int,int>> runs;
        for (size_t j{0}; j + THREAD_W <= r->kmers.size(); j++){
            auto it = owner.find(window_hash(r->kmers.data() + j, THREAD_W));
            if (it == owner.end() || it->second < 0) continue;
            if (!runs.empty() && runs.back().first == it->second) runs.back().second++;
            else runs.emplace_back(it->second, 1);
        }
        std::vector<int> path;
        for (auto& run : runs){
            if (run.second < THREAD_MIN_RUN) continue;
            if (path.empty() || path.back() != run.first) path.push_back(run.first);
        }
        if (!path.empty()) g_threads.paths.push_back(std::move(path));
    }
    ankerl::unordered_dense::map<int, int> n_reg, n_revisit;
    for (auto& p : g_threads.paths){
        ankerl::unordered_dense::map<int, int> seen;
        for (int x : p) seen[x]++;
        for (auto& kv : seen){ n_reg[kv.first]++; if (kv.second > 1) n_revisit[kv.first]++; }
    }
    ankerl::unordered_dense::set<int> detectable;
    for (Unitig* u : utgs){
        if (u == NULL || u->get_component() != comp_num || !g_threads.placeable.count(u->get_id())) continue;
        double need = std::max((double)THREAD_MIN_READS, THREAD_DETECT_FRAC * u->get_median_freq());
        int reg = n_reg[u->get_id()], rev = n_revisit[u->get_id()];
        bool repeat = rev >= THREAD_MIN_READS && rev >= THREAD_REVISIT_FRAC * reg;
        if (reg >= need && !repeat) detectable.insert(u->get_id());
        if (getenv("DK_THREAD_DUMP")) std::cout << "TDETECT utg" << u->get_id() << " MF=" << u->get_median_freq() << " reads=" << reg << " revisits=" << rev << (detectable.count(u->get_id()) ? "" : " excluded") << std::endl;
    }
    g_threads.placeable = detectable;
    std::vector<std::vector<int>> kept;
    for (auto& p : g_threads.paths){
        std::vector<int> q;
        for (int x : p) if (detectable.count(x) && (q.empty() || q.back() != x)) q.push_back(x);
        if (q.size() >= 2){ kept.push_back(std::move(q)); n_multi++; }
    }
    g_threads.paths.swap(kept);
    g_threads.active = true;
    if (verbose) std::cout << "(" << comp_num << ") read threading: " << owner.size() << " placement windows, " << n_multi << " reads span >= 2 unitigs" << std::endl;
    if (getenv("DK_THREAD_DUMP")){
        for (auto& p : g_threads.paths){
            std::cout << "TPATH";
            for (int x : p) std::cout << " " << x;
            std::cout << std::endl;
        }
    }
}

static int thread_support(const std::vector<int>& prefix, const ankerl::unordered_dense::set<int>& nexts){
    int n {0};
    size_t len = prefix.size() + 1;
    for (const std::vector<int>& p : g_threads.paths){
        if (p.size() < len) continue;
        bool found {false};
        for (size_t s{0}; s + len <= p.size() && !found; s++){
            bool fw {nexts.count(p[s + len - 1]) > 0}, rv {nexts.count(p[s]) > 0};
            for (size_t i{0}; i < prefix.size() && (fw || rv); i++){
                if (p[s + i] != prefix[i]) fw = false;
                if (p[s + len - 1 - i] != prefix[i]) rv = false;
            }
            found = fw || rv;
        }
        if (found) n++;
    }
    return n;
}

static void thread_signature(Node* n, int depth, ankerl::unordered_dense::set<int>& sig){
    Unitig* u = n->get_unitig();
    if (g_threads.placeable.count(u->get_id())){ sig.insert(u->get_id()); return; }
    if (depth <= 0) return;
    Node* opp = u->opposite_end(n);
    for (const Edge& e : opp->view_neighbors()){
        if (e.neighbor->get_unitig() == u) continue;
        thread_signature(e.neighbor, depth - 1, sig);
    }
}

static bool thread_filter(const Walk& walk, std::vector<Node*>& hap_neighs, std::vector<Node*>& other_neighs){
    if (!g_threads.active || walk.utgs.empty()) return false;
    std::vector<Node*> live;
    for (std::vector<Node*>* v : {&hap_neighs, &other_neighs}){
        for (Node* n : *v){
            Unitig* u = n->get_unitig();
            if (u->get_times_walked() < u->get_walk_limit()) live.push_back(n);
        }
    }
    if (live.size() < 2) return false;
    std::vector<ankerl::unordered_dense::set<int>> sigs(live.size());
    for (size_t c{0}; c < live.size(); c++) thread_signature(live[c], THREAD_SIG_DEPTH, sigs[c]);
    std::vector<int> placed;
    for (size_t i{walk.utgs.size()}; i-- > 0 && placed.size() < THREAD_MAX_SUFFIX; ){
        int id = walk.utgs[i]->get_id();
        if (g_threads.placeable.count(id) && (placed.empty() || placed.back() != id)) placed.push_back(id);
    }
    std::reverse(placed.begin(), placed.end());
    if (placed.empty()) return false;
    for (size_t m{placed.size()}; m >= 1; m--){
        std::vector<int> prefix(placed.end() - (long)m, placed.end());
        std::vector<int> sup(live.size());
        int total {0}, best {0};
        for (size_t c{0}; c < live.size(); c++){
            ankerl::unordered_dense::set<int> nexts;
            for (int x : sigs[c]) if (x != prefix.back()) nexts.insert(x);
            sup[c] = nexts.empty() ? 0 : thread_support(prefix, nexts);
            total += sup[c];
            best = std::max(best, sup[c]);
        }
        if (total < THREAD_MIN_READS) continue;
        {
            bool overlap {false};
            for (size_t a{0}; a < live.size() && !overlap; a++)
                for (size_t b{a + 1}; b < live.size() && !overlap; b++)
                    for (int x : sigs[a]) if (sigs[b].count(x)) { overlap = true; break; }
            if (overlap) return false;
        }
        int keep_min = std::max(1, (int)std::ceil(THREAD_MIN_FRAC * best));
        ankerl::unordered_dense::set<Node*> allowed;
        for (size_t c{0}; c < live.size(); c++) if (sup[c] >= keep_min) allowed.insert(live[c]);
        if (allowed.size() == live.size()) return false;
        auto keep = [&](std::vector<Node*>& v){
            std::vector<Node*> out;
            for (Node* n : v) if (allowed.count(n)) out.push_back(n);
            v.swap(out);
        };
        keep(hap_neighs);
        keep(other_neighs);
        if (getenv("DK_THREAD_DUMP")){
            std::cout << "TCHOICE at utg" << walk.utgs.back()->get_id() << " suffix";
            for (int x : prefix) std::cout << " " << x;
            std::cout << " :";
            for (size_t c{0}; c < live.size(); c++) std::cout << " utg" << live[c]->get_unitig()->get_id() << "=" << sup[c] << (allowed.count(live[c]) ? "" : "(x)");
            std::cout << std::endl;
        }
        return true;
    }
    return false;
}

void phase::build_walk(Walk& walk, Node* start_node, int hap_num){
    int other_hap_num = (hap_num == HAP1) ? HAP2 : HAP1;
    Node* curr = start_node;
    while (curr != NULL){
        Unitig* add_utg = curr->get_unitig();
        int times_walked = add_utg->get_times_walked();
        if (times_walked >= add_utg->get_walk_limit()) break;
        add_utg->set_times_walked(times_walked + 1);
        walk.utgs.push_back(add_utg);
        if (add_utg->get_assignment() == hap_num || add_utg->get_assignment() == HOM){
            walk.num_hap_nodes += add_utg->num_nodes();
        }
        Node* opp_node = add_utg->opposite_end(curr);
        const EdgeList& opp_neighbors = opp_node->view_neighbors();
        std::vector<Node*> hap_neighs;
        std::vector<Node*> other_neighs;
        // find neighbor haplotype assignments
        for (size_t i{0}; i < opp_neighbors.size(); i++){
            Node* neigh = opp_neighbors[i].neighbor;
            Unitig* neigh_utg = neigh->get_unitig();
            if (neigh_utg->get_id() != opp_node->get_unitig()->get_id()){
                if (neigh_utg->get_assignment() == hap_num) hap_neighs.push_back(neigh);
                if(neigh_utg->get_assignment() == HOM) hap_neighs.push_back(neigh);
                if(neigh_utg->get_assignment() == other_hap_num) other_neighs.push_back(neigh);
            }
        }
        thread_filter(walk, hap_neighs, other_neighs);
        // prioritize paths that are this haplotype or homozygous, then other if needed
        std::vector<Node*>* choices = NULL;

        if (!hap_neighs.empty()){
            choices = &hap_neighs;
        }else if(!other_neighs.empty()){
            choices = &other_neighs;
        }
        if (choices == NULL) break;
        Node* next = (*choices)[0];
        if (choices->size() > 1){
            int best_score {-1};
            for (size_t i{0}; i < choices->size(); i++){
                int budget {LOOKAHEAD_BUDGET};
                int score = phase::lookahead_score((*choices)[i], hap_num, budget);
                if (score > best_score){ best_score = score; next = (*choices)[i]; }
            }
        }
        curr = next;
    }
}

int phase::lookahead_score(Node* to_add, int hap_num, int& budget){
    int other_hap_num = (hap_num == HAP1) ? HAP2 : HAP1;
    Unitig* add_utg = to_add->get_unitig();
    int times_walked = add_utg->get_times_walked();
    if (budget <= 0 || times_walked >= add_utg->get_walk_limit()) return 0;
    budget--;

    add_utg->set_times_walked(times_walked + 1);
    int hap_weight = (add_utg->get_assignment() == hap_num || add_utg->get_assignment() == HOM) ? add_utg->num_nodes() : 0;

    Node* opp_node = add_utg->opposite_end(to_add);
    const EdgeList& opp_neighbors = opp_node->view_neighbors();
    std::vector<Node*> hap_neighs;
    std::vector<Node*> other_neighs;
    // find neighbor haplotype assignments
    for (size_t i{0}; i < opp_neighbors.size(); i++){
        Node* neigh = opp_neighbors[i].neighbor;
        Unitig* neigh_utg = neigh->get_unitig();
        if (neigh_utg->get_id() != opp_node->get_unitig()->get_id()){
            if (neigh_utg->get_assignment() == hap_num) hap_neighs.push_back(neigh);
            if(neigh_utg->get_assignment() == HOM) hap_neighs.push_back(neigh);
            if(neigh_utg->get_assignment() == other_hap_num) other_neighs.push_back(neigh);
        }
    }
    std::vector<Node*>* choices = NULL;
    if (!hap_neighs.empty()) {
        choices = &hap_neighs;
    } else if (!other_neighs.empty()) {
        choices = &other_neighs;
    }

    int best_ext {0};
    if (choices != NULL){
        for (size_t i{0}; i < choices->size() && budget > 0; i++){
            int ext = phase::lookahead_score((*choices)[i], hap_num, budget);
            if (ext > best_ext) best_ext = ext;
        }
    }
    add_utg->set_times_walked(times_walked);
    return hap_weight + best_ext;
}



void phase::phase_component(int comp_num, std::vector<std::string>& phased_walks, Graph* prev_graph, std::vector<SketchedSeq*>& corrected_seqs, bool verbose){
    // this will be filled by next function call
    std::vector<Unitig*> utg_list;
    phase::initialize_unitigs_for_tagging(prev_graph, corrected_seqs, utg_list, comp_num);
    if (utg_list.size() == 1){
        std::string utg_seq;
        utg_list[0]->get_hpc_sequence(utg_seq, false, false);
        phased_walks.emplace_back(utg_seq);
    }else{
        bool has_bubble {false};
        for (size_t i{0}; i < utg_list.size(); i++){
            if (utg_list[i]->get_parallel_utg() != NULL){ has_bubble = true; break; }
        }
        if (!has_bubble){
            for (size_t i{0}; i < utg_list.size(); i++) utg_list[i]->set_assignment(HOM);
            std::vector<std::string> path_seqs;
            std::vector<Unitig*> walked;
            build_threads(prev_graph, corrected_seqs, comp_num, verbose);
            bool got = phase::get_haplotype_sequences(prev_graph, path_seqs, comp_num, &walked);
            g_threads.active = false; g_threads.paths.clear();
            if (got && !path_seqs.empty()){
                phased_walks.emplace_back(path_seqs[0]);
            }
            std::unordered_set<Unitig*> visited(walked.begin(), walked.end());
            size_t n_extra {0};
            for (size_t i{0}; i < utg_list.size(); i++){
                if (visited.count(utg_list[i]) == 0){
                    std::string utg_seq;
                    utg_list[i]->get_hpc_sequence(utg_seq, false, false);
                    phased_walks.emplace_back(utg_seq);
                    n_extra++;
                }
            }
            if (verbose) std::cout << "(" << comp_num << ") no bubble: emitted the homozygous path once plus " << n_extra << " unvisited unitigs" << std::endl;
            return;
        }
        int round {0};
        while(utg_list.size() > 0){
            if (!phase::seed_initial_haplotypes(utg_list)){
                // no more bubbles, remaining unitigs are homozygous
                for (size_t i{0}; i < utg_list.size(); i++){
                    utg_list[i]->set_assignment(HOM);
                }
                break;
            }
            int it {0};
            while(phase::reassign_partitions(prev_graph, corrected_seqs, utg_list) && it < 100){
                if (verbose && it % 10 == 0) std::cout << "(" << comp_num << ") Phasing round " << round << ", iteration " << it << std::endl;
                it++;
            }
            phase::remove_tagged_unitigs(utg_list);
            round++;
        }
        if (verbose) std::cout << "(" << comp_num << ") Phasing assignments complete" << std::endl;
        build_threads(prev_graph, corrected_seqs, comp_num, verbose);
        phase::get_haplotype_sequences(prev_graph, phased_walks, comp_num);
        g_threads.active = false; g_threads.paths.clear();
    }
}

void phase::assign_reads_by_syncmers(Graph* graph, std::vector<SketchedSeq*>& reads, const std::vector<uint64_t>& neutral, bool verbose){
    const uint8_t IN_HAP1 = 1, IN_HAP2 = 2, IN_HOM = 4;
    ankerl::unordered_dense::map<uint64_t, uint8_t> where;
    const std::vector<Unitig*> utgs = graph->view_unitigs();
    for (size_t i{0}; i < utgs.size(); i++){
        Unitig* utg = utgs[i];
        if (utg == NULL) continue;
        int a = utg->get_assignment();
        uint8_t bit = (a == HAP1) ? IN_HAP1 : (a == HAP2) ? IN_HAP2 : IN_HOM;
        for (int n{0}; n < utg->num_nodes(); n++){
            Node* node = utg->node_at((size_t)n);
            int nk = node->num_kmers();
            if (n == 0){
                for (int k{0}; k < nk; k++) where[node->ith_kmer(k)] |= bit;
            }else{
                where[node->ith_kmer(0)] |= bit;
                where[node->ith_kmer(nk - 1)] |= bit;
            }
        }
    }
    for (uint64_t k : neutral) where[k] |= IN_HOM;
    size_t n_set {0}, n_changed {0}, n_conflict {0};
    for (size_t i{0}; i < reads.size(); i++){
        SketchedSeq* r = reads[i];
        if (!r->is_read) continue;
        int num1 {0}, num2 {0};
        for (size_t j{0}; j < r->kmers.size(); j++){
            auto it = where.find(r->kmers[j]);
            if (it == where.end()) continue;
            if (it->second == IN_HAP1) num1++;
            else if (it->second == IN_HAP2) num2++;
        }
        int winner = 0;
        if (num1 >= 2 && num1 > 2 * num2) winner = HAP1;
        else if (num2 >= 2 && num2 > 2 * num1) winner = HAP2;
        if (winner == 0) continue;
        if (r->assignment == HAP1 || r->assignment == HAP2){
            if (r->assignment != winner) n_conflict++;
            continue;
        }
        r->assignment = winner;
        n_set++;
    }
    (void)n_changed;
    if (verbose) std::cout << "syncmer read assignment: " << n_set << " reads newly assigned, " << n_conflict << " disagreed with the phasing vote (kept)" << std::endl;
}
