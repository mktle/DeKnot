#include "phase.h"

#include <iostream>

void phase::initialize_unitigs_for_tagging(Graph* graph, std::vector<SketchedSeq*>& reads, std::vector<Unitig*>& sorted_unitigs){
    /************************
     * Map reads to unitigs *
     ************************/
    for (size_t i{0}; i < reads.size(); i++){
        if (reads[i]->is_read){
            ankerl::unordered_dense::set<Unitig*> alned_to;
            Unitig* track_utg {NULL};
            
            size_t check_thresh = reads[i]->kmers.size() - (size_t)graph->get_k_prime() + 1;
            for (size_t j{0}; j < check_thresh; j++){
                Span<uint64_t> curr_kmers(reads[i]->kmers.data() + j, (size_t)graph->get_k_prime());
                char curr_orient = syncmer::canonical_syncmer(curr_kmers);
                KmerView curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION};
                Node* curr_node = graph->get_node_from_lookup_table(curr_canon);
                if (curr_node != NULL && track_utg != curr_node->get_unitig()){
                    alned_to.insert(curr_node->get_unitig());
                    track_utg = curr_node->get_unitig();
                }
            }
            // now add read record to all unitigs
            for (const auto& aln_utg : alned_to) {
                aln_utg->add_read(reads[i]);
            }
        }
    }

    /*********************************************************
     * Tag unitigs as being part of bubbles or heterozygous  *
     *********************************************************/
    const std::vector<Unitig*> unitig_list = graph->view_unitigs();
    for (size_t i{0}; i < unitig_list.size(); i++){
        Unitig* utg = unitig_list[i];
        if (utg != NULL){
            utg->update_median_freq();
            sorted_unitigs.push_back(utg);
            // set bubble information; should not change now that the graph is fixed
            if (!utg->get_checked_bubble()){
                assert(utg->get_front() != utg->get_back());
                utg->determine_parallel();
                if (utg->get_in_bubble()){
                    utg->set_assignment(UNASSIGNED);
                }else{
                    utg->set_assignment(HET);
                }
            }
        }
    }
    std::sort(sorted_unitigs.begin(), sorted_unitigs.end(), compare_unitigs);
}

void phase::choose_initial_bubble(std::vector<Unitig*>& sorted_unitigs){
    /***************************************
     * Choose bubble to start phasing from *
     ***************************************/
    Unitig* initial_utg = sorted_unitigs[0];
    size_t i{1};
    while(!initial_utg->get_in_bubble()){
        assert(i < sorted_unitigs.size());
        initial_utg = sorted_unitigs[i];
        i++;
    }
    initial_utg->set_assignment(HAP1);
    initial_utg->assign_reads();
    initial_utg->get_parallel_utg()->set_assignment(HAP2);
    initial_utg->get_parallel_utg()->assign_reads();
}

// return true if any assignments changed
bool phase::reassign_partitions(Graph* graph, std::vector<SketchedSeq*>& reads){
    bool changed {false};
    /**************************
     * Reassign unitig labels *
     **************************/
    const std::vector<Unitig*> unitig_list = graph->view_unitigs();
    for (size_t i{0}; i < unitig_list.size(); i++){
        Unitig* utg = unitig_list[i];
        if (utg != NULL && utg->get_in_bubble()){
            bool flipped = utg->tally_assignment();
            if (flipped) changed = true;
        }
    }

    /************************
     * Reassign read labels *
     ************************/
    for (size_t i{0}; i < reads.size(); i++){
        if (reads[i]->is_read){
            int num1 {0};
            int num2 {0};

            size_t check_thresh = reads[i]->kmers.size() - (size_t)graph->get_k_prime() + 1;
            for (size_t j{0}; j < check_thresh; j++){
                Span<uint64_t> curr_kmers(reads[i]->kmers.data() + j, (size_t)graph->get_k_prime());
                char curr_orient = syncmer::canonical_syncmer(curr_kmers);
                KmerView curr_canon = {curr_kmers, curr_orient == REVERSE_ORIENTATION};
                Node* curr_node = graph->get_node_from_lookup_table(curr_canon);
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
