#include "hash.h"
#include "io.h"

#include <fstream>

void io::parse_fq(ArgumentParser& input, std::vector<std::string>& seq_list, std::vector<std::string>& names_list){
    std::ifstream in_file(input.args["--input"]);
    std::string name, orig_seq, ignore;
    while (in_file >> name){
        in_file >> orig_seq >> ignore >> ignore;
        seq_list.push_back(orig_seq);
        names_list.push_back(name.substr(1));
    }
}

void io::write_sequences(std::vector<SketchedSeq*>& seqs, std::string out_path){
    std::ofstream out_file(out_path);
    for (size_t i{0}; i < seqs.size(); i++){
        out_file << ">" << seqs[i]->name << '\n';
        std::string hpc_seq {""};
        hash::hash_to_kmer(seqs[i]->kmers[0], seqs[i]->flipped[0], hpc_seq);
        for (size_t j{1}; j < seqs[i]->kmers.size(); j++){
            hash::add_to_seq_from_hash(seqs[i]->kmers[j], seqs[i]->flipped[j], seqs[i]->unique_seq_len[j], hpc_seq);
        }
        out_file << hpc_seq << '\n';
    }
    out_file.close();
}

void io::write_walks(std::vector<std::string>& seqs, std::string out_path){
    std::ofstream out_file(out_path);
    for (size_t i{0}; i < seqs.size(); i++){
        out_file << ">unitig" << i << '\n';
        out_file << seqs[i] << '\n';
    }
}
