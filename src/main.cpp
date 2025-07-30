#include "argument_parser.h"
#include "deknot.h"
#include "graph.h"
#include "hash.h"
#include "hpc.h"
#include "syncmer.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

void parse_fq(ArgumentParser& input, std::vector<std::string>& seq_list);

int main(int argc, char* argv[]){
    ArgumentParser input(argc, argv);
    std::cout << '\n';
    if (input.args["command"] == "--help"){
        print_help();
    }else if (input.args["command"] == "assemble"){
        if (!input.check_arguments()) return 1;
        std::vector<std::string> read_seqs;
        parse_fq(input, read_seqs);

        char** reads_arr = new char*[read_seqs.size()];
        for (size_t i{0}; i < read_seqs.size(); i++){
            reads_arr[i] = new char[read_seqs[i].size() + 1];
            std::strcpy(reads_arr[i], read_seqs[i].c_str());
        }

        dk_opt_t asm_opt;
        dk_opt_init(&asm_opt);
        input.set_options(asm_opt);
        asm_opt.out_dir = input.args["-o"].c_str();

        dk_graph_t* final_graph = dk_assemble_reads(&asm_opt, reads_arr, (int)read_seqs.size());
        if (final_graph == NULL) return 1;

        std::string out_gfa = input.args["-o"] + "/final_graph.gfa";
        dk_write_gfa(out_gfa.c_str(), final_graph);
    
        // free memory
        dk_graph_free(final_graph);
        for (size_t i{0}; i < read_seqs.size(); i++){
            delete[] reads_arr[i];
        }
        delete[] reads_arr;
    }
    return 0;
}


void parse_fq(ArgumentParser& input, std::vector<std::string>& seq_list){
    std::ifstream in_file(input.args["--input"]);
    std::string name, orig_seq, ignore;
    while (in_file >> name){
        in_file >> orig_seq >> ignore >> ignore;
        seq_list.push_back(orig_seq);
    }
}
