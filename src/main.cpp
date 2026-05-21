#include "argument_parser.h"
#include "deknot.h"
#include "graph.h"
#include "hash.h"
#include "hpc.h"
#include "io.h"
#include "syncmer.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[]){
    ArgumentParser input(argc, argv);
    std::cout << '\n';
    if (input.args["command"] == "--help"){
        print_help();
    }else if (input.args["command"] == "assemble"){
        if (!input.check_arguments()) return 1;
        std::vector<std::string> read_seqs;
        std::vector<std::string> name_seqs;
        io::parse_fq(input, read_seqs, name_seqs);
        assert(read_seqs.size() == name_seqs.size());
        
        char** reads_arr = new char*[read_seqs.size()];
        char** names_arr = new char*[read_seqs.size()];
        for (size_t i{0}; i < read_seqs.size(); i++){
            reads_arr[i] = new char[read_seqs[i].size() + 1];
            std::strcpy(reads_arr[i], read_seqs[i].c_str());
            names_arr[i] = new char[name_seqs[i].size() + 1];
            std::strcpy(names_arr[i], name_seqs[i].c_str());
        }

        dk_opt_t asm_opt;
        dk_opt_init(&asm_opt);
        input.set_options(asm_opt);
        asm_opt.out_dir = input.args["-o"].c_str();

        dk_graph_t* final_graph = dk_assemble_reads(&asm_opt, reads_arr, names_arr, (int)read_seqs.size());
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
