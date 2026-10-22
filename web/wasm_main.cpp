// Entry point for the WebAssembly build of DeKnot (see web/README.md).
//
// The browser page writes the uploaded FASTQ into Emscripten's in-memory
// filesystem, calls dk_wasm_assemble(), and reads the output files back.
// This deliberately bypasses ArgumentParser (which shells out to mkdir) and
// runs with out_dir == NULL so no per-iteration checkpoint files are written.

#include "deknot.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

// Same parsing rules as io::parse_fq: whitespace-delimited name/seq/+/qual records,
// read name is the first token after '@'.
static void parse_fastq(const char* path, std::vector<std::string>& seqs, std::vector<std::string>& names){
    std::ifstream in_file(path);
    std::string name, seq, ignore;
    while (in_file >> name){
        in_file >> seq >> ignore >> ignore;
        seqs.push_back(seq);
        names.push_back(name.substr(1));
    }
}

extern "C" {

/**
 * Assembles the reads in fq_path and writes:
 *   walks_path   FASTA of expanded haplotype walks (same content as final_walks.fa)
 *   assign_path  TSV of read name -> haplotype (1/2, 0 = unassigned)
 *
 * @return number of walks written, -1 if the input had no reads, -2 if assembly failed
 */
EMSCRIPTEN_KEEPALIVE
int dk_wasm_assemble(const char* fq_path, const char* walks_path, const char* assign_path,
                     int s, int k_prime, int freq_thresh, int max_iterations, int check_strand, int verbose){
    std::vector<std::string> seqs;
    std::vector<std::string> names;
    parse_fastq(fq_path, seqs, names);
    if (seqs.empty()) return -1;

    std::vector<char*> reads_arr(seqs.size());
    std::vector<char*> names_arr(seqs.size());
    for (size_t i{0}; i < seqs.size(); i++){
        reads_arr[i] = const_cast<char*>(seqs[i].c_str());
        names_arr[i] = const_cast<char*>(names[i].c_str());
    }

    dk_opt_t opt;
    dk_opt_init(&opt);
    opt.s = s;
    opt.init_k_prime = k_prime;
    opt.compact_thresh = freq_thresh;
    opt.max_iterations = max_iterations;
    opt.check_strand = check_strand;
    opt.verbose = verbose;
    opt.out_dir = NULL;

    dk_walks_t* walks = dk_assemble_walks(&opt, reads_arr.data(), names_arr.data(), (int)seqs.size());
    if (walks == NULL) return -2;

    {
        std::ofstream out_file(walks_path);
        for (int i{0}; i < walks->n_walks; i++){
            out_file << ">unitig" << i << '\n' << walks->walks[i] << '\n';
        }
    }
    {
        std::ofstream out_file(assign_path);
        out_file << "read\thaplotype\n";
        for (size_t i{0}; i < names.size(); i++){
            out_file << names[i] << '\t' << walks->read_assignments[i] << '\n';
        }
    }
    std::cout.flush();
    int n_walks = walks->n_walks;
    dk_walks_free(walks);
    return n_walks;
}

}
