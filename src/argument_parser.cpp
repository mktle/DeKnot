#include "argument_parser.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

// typical call should have at least 3 argv members: "assembler", input file, and output directory
const int BASE_NUM_ARGS = 3;

ArgumentParser::ArgumentParser(int &argc, char** argv){
    if(argc == 1){
        this->args.insert({"command", "--help"});
    }else if(argc==2 && !std::strcmp(*(argv + 1), "--help")){
        this->args.insert({"command", "--help"});
    }else if (argc < BASE_NUM_ARGS){
        throw std::invalid_argument("Missing at least one required argument (input file output directory)");
    }else{
        // add the rest of the command-line options to map of args
        this->args.insert({"command", "assemble"});
        this->args.insert({"--input", *(argv + 1)});
        this->args.insert({"-o", *(argv + 2)});

        int i{BASE_NUM_ARGS};
        int prev_opt_index{-1};
        while (*(argv + i) != 0){
            // check if this is an option by looking for '-' prefix
            if (*(argv + i)[0] == '-'){
                // if this is not the first option, add previous option to map
                if (prev_opt_index != -1){
                    // check whether they specified arguments after flag
                    if (prev_opt_index == i - 1){
                        this->args.insert({*(argv + prev_opt_index), "true"});
                    }else{
                        // if valid, add options to argument map
                        std::string opts{*(argv + prev_opt_index + 1)};
                        for (int k{BASE_NUM_ARGS}; prev_opt_index + k < i; k++){
                            opts += ",";
                            opts += *(argv + prev_opt_index + k);
                        }
                        this->args.insert({*(argv + prev_opt_index), opts});
                    }
                }
                prev_opt_index = i;
            }
            i++;
        }
        if (argc > BASE_NUM_ARGS){
            // push final flag
            // check whether they specified arguments after flag
            if (prev_opt_index == i - 1){
                this->args.insert({*(argv + prev_opt_index), "true"});
            }else{
                std::string opts{*(argv + prev_opt_index + 1)};
                for (int k{2}; prev_opt_index + k < i; k++){
                    opts += ",";
                    opts += *(argv + prev_opt_index + k);
                }
                this->args.insert({*(argv + prev_opt_index), opts});
            }
        }
    }
}

bool ArgumentParser::check_required_flags(std::list<std::string>& r_flags){
    for (auto it = r_flags.begin(); it != r_flags.end(); it++){
        if (this->args.count(*it) == 0){
            std::cout << "[ArgumentParser::check_required_flags][ERROR] missing required flag " << *it << '\n';
            return false;
        }
    }
    return true;
}

bool ArgumentParser::check_file(std::string opt, std::string ext){
    struct stat buffer;
    size_t last_delim {this->args[opt].find_last_of('.')};
    if (last_delim == std::string::npos || this->args[opt].substr(last_delim) != ext){
        return false;
    }

    if (stat(this->args[opt].c_str(), &buffer) != 0){
        return false;
    }

    return true;
}

bool ArgumentParser::check_arguments(){
    if (!this->check_file("--input", ".fq") && !this->check_file("--input", ".fastq")){
        std::cout << "[ArgumentParser::check_file][ERROR] invalid input file specified\n";
        return false;
    }

    std::string cmd {"mkdir -p " + this->args["-o"]};
    if(system(cmd.c_str()) != 0){
        std::cout << "[ArgumentParser::check_arguments][ERROR] unable to make output directory " << this->args["-o"] << "\n";
        return false;
    }

    if (access(this->args["-o"].c_str(), W_OK)){
        std::cout << "[ArgumentParser::check_arguments][ERROR] output directory " << this->args["-o"] << " is not writable\n";
        return false;
    }
    return true;
}

void ArgumentParser::set_options(dk_opt_t& asm_opt){
    if (this->args.count("--k-prime") != 0){
        asm_opt.init_k_prime = std::stoi(this->args["--k-prime"]);
    }

    if (this->args.count("-s") != 0){
        asm_opt.s = std::stoi(this->args["-s"]);
    }

    if (this->args.count("--freq-thresh") != 0){
        asm_opt.compact_thresh = std::stoi(this->args["--freq-thresh"]);
    }

    if (this->args.count("--max-iterations") != 0){
        asm_opt.max_iterations = std::stoi(this->args["--max-iterations"]);
    }

    if (this->args.count("-v") != 0){
        asm_opt.verbose = true;
    }
}

void print_help(){
    std::cout << "Usage: [input.fq] [output directory] <flags>" << std::endl;
    std::cout << "Options:\n";
    std::cout << "  General:\n";
    std::cout << "    -v                verbose mode with logging statements\n";
    std::cout << "  Sketching:\n";
    std::cout << "    -s                s value when selecting syncmers\n";
    std::cout << "  Graph growth:\n";
    std::cout << "    --k-prime         initial number of syncmers per k-syncmer\n";
    std::cout << "    --freq-thresh     unitigs under this frequency will be dropped during graph compaction\n";
    std::cout << "    --max-iterations  maximum number of growth iterations before forcing DeKnot to terminate\n";
}

