#include "io.h"

#include <fstream>

void io::parse_fq(ArgumentParser& input, std::vector<std::string>& seq_list, std::vector<std::string>& names_list){
    std::ifstream in_file(input.args["--input"]);
    std::string name, orig_seq, ignore;
    while (in_file >> name){
        in_file >> orig_seq >> ignore >> ignore;
        seq_list.push_back(orig_seq);
        names_list.push_back(name);
    }
}

void io::write_sequences(){


}
