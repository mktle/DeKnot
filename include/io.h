#ifndef IO_H
#define IO_H

#include "argument_parser.h"
#include "syncmer.h"

#include <string>
#include <vector>

namespace io{
    void parse_fq(ArgumentParser& input, std::vector<std::string>& seq_list, std::vector<std::string>& names_list);
    void write_sequences(std::vector<SketchedSeq*>& seqs, std::string out_path);
    void write_walks(std::vector<std::string>& seqs, std::string out_path);
}

#endif
