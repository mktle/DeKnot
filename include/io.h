#ifndef IO_H
#define IO_H

#include "argument_parser.h"

namespace io{
    void parse_fq(ArgumentParser& input, std::vector<std::string>& seq_list, std::vector<std::string>& names_list);
    void write_sequences();
}

#endif
