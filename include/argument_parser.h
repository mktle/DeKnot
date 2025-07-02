#ifndef ARGUMENT_PARSER_H
#define ARGUMENT_PARSER_H
#include "deknot.h"

#include <list>
#include <map>
#include <string>

class ArgumentParser{
    public:
        ArgumentParser(int &argc, char**argv);
        bool check_required_flags(std::list<std::string>& r_flags);
        bool check_file(std::string opt, std::string ext);
        bool check_arguments();
        void set_options(dk_opt_t& asm_opt);

        std::map<std::string, std::string> args;
};

void print_help();

#endif
