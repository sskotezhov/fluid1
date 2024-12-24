#include "simulator.hpp"
#define TYPES FIXED(32, 16),DOUBLE,FLOAT
#define SIZES S(36, 85)
#ifndef TYPES
#error "TYPES is not defined"
#endif
#ifdef SIZES

using Simulator = types::Simulator<types::TypesList<TYPES>, SIZES>;

#else

using Simulator = types::Simulator<types::TypesList<TYPES>>;

#endif
//./main --p-type="FIXED(32,16)" --v-type="FIXED(32,16)" --v-flow-type="FIXED(32,16)" blep.txt
int main(int argc, char **argv) {
    if (argc < 5)
    {
        throw std::invalid_argument("Not enought arguments");
    }
    std::unordered_map<std::string, std::string> args;
    std::string filename;
    for (int i = 0; i < argc; i++)
    {
        std::string tmp = argv[i];
        char *t;
        if (tmp.find("=") != std::string::npos)
        {
            t = strtok(argv[i], "=");
            std::string tmp = t;
            t = strtok(NULL, "=");
            args[tmp] =  std::string(t);
        }
        else if (tmp.find(".") != std::string::npos)
        {
            filename = argv[i];
        }
    }
    const Simulator simulator = Simulator::from_params(types::SimulationParams{
        .p_type_name      = args["--p-type"],
        .v_type_name      = args["--v-type"],
        .v_flow_type_name = args["--v-flow-type"],
    });
    std::ifstream in(filename); 
    std::vector<std::string> fl;
    size_t rows = 0, columns = 0;
    if (in.is_open())
    {
        fl.push_back(std::string());
        std::getline(in, fl.back());
        while (&fl.back() == &fl[0] || fl[0] != fl.back())
        {
            fl.push_back(std::string());
            std::getline(in, fl.back());
        }
        rows = fl.size();
        columns = fl.back().size();
    }
    else 
    {
        exit(-1);
    }
    types::data_from_file tmp(rows, columns, in);
    in.close();
    std::vector<std::vector<char>> input;
    for (auto i : fl)
    {
        input.push_back(std::vector<char>());
        for (auto j : i)
        {
            input.back().push_back(j);
        }
    }
    simulator.start_on_field(types::Context{
        .field = input,
        .data = tmp,
    });
}