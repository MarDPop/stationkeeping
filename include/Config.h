#pragma once

#include <vector>

struct NestedCharTable {
    static const char VALS[38] = {  'A','B','C','D','E','F','G','H','I','J','K','L',
                                    'M','N','O','P','Q','R','S','T','U','V','W','X',
                                    'Y','Z','0','1','2','3','4','5','6','7','8','9','_',''};

    std::string name = "";

    int idx = -1;

    char type = -1;

    std::array<NestedCharTable*,38> chars = {NULL}; 

    NestedCharTable();
    NestedCharTable(std::string NAME);
    ~NestedCharTable();
}

class Config {

    Config(){}
    Config(const Config &){}
    ~Config(){}

    static Config* _instance = 0x0;

    NestedCharTable * const index = new NestedCharTable();

    std::vector<std::string> strings;

    std::vector<double> doubles;

    std::vector<int> integers;

    std::vector<bool> flags;

    NestedCharTable* get_table(const std::string& name);

public: 

    static Config* get_instance();

    static void convert_name(std::string& name);

    std::string check_param(std::string name) const;

    static void reset();

    static void load(const std::string& location);

    template<typename T>
    static void add(const std::string& name, T val);

    static void remove(const std::string& name);

    template<typename T>
    static T get(const std::string& name);
};