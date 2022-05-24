#include "../include/Config.h"

#include <fstream>

NestedCharTable::NestedCharTable(){
}

NestedCharTable::NestedCharTable(std::string NAME) : name(NAME){
}

NestedCharTable::~NestedCharTable(){
    for(int i = 0; i < 38;i++){
        delete this->chars[i];
    }
}

Config* Config::get_instance() {
    if(_instance = 0x0){
        _instance = new Config();
    }

    return _instance;
}

void Config::convert_name(std::string& name){
    const int nChar = name.size();

    for(int i = 0; i < nChar; i++){
        char& c = name[i];
        if( c > 96 && c < 123 ) {
            c -= 32;
        }

        if( c > 64 && c < 91 ) {
            c -= 65;
        }

        if( c > 47 && c < 58 ) {
            c -= 23;
        }

        if( c == 95 ) {
            c = 36;
        }

        if( c > 37 ) {
            throw "Invalid Parameter Name.";
        }
    }
}

void Config::reset(){
    Config* config = Config::get_instance();

    delete config->index;

    config->index = new NestedCharTable();

    config->doubles.clear();
    config->strings.clear();
    config->integers.clear();
    config->flags.clear();
}

void Config::load(const std::string& location) {
    std::vector<std::string> lines;
    std::string line;

    std::ifstream input_file(location);
    if (!input_file.is_open()) {
        throw "Couldn't open: " + location;
    }

    while (getline(input_file, line)){
        lines.push_back(line);
    }

    input_file.close();
}

std::string Config::check_param(std::string name) const {
    const int nChar = name.size();

    if ( nChar < 2 ) {
        throw "invalid parameter length.";
    }

    Config::convert_name(name);

    NestedCharTable* table = this->index;
    bool found = true;
    for(int i = 0; i < nChar; i++){
        if( table->chars[name[i]] == NULL ) {
            found = false;
        } else {
            table = table->chars[name[i]];
        }
    }
    
    if(found) {
        if (table->idx >= 0){
            throw "parameter already exists.";
        }
    }
    return name;
}

NestedCharTable* Config::get_table(const std::string& name) {
    NestedCharTable* table = this->index;
    const int nChar = name.size();
    std::string nestName = "";
    for(int i = 0; i < nChar; i++){
        nestName += name[i];
        if( table->chars[name[i]] == NULL ) {
            table->chars[name[i]] = new NestedCharTable(nestName);
        } 
        table = table->chars[name[i]];
    }
    return table;
}

template<>
void Config::add(const std::string& name, double val) {

    Config* config = Config::get_instance();

    NestedCharTable* table = config->get_table(config->check_param(name));

    table->idx = config->doubles.size();

    table->type = 1;

    config->doubles.push_back(val);
}

void Config::remove(const std::string& name){
    Config* config = Config::get_instance();

    NestedCharTable* table = config->get_table(name);

    switch(table->type) {
        case 0:
            config->strings.erase(config->strings.begin() + table->idx);
            break;
        case 1:
            config->doubles.erase(config->doubles.begin() + table->idx);
            break;
        case 2:
            config->integers.erase(config->integers.begin() + table->idx);
            break;
        case 3:
            config->flags.erase(config->flags.begin() + table->idx);
            break;
    }

    table->idx = -1;

    table->type = -1;
}

template<typename T>
T Config::get(const std::string& name){

}