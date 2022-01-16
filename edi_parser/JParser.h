//
// Created by eric on 16/04/2021.
//

#ifndef MTPARSER_PARSER_H
#define MTPARSER_PARSER_H
#include "JTokenizer.h"
#include <map>
#include <tuple>
#include <string>
namespace MTparser {
    using namespace std;
    typedef map<string,double> Information_Block;
    typedef array<vector<double>,5> Data_Table;
    typedef map<string,Data_Table> Data_Block; // map from string to a table of 5 columns
    enum dataMap{begin = 0, period = 0, real = 1, imag = 2, error = 3, weight = 4, end = 5};

    class Parser {
    public:
        void parse(string const&);
        void printInfoBlock();
        Data_Table getDataFor(std::string const &);
        unsigned getNfreq() const;


    private:
        Data_Block dataBlock;
        Information_Block informationBlock;
        void parse(vector<Token> &);
        Tokenizer tokenizer;
        unsigned nfreq{0};
//
//        void data_block(vector<Token>::iterator &it);
//        string station_name{};
    };
}

#endif //MTPARSER_PARSER_H
