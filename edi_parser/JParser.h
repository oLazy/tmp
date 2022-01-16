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
    typedef map<string,vector<double>> Data_Block;


    class Parser {
    public:
        void parse(string const&);
        void printDataBlock();

    private:
        Information_Block informationBlock;
        void parse(vector<Token> &);
        Tokenizer tokenizer;
//
//        void data_block(vector<Token>::iterator &it);
//        string station_name{};
    };
}

#endif //MTPARSER_PARSER_H
