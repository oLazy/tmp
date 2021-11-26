//
// Created by eric on 16/04/2021.
//

#include "JParser.h"
#include <cmath>
#include <exception>
namespace MTparser {
    void Parser::parse(const string &fileContents) {
        vector<Token> tokens = tokenizer.parse(fileContents);
        parse(tokens);
    }
    void Parser::parse(vector<Token> &t){
        for(auto ti = t.begin(); ti != t.end();ti++){
            data_block(ti);
        }
    }
    void Parser::data_block(vector<Token>::iterator &it) {
        if(it->type == KEYWORD && next(it)->type==OPERATOR && next(it,2)->type != KEYWORD){ // info block
            auto keyword = it->text;
            it++;
            it++;
            auto value = stod(it->text);
        } else if (it->type == STRING){ //station name, data block
            station_name = it->text;
            it++;
            if (it->type == STRING){
                if (it->text[0] == 'Z'){//read file

                } else {
                    throw std::invalid_argument("not ready to read apparent resistivity and phases.");
                }
            } else {
                throw logic_error("MALFORMED J-FILE (No string after station name).");
            }
        } else {
            throw logic_error("MALFORMED J-FILE.");
        }
    }
/*
    void read_impedence_data(vector<Token>::iterator &it){
        if(it->type == STRING && (it->text[0] == 'Z' || it->text[0] == 'z')){
            auto block_name = it->text;
            it++;
            int n_freq = stoi(it->text);
            it++;
            for (auto i = 0; i<n_freq; i++){
                for (auto j=0; i<)
            }
        }
    }
    */
}