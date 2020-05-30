//
// Created by Eric Mandolesi on 30/05/2020.
//
#include "objects.h"
#include "global.h"
#include <sstream>
#include <string>
using namespace mtobj;
int main(int argn, char* argv[]){
    if (argn > 2){
        std::cout << "Error! this program requires a single argument.\n";
        return 1;
    }
    std::string out_file_name{"model.dat"}; // TODO make this custom
    std::ifstream is(argv[1]);
    std::string line;
    int ln = 0;
    mtobj::model m;
    while (std::getline(is,line))
    {
        std::istringstream iss(line);
        double z, h, l, b;
        if((iss >> z >> h)){
            //process
            m.nodes.push_back(mtobj::node(z, log10(h)));
        }
        else if((iss >> z >> h >> l >> b)){
            //process
            m.nodes.push_back(mtobj::node(z, log10(0.5*(h+l)), log10(h/l), b));
        }
        else{
            std::cerr << "error in input file, line " << ln << ". Please check.\n";
            return -ln;
        }

    }
    // check node 0 has depth 0
    if(m.nodes[0].params[mtobj::paramType::depth].getValue() == 0){
        std::ofstream os(out_file_name);
        boost::archive::text_oarchive oa(os);
        oa << m ;
    }else{
        return 2;
    }
    return 0;
}