//
// Created by Eric Mandolesi on 30/05/2020.
//

#include "objects.h"
#include "global.h"
using namespace mtobj;
int main (int argc, char *argv[]){
    std::ifstream is(argv[1]);
    boost::archive::text_iarchive ia(is);
    mtobj::model model;
    ia >> model;
    model.calc_params();

    std::cout << model(0.0001) <<std::endl;
    std::cout << model(0.0001).maxAbsImpedance() <<std::endl;
    return 0;
}