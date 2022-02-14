//
// Created by eric on 11/02/2022.
//

#ifndef MT1DANISMODELPARAMS_SAMPLER_INIT_H
#define MT1DANISMODELPARAMS_SAMPLER_INIT_H
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "objects.h"

void init_sampler(boost::program_options::variables_map const& input){
    // load data
    std::string dfile;
    if(input.count("data-filename")){
        dfile = input["data-filename"].as<std::string>();
    }else{
        dfile = input["base-filename"].as<std::string>() + "_rep.bin";
        std::cerr << "data-filename not provided. Trying " + dfile +"\n";
    }
    if( !boost::filesystem::exists(dfile))throw std::runtime_error(dfile + " file not found.");
    auto readData = mtobj::io::loadDataset(dfile);

//    if(readData.cov_type!=mtobj::io::cov_code::real)throw std::runtime_error("covariance type not implemented yet");
// implementing only the Cov1 case
    mtobj::Cov1 cov = (readData.cov_type==mtobj::io::cov_code::real)? mtobj::initFrom(readData.c0) : readData.c1;
    auto d=readData.d;

    // load model
    std::string mfile;
    if(input.count("m0-filename")){
        mfile = input["m0-filename"].as<std::string>();
    }else{
        mfile = input["base-filename"].as<std::string>() + ".bin";
        std::cerr << "m0-filename not provided. Trying " + mfile +"\n";
    }
    mtobj::model m;
    if( !boost::filesystem::exists(mfile)) {
        std::cerr << mfile + " file not found.\n"
                             "Setting up uniform isotropic half-space with resistivity 10 Ohm per meter.\n";
        m.nodes.emplace_back(0,-1);
    } else {
        m = mtobj::io::load(mfile);
    }
    m.calc_params();
    if(!m.isInPrior())throw std::runtime_error("initial model not in prior.\n");
    m.setLogL(mtobj::logL(m, d, cov));
}
#endif //MT1DANISMODELPARAMS_SAMPLER_INIT_H
