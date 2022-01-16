//
// Created by eric on 06/06/2021.
//
#include <iostream>
#include <fstream>
#include <strstream>
//#include "edi_parser/JTokenizer.h"
#include "edi_parser/JParser.h"
#include "MTTensor.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <cmath>

std::string getFileContents(std::ifstream&);
boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm);
typedef std::map<double, MTTensor> cov1; // diagonal covariance matrix which allows different errors for different Z components


using namespace MTparser;
int main(int argn, char* argv[]){
    boost::program_options::variables_map vm;
    auto desc = parse_cmdline(argn, argv, vm);
    if(vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }
    std::map<double, MTTensor> dataset; // period, tensor
    std::map<double, MTTensor> covariance; // period, tensor, only real part of tensor is used
    auto fileName = vm["in-path"].as<std::string>();
    bool change_units = vm["field-units"].as<bool>();
    auto pathObject = boost::filesystem::path(fileName);
    auto outPathObject = boost::filesystem::path(vm["out-path"].as<std::string>());

    auto fname = pathObject.filename().string();

    auto out_f = boost::filesystem::change_extension(pathObject.filename(),".bin");
    auto dest = boost::filesystem::change_extension(pathObject.filename(),".bin");
    auto out_file_name = outPathObject.relative_path().append(dest.string()).string();;
    std::cout << "output path: " << outPathObject.string() << std::endl;
    std::cout << "handling file: " << fname << "\n";
    std::cout << "output file: " << out_file_name << std::endl;

    std::ifstream is;
    is.open(fileName);
    if(!is.is_open())throw std::runtime_error("Cannot open file " + fileName+ "\n");
    std::string fileContents = getFileContents(is);

    is.close();
//    Tokenizer tokenizer;

//    auto tokens = tokenizer.parse(fileContents);
    MTparser::Parser p;
    p.parse(fileContents);
    p.printInfoBlock();

    // build an impedance tensor from the parsed file
    Data_Table zxxData = p.getDataFor("ZXX");
    Data_Table zxyData = p.getDataFor("ZXY");
    Data_Table zyxData = p.getDataFor("ZYX");
    Data_Table zyyData = p.getDataFor("ZYY");
    for (int i=0; i<p.getNfreq(); i++) {
        double freq;
        if (
                zxxData[MTparser::dataMap::period][i] == zxyData[MTparser::dataMap::period][i] &&
                zxyData[MTparser::dataMap::period][i] == zyxData[MTparser::dataMap::period][i] &&
                zyxData[MTparser::dataMap::period][i] == zyyData[MTparser::dataMap::period][i]
                ) {
            freq = zxxData[MTparser::dataMap::period][i];
        } else {
            throw std::runtime_error("freqs difference problem.");
        }
        MTTensor z{{zxxData[MTparser::dataMap::real][i], zxxData[MTparser::dataMap::imag][i]}, //xx
                   {zxyData[MTparser::dataMap::real][i], zxyData[MTparser::dataMap::imag][i]}, //xy
                   {zyxData[MTparser::dataMap::real][i], zyxData[MTparser::dataMap::imag][i]}, //yx
                   {zyyData[MTparser::dataMap::real][i], zyyData[MTparser::dataMap::imag][i]}};

        MTTensor zv{
            {pow(zxxData[MTparser::dataMap::error][i],2), 0}, //xx
            {pow(zxyData[MTparser::dataMap::error][i],2), 0}, //xy
            {pow(zyxData[MTparser::dataMap::error][i],2), 0}, //yx
            {pow(zyyData[MTparser::dataMap::error][i],2), 0}};

        dataset[freq] = z;
        covariance[freq] = zv;
    }

    // save tensors and periods
    std::ofstream os(out_file_name);
    boost::archive::text_oarchive oa(os);
    oa << dataset << covariance;
    os.close();
    return 0;
}

std::string getFileContents(std::ifstream& input){
    std::ostrstream sstr;
    sstr << input.rdbuf();
    return sstr.str();
}

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "display this message and exit.")
            ("in-path,i", po::value<std::string>(), "PATH to the input .edi file.")
            ("out-path,o", po::value<std::string>()->default_value("./"), "Where I save .bin file.")
            ("field-units,u", po::value<bool>()->default_value(false), "Are tensor element in the edi file provided in field units?");
    po::store(po::parse_command_line(argc, argv, generic), p_vm);
    po::notify(p_vm);
    return generic;
}