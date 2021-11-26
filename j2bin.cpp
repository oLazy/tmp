//
// Created by eric on 06/06/2021.
//
#include <iostream>
#include <fstream>
#include <strstream>
//#include "edi_parser/Parser.h"
#include "edi_parser/JTokenizer.h"
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
    Tokenizer tokenizer;

    auto tokens = tokenizer.parse(fileContents);
    for (auto t: tokens){
        std::cout << t.type << "; " << t.text << "\n";
    }

//    Parser p;
//    p.parse(fileContents);
//    auto skip_string = p.get_option_list_for(">HEAD")["EMPTY"];
//    auto freq = MTparser::dataset2double(p.get_data_set_for(">FREQ"), skip_string);
//    // OFF-DIAG
//    auto zxyr = MTparser::dataset2double(p.get_data_set_for(">ZXYR"), skip_string);
//    auto zyxr = MTparser::dataset2double(p.get_data_set_for(">ZYXR"), skip_string);
//    auto zxyi = MTparser::dataset2double(p.get_data_set_for(">ZXYI"), skip_string);
//    auto zyxi = MTparser::dataset2double(p.get_data_set_for(">ZYXI"), skip_string);
//    // MAIN-DIAG
//    auto zxxr = MTparser::dataset2double(p.get_data_set_for(">ZXXR"), skip_string);
//    auto zyyr = MTparser::dataset2double(p.get_data_set_for(">ZYYR"), skip_string);
//    auto zxxi = MTparser::dataset2double(p.get_data_set_for(">ZXXI"), skip_string);
//    auto zyyi = MTparser::dataset2double(p.get_data_set_for(">ZYYI"), skip_string);
//
//    auto zxyv = MTparser::dataset2double(p.get_data_set_for(">ZXY.VAR"), skip_string);
//    auto zyxv = MTparser::dataset2double(p.get_data_set_for(">ZYX.VAR"), skip_string);
//    auto zxxv = MTparser::dataset2double(p.get_data_set_for(">ZXX.VAR"), skip_string);
//    auto zyyv = MTparser::dataset2double(p.get_data_set_for(">ZYY.VAR"), skip_string);
//
//    int i = 0;
//    double conversion_factor{1};
//    if(change_units){
//        conversion_factor = 4*M_PI*0.0001;
//    }
//    for (auto f:freq){
//        bool skip{true};
//        auto T = pow(f,-1);
//        std::array<double,8> thisT {zxxr[i],
//                                    zxxi[i],
//                                    zxyr[i],
//                                    zxyi[i],
//                                    zyxr[i],
//                                    zyxi[i],
//                                    zyyr[i],
//                                    zyyi[i]};
//        for (auto e:thisT){
//            if (std::isnan(e)) {
//                skip=true;
//                break;
//            }else{
//                skip=false;
//            }
//        }
//        if(!skip) {
//            MTTensor z{{conversion_factor*zxxr[i], conversion_factor*zxxi[i]},
//                       {conversion_factor*zxyr[i], conversion_factor*zxyi[i]},
//                       {conversion_factor*zyxr[i], conversion_factor*zyxi[i]},
//                       {conversion_factor*zyyr[i], conversion_factor*zyyi[i]}};
//            dataset[T] = z;
//            MTTensor vz{{conversion_factor*conversion_factor*zxxv[i], 0},
//                       {conversion_factor*conversion_factor*zxyv[i], 0},
//                       {conversion_factor*conversion_factor*zyxv[i], 0},
//                       {conversion_factor*conversion_factor*zyyv[i], 0}};
//            covariance[T] = vz;
//
//            std::cout << T << ": " << z << "\n" << vz <<"\n";
//        }
//        i++;
//    }
//
//    // save tensors and periods
//    std::ofstream os(out_file_name);
//    boost::archive::text_oarchive oa(os);
//    oa << dataset << covariance;
////    oa << cov_0;
//    os.close();
//
//    std::cout << "skip string:" << skip_string << "\n";
//
////    for (i = 0; i< zxxv.size();i++) {
////        std::cout << "T: " << std::pow(freq[i],-1) << "\n"
////                  << "std(xx): " << std::sqrt(zxxv[i]) << "\n"
////                  << "std(xy): " << std::sqrt(zxyv[i]) << "\n"
////                  << "std(yx): " << std::sqrt(zyxv[i]) << "\n"
////                  << "std(yy): " << std::sqrt(zyyv[i]) << "\n";
////    }
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