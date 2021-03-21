//
// Created by eric on 20/03/2021.
//
#include "objects.h"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/program_options.hpp>
#include "global.h"

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm);
boost::program_options::options_description parse_config(boost::program_options::variables_map& p_vm);

int main(int argn, char* argv[]) {
    // parse command line and config file
    auto desc = parse_cmdline(argn, argv, vm);
    auto config_desc = parse_config(vm);
    if(vm.count("help")){
        std::cout << desc <<std::endl;
        std::cout << config_desc << std::endl;
        return 0;
    }
    bool done{false};
// load inverted dataset
    std::string base_filename{vm["base-filename"].as<std::string>().c_str()};
    mtobj::Cov0 cov;
    mtobj::Dataset d;
    std::ifstream is(base_filename+"_rep.dat");
    boost::archive::text_iarchive ia(is);
    ia >> d >> cov;
    is.close();
//    for (auto D : d){
//        std::cout << "T= " << D.first << ": Z= " << D.second << std::endl;
//    }
// create periods vector
    std::vector<double> periods;
    for (auto D : d){
        periods.push_back(D.first);
    }

    // create histograms
    auto ax_h = boost::histogram::axis::regular<>(periods.size()*2,
                                                  log10(periods[0]*0.9),
                                                  log10(periods[periods.size()-1]*1.1), "T");
    auto ax_v = boost::histogram::axis::regular<>(2001,
                                                  -0.05,
                                                  0.15, "component/absolute");
    auto zhist = boost::histogram::make_histogram(ax_h, ax_v);

    //load the first 1000 samples
    std::ifstream sample_output(base_filename+"_sample_out.bin");
    boost::archive::binary_iarchive sample_ia(sample_output);
    mtobj::model m;
    boost::circular_buffer<mtobj::model> s{1000};
//    for (auto i=0; i<1000; i++) {
    int iBuffer{0};
    int n_m_processed{0};
    while (!done){
        try {
            sample_ia >> m;
        }
        catch (boost::archive::archive_exception &e){
            std::cout << "catch exception " << e.what() <<"\n";
            done=true;
            break;
        }
        m.calc_params();
        s.push_back(m);
        iBuffer++;
        if(iBuffer==s.capacity()){
            iBuffer=0;
            for (auto mb=s.begin();mb!=s.end();mb++){
                n_m_processed++;
                for (auto T : periods) {
                    auto this_model = *mb;
                    this_model.calc_params();
                    auto z = this_model(T);
//                        zhist(log10(T),std::real(z.xy)/std::abs(z.xy));
//                    zhist(log10(T),std::real(z.xy)/std::abs(d[T].xy));
                    zhist(log10(T),std::real(z.xy));
                }
            }
        }
    }

    sample_output.close();
//    for (auto i=0; i<1000; i++) {
////        fill the histogram
//        s[i].calc_params();
//        for (auto T : periods) {
//            auto z = s[i](T);
//            zhist(log10(T),std::real(z.xy)/std::abs(z.xy));
//        }
//    }
    std::ofstream histogram_file;
    histogram_file.open(base_filename + "data_fit.res");
    int linecount = 0;
    for (auto &&x : boost::histogram::indexed(zhist)) {
        auto t_hist = (x.bin(0).upper() - x.bin(0).lower()) * 0.5 + x.bin(0).lower();
        auto z_hist = (x.bin(1).upper() - x.bin(1).lower()) * 0.5 + x.bin(1).lower();
        histogram_file << t_hist << " " << z_hist << " " << *x << "\n";
        linecount++;
        if (linecount % 2001 == 0) {
            histogram_file << "\n";
        }
    }
    histogram_file.close();
    histogram_file.open("real_xy.txt");
    for (auto D : d){
        histogram_file << std::setprecision(3) << std::setw(15) << log10(D.first) << std::setw(15)
                       //                       << std::real(D.second.xy)/std::abs(D.second.xy) << "\n";
                       << std::real(D.second.xy) << std::setw(15) << cov[D.first] << "\n";
    }
    histogram_file.close();
    std::cout << "number of model processed: " << n_m_processed <<"\n";
//
//
//    std::cout << "=========================================================================\n";
//    std::cout << s[0]<<std::endl;
//    std::cout << "=========================================================================\n";
//    std::cout << s[1]<<std::endl;
//    std::cout << "=========================================================================\n";

    return 0;
}

















/* parsing program options and config */

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "Display this message and exit.")
            ("config,C", po::value<std::string>()->default_value("MTd1ASampler.cfg"), "Specify configuration file path.");
    po::store(po::parse_command_line(argc, argv, generic), p_vm);
    po::notify(p_vm);
    return generic;
}
boost::program_options::options_description parse_config(boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description config("Configuration");
    config.add_options()
            // Algorithm section
            ("n-interface-max", po::value<int>(), "Maximum number of interfaces in the partition model. It MUST include the earth/air interface.")
            ("n-sigma-bins", po::value<int>(), "Number of bins to be used in sigma-mean discretization, in log space (for integrals and histograms).")
            ("n-z-bins", po::value<int>(), "Number of bins to be used in depth discretization (for integrals and histograms).")
            ("n-temperatures", po::value<int>(), "Number of temperatures to be used for parallel tempering."
                                                 "WARNING: the sampler uses two independent chains at temperature T=1 to check for MC convergence, "
                                                 "therefore the total number of parallel chains used in a sampler run will be n-temperatures + 1.")
            ("max-temperature", po::value<double>(), "Value for T_max, must be greater than 1.")
            ("max-depth", po::value<double>(), "Value for the max depth for the deeper interface.")
            ("n-max-iterations", po::value<int>(), "Maximum number of iterations. This number includes burn-in. "
                                                   "If the flag --Converge is used, this value will be ignored. WARNING: convergence might be slow or impossible.")
            ("Infer-burn-in", "Burn-in phase terminates when all the chains are producing logL close to their expected values.")
            ("n-burn-in-iterations", po::value<int>(), "Maximum number of iterations to be used in the burn-in phase. If the flag --Infer-burn-in is used than this value will be ignored.")
            ("random-seed", po::value<int>(), "Seed to initialize random engine.")
            ("n-iterations-between-pt-swaps", po::value<int>(), "Number of iterations between two subsequent parallel tempering swaps.")
            ("n-iter-between-convergence-checks", po::value<int>(), "Number of iterations between two subsequent tests for convergence.")
            ("significance", po::value<double>(), "Significance value for the convergence test.")
            ("n1-subsample", po::value<int>(), "Sub-sample interval for Sample 1")
            ("n2-subsample", po::value<int>(), "Sub-sample interval for Sample 2")
            ("min-convergence-ratio", po::value<double>()->default_value(0.90), "Ratio of histograms that must results compatible according to the ks statistics to ensure convergence.\nUnused if the code runs without the -v flag." )
            // Distribution section
            ("prior-min-sigma-mean",po::value<double>(), "lower bound for log(sigma-mean).")
            ("prior-max-sigma-mean",po::value<double>(), "upper bound for log(sigma-mean).")
            ("prior-min-sigma-ratio",po::value<double>(), "lower bound for log(sigma_low/sigma_high). Must be negative.")
            ("prior-max-sigma-ratio",po::value<double>()->default_value(0), "upper bound for log(sigma_low/sigma_high). Must be less or equal to zero. Ideally exactly 0.")
            ("prior-min-beta-strike",po::value<double>()->default_value(-90), "lower bound for beta_strike, in Deg.")
            ("prior-max-beta-strike",po::value<double>()->default_value(90), "upper bound for beta_strike, in Deg.")
            ("proposal-scale", po::value<double>(), "The proposal is defined as a normal distribution with mean = current-value and std = (prior-upper-bound - prior-lower-bound) / scale")
            // Receipt weights
            ("perturb", po::value<double>(), "Probability weight for perturbation swap move.")
            ("birth", po::value<double>(), "Probability weight for birth move.")
            ("death", po::value<double>(), "Probability weight for death move.")
            ("iso-switch", po::value<double>(), "Probability weight for isotropy/anisotropy switch swap move.")
            // Filenames
            ("base-filename", po::value<std::string>(), "Base input/output file name (with path).")
        // initial model

            ;
    if(not p_vm.count("help")){
        po::store(po::parse_config_file(p_vm["config"].as<std::string>().c_str(), config),p_vm);
    }
    po::notify(p_vm);
    return config;
}