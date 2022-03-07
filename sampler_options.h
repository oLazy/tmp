//
// Created by eric on 11/02/2022.
//

#ifndef MT1DANISMODELPARAMS_SAMPLER_OPTIONS_H
#define MT1DANISMODELPARAMS_SAMPLER_OPTIONS_H




/* parsing program options and config */
#include <boost/program_options.hpp>
boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "Display this message and exit.")
            ("config,C", po::value<std::string>()->default_value("MTd1ASampler.cfg"), "Specify configuration file path.")
            ("init-config,c", po::value<bool>()->default_value(false), "Create simple configuration file and exit.")
            ("Converge,v","Ignore the maximum number of iteration set and terminate the sampling once convergence can be defended (see paper for convergence criterion).")
            ("Infer-burn-in,b","Ignore number of iterations to be used in burn-in phase and infer it from the expected log-likelihood (see paper for end of burn-in criterion).")
            ;
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
            ("m0-filename", po::value<std::string>(), "Initial model input file name (with path).")
            ("data-filename", po::value<std::string>(), "Data file name (with path).")
        // initial model

            ;
    auto isGeneratingConfig = p_vm["init-config"].as<bool>();
    if(not isGeneratingConfig and not p_vm.count("help")){
        po::store(po::parse_config_file(p_vm["config"].as<std::string>().c_str(), config),p_vm);
    }
    po::notify(p_vm);
    return config;
}
int generate_configuration_file(boost::program_options::variables_map& p_vm){
    try{
        auto conf_filename = p_vm["config"].as<std::string>();
        std::ofstream os;
        os.open(conf_filename, std::ios::trunc);
        // Algorithm section
        os << "n-interface-max=8" << "\n";
        os << "n-sigma-bins=512" << "\n";
        os << "n-z-bins=1024" << "\n";
        os << "n-temperatures=15" << "\n";
        os << "max-temperature=1000." << "\n";
        os << "max-depth=10000." << "\n";
        os << "n-max-iterations=300000" << "\n";
        os << "n-burn-in-iterations=30000" << "\n";
        os << "random-seed=23" << "\n";
        os << "n-iterations-between-pt-swaps=100" << "\n";
        os << "n-iter-between-convergence-checks=1000" << "\n";
        os << "significance=0.05" << "\n";
        os << "n1-subsample=2" << "\n";
        os << "n2-subsample=3" << "\n";
        os << "min-convergence-ratio=0.90" << "\n";
        // Distribution section
        os << "prior-min-sigma-mean=-5" << "\n";
        os << "prior-max-sigma-mean=2" << "\n";
        os << "prior-min-sigma-ratio=-3" << "\n";
        os << "prior-max-sigma-ratio=0" << "\n";
        os << "prior-min-beta-strike=-90" << "\n";
        os << "prior-max-beta-strike=90" << "\n";
        os << "proposal-scale=20" << "\n";
        // Receipt weights
        os << "perturb=0.7" << "\n";
        os << "birth=0.1" << "\n";
        os << "death=0.1" << "\n";
        os << "iso-switch=0.1" << "\n";
        // File names
        os << "base-filename=test\n";
        os << "m0-filename=test.bin\n";
        os << "data-filename=test_rep.bin\n";

        return 0;}
    catch (...){
        return -1;
    }
}
int handle_blocking_program_options(int argn, char* argv[]){

    // parse command line and config file
    auto desc = parse_cmdline(argn, argv, vm);
    auto config_desc = parse_config(vm);
    if(vm.count("help")){
        std::cout << desc <<std::endl;
        std::cout << config_desc << std::endl;
        return 0;
    }
    auto isBuildingConfig = vm["init-config"].as<bool>();
    if(isBuildingConfig) {
        return generate_configuration_file(vm);
    }
    return 1;
}
#endif //MT1DANISMODELPARAMS_SAMPLER_OPTIONS_H
