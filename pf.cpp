//
// Created by eric on 20/03/2021.
//
#include "objects.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include "global.h"

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm);
boost::program_options::options_description parse_config(boost::program_options::variables_map& p_vm);

int main(int argn, char* argv[]) {

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
    auto isGeneratingConfig = p_vm["init-config"].as<bool>();
    if(not isGeneratingConfig and not p_vm.count("help")){
        po::store(po::parse_config_file(p_vm["config"].as<std::string>().c_str(), config),p_vm);
    }
    po::notify(p_vm);
    return config;
}