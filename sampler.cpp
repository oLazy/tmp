//
// Created by Eric Mandolesi on 03/06/2020.
//
#include "objects.h"
#include <string>
#include <chrono>
#include <boost/timer/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <iomanip>
#include <boost/program_options.hpp>
#include "global.h"

typedef std::chrono::high_resolution_clock Clock;
using namespace mtobj;
enum class move{perturb,birth,death,iso_switch};

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm);
boost::program_options::options_description parse_config(boost::program_options::variables_map& p_vm);

int main(int argn, char* argv[]) {


    auto desc = parse_cmdline(argn, argv, vm);
    auto config_desc = parse_config(vm);
    if(vm.count("help")){
        std::cout << desc <<std::endl;
        std::cout << config_desc << std::endl;
        return 0;
    }
    auto conf_filename = vm["config"].as<std::string>();
    auto isBuildingConfig = vm["init_config"].as<bool>();
    if(isBuildingConfig) {
        std::ofstream os;
        os.open(conf_filename, std::ios::trunc);
// Algorithm section
        os << "n-interface-max=8" << "\n";
        os << "n-sigma-bins=512" << "\n";
        os << "n-z-bins=1024" << "\n";
        os << "n-temperatures=7" << "\n";
        os << "max-temperature=1000." << "\n";
        os << "max-depth=10000." << "\n";
        os << "n-max-iterations=300000" << "\n";
        os << "n-burn-in-iterations=30000" << "\n";
        os << "random-seed=23" << "\n";
       // Distribution section
        os << "prior-min-sigma-mean=-5" << "\n";
        os << "prior-max-sigma-mean=2" << "\n";
        os << "prior-min-sigma-ratio=-3" << "\n";
        os << "prior-max-sigma-ratio=0" << "\n";
        os << "prior-min-beta-strike=-90" << "\n";
        os << "prior-max-beta-strike=90" << "\n";
        os << "proposal-scale=20" << "\n";
        // Receipt weights
        os <<"perturb=0.7" << "\n";
        os <<"birth=0.1" << "\n";
        os <<"death=0.1" << "\n";
        os <<"iso-switch=0.1" << "\n";
        return 0;
    }
    std::cout << vm["iso-switch"].as<double>() << "\n";
    std::cout << "running\n";
    return 0;

    if (argn != 2) {
        std::cout << "Error! this program requires a single argument.\n";
        return 1;
    }

    std::string base_filename{argv[1]};
    // init
    const int max_interfaces{8};
    const int n_sigma_bins{512};
    const int n_z_bins{1024};
    const int n_temp{8};
    const double max_temp{1000.};
    const double max_depth{8000.};
    const int itern_max{300000};
    const int burn_in_n{30000};

    std::array<std::map<std::string , unsigned long>,n_temp> proposed, accepted;
    for (int iic=0; iic<n_temp;iic++) {
        proposed[iic]["perturb"] = 0;
        proposed[iic]["birth"] = 0;
        proposed[iic]["death"] = 0;
        proposed[iic]["iso_switch"] = 0;
        accepted[iic]["perturb"] = 0;
        accepted[iic]["birth"] = 0;
        accepted[iic]["death"] = 0;
        accepted[iic]["iso_switch"] = 0;
    }

    initPrior({0.,max_depth}, {-5,2}, {-3,0}, {-90,90});
    initProposal({0.,max_depth}, {-5,2}, {-3,0}, {-90,90},100.);
    Cov0 cov;
    Dataset d;
    std::ifstream is(base_filename+"_rep.dat");
    boost::archive::text_iarchive ia(is);
//    boost::archive::binary_iarchive ia(is);
    ia >> d >> cov;
    is.close();
    // init model zero. it is a conductive crust lying over a resistive basement
    model m;
    m.nodes.push_back({0,-3});
    m.nodes.push_back({3000,1});


    if(!m.isInPrior()){
        std::cerr << "model not in prior\n";
        return 17;
    }
    m.calc_params();
    m.setLogL(mtobj::logL(m, d, cov));
    std::vector<model> chains; // chain set
    calc_beta(n_temp, max_temp, m, chains);
//    { // make tempering schedule
//        double length = log10(max_temp);
//        double delta = length / static_cast<double>(n_temp-1);
//        for (int i = 0; i < n_temp; i++) {
//
//
//            double temperature = pow(10,static_cast<double>(i)*delta);
//            double beta = pow(temperature,-1);
//            std::cout << "t[" << i << "]: " << temperature << "\n";
//            m.setBeta(beta);
//            chains.push_back(m);
//        }
//    }
    m = chains[0]; // t=1
    // log dataset statistics

    auto el = mtobj::expectedLogL(d, cov);
    auto sl = mtobj::stdLogL(d);
    std::cout <<
              "ELog(L): " << mtobj::expectedLogL(d, cov) <<
              "\nLog[L(m0)]: " << mtobj::logL(m, d, cov) <<
              "\nstd(ELogL): " << mtobj::stdLogL(d) <<
              "\nvar(ELogL): " << mtobj::varLogL(d) << "\n";
    // create histogram
    auto sm_ax1 = boost::histogram::axis::regular<>(n_sigma_bins,
                                                    prior[paramType::sigmaMean].first,
                                                    prior[paramType::sigmaMean].second, "sigma");
    auto sm_ax2 = boost::histogram::axis::regular<>(n_z_bins,
                                                    prior[paramType::depth].first,
                                                    prior[paramType::depth].second, "depth");
    auto hg = boost::histogram::make_histogram(sm_ax1,sm_ax2);

    auto logl_ax = boost::histogram::axis::regular<>(13,el-3*sl,el+3*sl,"logL");
    auto hll = boost::histogram::make_histogram(logl_ax);
    auto ni_ax = boost::histogram::axis::regular<>(16,1,16,"interfaces");
    auto h_n_inter = boost::histogram::make_histogram(ni_ax);

    auto t0 = Clock::now();
    boost::timer::cpu_timer timer;
    timer.stop();
    for (auto itern=0; itern<itern_max;itern++) {

        for (int iic = 0; iic < n_temp; iic++) {
            m = chains[iic];
            move mt;
            double move_n = urn(gen);
            //====================================================
            //receipt//
            //====================================================
            if (move_n <= 0.7) mt = move::perturb;
            if (move_n > 0.7 && move_n <= 0.8) mt = move::birth;
            if (move_n > 0.8 && move_n <= 0.9) mt = move::death;
            if (move_n > 0.9 && move_n <= 1.0) mt = move::iso_switch;
            //====================================================

            switch (mt) {
                case move::perturb: {
//                    boost::timer::auto_cpu_timer tm;
                    for (int n = 0; n < m.nodes.size(); n++) { // perturb each parameter independently
                        for (int pt = paramType::begin; pt != paramType::end; pt++) {
                            if (m.nodes[n].params[pt].isActive()) {
                                proposed[iic]["perturb"]++;
                                auto m1 = perturb(m, n, static_cast<paramType>(pt));
                                if (m1.isInPrior() && m1.isValid()) {
                                    m1.calc_params();
                                    timer.resume();
                                    m1.setLogL(logL(m1, d, cov));
                                    timer.stop();
                                    auto u = urn(gen);
                                    auto l0 = m.logL;
                                    auto l1 = m1.logL;
                                    if (u < pow(exp(l1 - l0), m.beta)) {
                                        m = m1;
                                        chains[iic] = m1;
                                        accepted[iic]["perturb"]++;
                                    }
                                }
                            }
                        }
                    }
//                    std::cout << "perturb: ";
                }
                    break;

                case move::birth: {
                    proposed[iic]["birth"]++;
//                    boost::timer::auto_cpu_timer tm;
                    if (m.nodes.size() < max_interfaces) {
                        auto m1 = birth(m, birthType::any);
                        if (m1.isInPrior()) {
                            m1.calc_params();
                            timer.resume();
                            m1.setLogL(logL(m1, d, cov));
                            timer.stop();
                            auto u = urn(gen);
                            auto l0 = m.logL;
                            auto l1 = m1.logL;
                            if (u < pow(exp(l1 - l0), m.beta)) {
                                m = m1;
                                chains[iic] = m1;
                                accepted[iic]["birth"]++;
                            }
                        }

                    }
//                else {
//                    continue;
//                }
//                    std::cout << "birth: ";
                }
                    break;
                case move::death: {
                    proposed[iic]["death"]++;
//                    boost::timer::auto_cpu_timer tm;
                    if (m.nodes.size() > 1) {
                        auto m1 = death(m);
                        if (m1.isInPrior()) {
                            m1.calc_params();
                            timer.resume();
                            m1.setLogL(logL(m1, d, cov));
                            timer.stop();
                            auto u = urn(gen);
                            auto l0 = m.logL;
                            auto l1 = m1.logL;
                            if (u < pow(exp(l1 - l0), m.beta)) {
                                m = m1;
                                chains[iic] = m1;
                                accepted[iic]["death"]++;
                            }
                        }
                    }
//                else{
//                    continue;
//                }
//                    std::cout << "death:";
                }
                    break;
                case move::iso_switch: {
//                    boost::timer::auto_cpu_timer tm;
                    for (int n = 0; n < m.nodes.size(); n++) { // try to switch each node independently
                        proposed[iic]["iso_switch"]++;
                        auto m1 = iso_switch(m, n);
                        if (m1.isInPrior() && m1.isValid()) {
                            m1.calc_params();
                            timer.resume();
                            m1.setLogL(logL(m1, d, cov));
                            timer.stop();
                            auto u = urn(gen);
                            auto l0 = m.logL;
                            auto l1 = m1.logL;
                            if (u < pow(exp(l1 - l0), m.beta)) {
                                m = m1;
                                chains[iic] = m1;
                                accepted[iic]["iso_switch"]++;
                            }
                        }
                    }
//                    std::cout << "perturb: ";
                }
                    break;
            }

            // fill the histogram
            if ((itern > burn_in_n) && (m.beta == 1)) { // only collect samples from t0 chain
                h_n_inter(m.nodes.size());
                hll(m.logL);
                for (int i = 0; i < n_z_bins; i++) {
                    auto dz = (prior[paramType::depth].second - prior[paramType::depth].first) /
                              static_cast<double>(n_z_bins);
                    auto z_hist = static_cast<double>(i) * dz + (dz * 0.5);
//            auto z_hist = prior[paramType::depth].first;
                    //            auto z_hist = (x.bin(1).upper() - x.bin(1).lower())*0.5 + x.bin(1).lower();
                    auto sm_hist = m.getNode(z_hist).params[paramType::sigmaMean].getValue();
                    hg(sm_hist, z_hist);
//            std::cout << z_hist << "\t" << sm_hist << "\n";
                }
            }

//        return 0;

            if (((itern + 1) % 1000 == 0) && (m.beta == 1)) {
                auto t1 = Clock::now();
                std::cout << "reached iter: " << itern + 1 << "in " <<
                          std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count() << " seconds. " <<
                          "log(L) = " << m.logL << "\n";

            }

        }

        /*===============================================================
         * PARALLEL TEMPERING SECTION
        =================================================================*/
        if ((itern + 1) % 100 == 0) { // propose exchange between chains
            parallel_tempering_swap(n_temp, chains);
//            boost::random::uniform_int_distribution<int> chain_picker(0,n_temp-1);
//            for (auto i = 0; i < n_temp * n_temp; i++){ // propose n^2 switches
//                int j = chain_picker(gen);
//                int k = chain_picker(gen);
//                if(j!=k){
//                    auto bj = chains[j].beta;
//                    auto bk = chains[k].beta;
//                    if(urn(gen) <  pow(exp(chains[j].logL - chains[k].logL),(bk-bj))){
//                        chains[j].setBeta(bk);
//                        chains[k].setBeta(bj);
//                        std::swap(chains[j],chains[k]); // i-th temperature remains associated to the i-th chain
//                    }
//                }
//            }
//        }
        }
    }
// print dataset statistics
    std::cout << "\n===================================\n" << m << "\n===================================\n";
    std::cout <<
              "ELog(L): " << mtobj::expectedLogL(d, cov) <<
              "\nLog[L(m*)]: " << mtobj::logL(chains[0], d, cov) <<
              "\nstd(ELogL): " << mtobj::stdLogL(d) <<
              "\nvar(ELogL): " << mtobj::varLogL(d) << "\n";

    // print histogram data
    std::ofstream histogram_file("hist_gp_data.res");
    int linecount = 0;
    for(auto &&x : boost::histogram::indexed(hg)){
        auto sm_hist = (x.bin(0).upper() - x.bin(0).lower())*0.5 + x.bin(0).lower();
        auto z_hist = (x.bin(1).upper() - x.bin(1).lower())*0.5 + x.bin(1).lower();
        histogram_file << sm_hist << " " << z_hist << " " << *x << "\n";
        linecount++;
        if(linecount%n_sigma_bins==0) histogram_file << "\n";
    }
    histogram_file.close();
    std::ofstream hll_file("logL.res");
    std::ofstream hin_file("hInter.res");
    for (auto &&x : boost::histogram::indexed(hll)){
        hll_file << 0.5*(x.bin(0).lower()+x.bin(0).upper()) << " " << *x << "\n";
    }
    hll_file.close();
    for (auto &&x : boost::histogram::indexed(h_n_inter)){
        hin_file << 0.5*(x.bin(0).lower()+x.bin(0).upper()) << " " << *x << "\n";
    }
    hin_file.close();
    std::cout << "time spent in log(L) subroutine: " << timer.format() << "\n";
#ifdef _OMP
    std::cout << "program run in parallel with OMP using a team of 4 threads.\n";
#else
    std::cout << "program run serial on single thread.\n";
#endif
    std::ofstream faccst(base_filename +"_acc_stat.res");
    for (int iic=0; iic<n_temp; iic++) {
        faccst << "Acceptance statistics for temperature T = " << pow(chains[iic].beta,-1) << "\n";
        faccst << std::setw(12) << "move" << std::setw(12) << "proposed" << std::setw(12) << "accepted"
               << std::setw(12) << "ratio" << std::endl;
        for (const auto& k: proposed[iic]) {
            auto key = k.first;
            auto ratio = static_cast<double>(accepted[iic][key]) / static_cast<double>(proposed[iic][key]);
            faccst    << std::setw(12) << key << std::setw(12) << proposed[iic][key] << std::setw(12)
                      << accepted[iic][key] << std::setw(12) << ratio << std::endl;
        }
        faccst << "\n\n";
    }
    faccst.close();
    return 0;
}

/* parsing program options and config */

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "Display this message and exit.")
            ("config,C", po::value<std::string>()->default_value("MTd1ASampler.cfg"), "Specify configuration file path.")
            ("init_config,c", po::value<bool>()->default_value(false), "Create simple configuration file and exit.")
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
            ("n-burn-in-iterations", po::value<int>(), "Maximum number of iterations to be used in the burn-in phase. If the flag --Infer-burn-in is used than this value will be ignored.")
            ("random-seed", po::value<int>(), "Seed to initialize random engine.")
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

        // initial model

            ;
    auto isGeneratingConfig = p_vm["init_config"].as<bool>();
    if(not isGeneratingConfig and not p_vm.count("help")){
        po::store(po::parse_config_file(p_vm["config"].as<std::string>().c_str(), config),p_vm);
    }
    po::notify(p_vm);
    return config;
}