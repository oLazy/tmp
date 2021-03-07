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
enum class Target{sample1, sample2, none};
enum class SamplerStatus{burn_in, sampling, convergence};
boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm);
boost::program_options::options_description parse_config(boost::program_options::variables_map& p_vm);
int generate_configuration_file(boost::program_options::variables_map& p_vm);

int main(int argn, char* argv[]) {

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
    std::string base_filename{vm["base-filename"].as<std::string>().c_str()};

//                      ▄▄                                         ▄▄
//    ▀████▀            ██   ██                              ██    ██
//      ██                   ██                              ██
//      ██ ▀████████▄ ▀███ ██████     ▄██▀███ ▄▄█▀██ ▄██▀████████▀███   ▄██▀██▄▀████████▄
//      ██   ██    ██   ██   ██       ██   ▀▀▄█▀   ███▀  ██  ██    ██  ██▀   ▀██ ██    ██
//      ██   ██    ██   ██   ██       ▀█████▄██▀▀▀▀▀▀█       ██    ██  ██     ██ ██    ██
//      ██   ██    ██   ██   ██       █▄   ████▄    ▄█▄    ▄ ██    ██  ██▄   ▄██ ██    ██
//    ▄████▄████  ████▄████▄ ▀████    ██████▀ ▀█████▀█████▀  ▀████████▄ ▀█████▀▄████  ████▄
//
//
    SamplerStatus status{SamplerStatus::burn_in}; // begin from random walk optimization
    std::cout << base_filename << std::endl;
    const int max_interfaces{vm["n-interface-max"].as<int>()};
    const int n_sigma_bins{vm["n-sigma-bins"].as<int>()};
    const int n_z_bins{vm["n-z-bins"].as<int>()};
    const int n_temp{vm["n-temperatures"].as<int>()};
    const double max_temp{vm["max-temperature"].as<double>()};
    const double max_depth{vm["max-depth"].as<double>()};
    const int itern_max{vm["n-max-iterations"].as<int>()};
    const int burn_in_n{vm["n-burn-in-iterations"].as<int>()};
    const int n_iter_in_pt{vm["n-iterations-between-pt-swaps"].as<int>()};
    const int n_iter_between_convergence_checks = n_iter_in_pt; // TODO: update this parameter in the parameter file
    const double significance = 0.1; // TODO: update this parameter in the parameter file
    std::vector<std::map<std::string , unsigned long>> proposed(n_temp+1), accepted(n_temp+1);
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

    initPrior({0.,max_depth},
              {vm["prior-min-sigma-mean"].as<double>(),vm["prior-max-sigma-mean"].as<double>()},
              {vm["prior-min-sigma-ratio"].as<double>(),vm["prior-max-sigma-ratio"].as<double>()},
              {vm["prior-min-beta-strike"].as<double>(),vm["prior-max-beta-strike"].as<double>()});
    initProposal({0.,max_depth},
                 {vm["prior-min-sigma-mean"].as<double>(),vm["prior-max-sigma-mean"].as<double>()},
                 {vm["prior-min-sigma-ratio"].as<double>(),vm["prior-max-sigma-ratio"].as<double>()},
                 {vm["prior-min-beta-strike"].as<double>(),vm["prior-max-beta-strike"].as<double>()},
                 vm["proposal-scale"].as<double>());
    Cov0 cov;
    Dataset d;
    std::ifstream is(base_filename+"_rep.dat");
    boost::archive::text_iarchive ia(is);
    ia >> d >> cov;
    is.close();
    // init model zero.
    model m;
    m.nodes.push_back({0,-3});
    m.nodes.push_back({3000,-3});

    if(!m.isInPrior()){
        std::cerr << "model 0 not in prior\n";
        return 17;
    }
    m.calc_params();
    m.setLogL(mtobj::logL(m, d, cov));
    std::vector<model> chains; // chain set. After initialization chain set dimension shall be n_temp + 1
    calc_beta(n_temp, max_temp, m, chains, true);

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
    //================================================================================================================//
    // SIGMA MEAN
    auto sm_ax1 = boost::histogram::axis::regular<>(n_sigma_bins,
                                                    prior[paramType::sigmaMean].first,
                                                    prior[paramType::sigmaMean].second, "sigma");
    auto sm_ax2 = boost::histogram::axis::regular<>(n_z_bins,
                                                    prior[paramType::depth].first,
                                                    prior[paramType::depth].second, "depth");
    auto h_sigmaMean1 = boost::histogram::make_histogram(sm_ax1, sm_ax2);
    auto h_sigmaMean2 = boost::histogram::make_histogram(sm_ax1, sm_ax2);

    // SIGMA RATIO
    auto sr_ax1 = boost::histogram::axis::regular<>(n_sigma_bins,
                                                    prior[paramType::sigmaRatio].first,
                                                    prior[paramType::sigmaRatio].second, "ratio");
    auto sr_ax2 = boost::histogram::axis::regular<>(n_z_bins,
                                                    prior[paramType::depth].first,
                                                    prior[paramType::depth].second, "depth");
    auto h_sigmaRatio = boost::histogram::make_histogram(sr_ax1,sr_ax2);
    // BETA STRIKE
    auto bs_ax1 = boost::histogram::axis::regular<>(n_sigma_bins,
                                                    prior[paramType::beta].first,
                                                    prior[paramType::beta].second, "beta");
    auto bs_ax2 = boost::histogram::axis::regular<>(n_z_bins,
                                                    prior[paramType::depth].first,
                                                    prior[paramType::depth].second, "depth");
    auto h_betaStrike = boost::histogram::make_histogram(bs_ax1,bs_ax2);

    auto logl_ax = boost::histogram::axis::regular<>(101,el-7*sl,el+7*sl,"logL");
    auto hll = boost::histogram::make_histogram(logl_ax);
    // anis probability===============================================================================================
    auto anisax = boost::histogram::axis::regular<>(n_z_bins,
                                                    prior[paramType::depth].first,
                                                    prior[paramType::depth].second, "depth");
    auto h_anis1 = boost::histogram::make_histogram(anisax);
    auto h_anis2 = boost::histogram::make_histogram(anisax);

    //================================================================================================================
    auto ni_ax = boost::histogram::axis::regular<>(max_interfaces,1,max_interfaces,"interfaces");
    auto h_n_inter = boost::histogram::make_histogram(ni_ax);
    //================================================================================================================//


    auto t0 = Clock::now();
    boost::timer::cpu_timer timer;
    timer.stop();
    auto receipt = balanceReceiptWeights(vm);
    auto pl = std::get<0>(receipt);
    auto bl = pl + std::get<1>(receipt);
    auto dl = bl + std::get<2>(receipt);
    auto il = dl + std::get<3>(receipt);
    std::cout << "iso-switch probability limit shall be 1. In this run is " << il << "\n";

    for (auto itern=0; itern<itern_max;itern++) {
        Target target = Target::none;
        for (int iic = 0; iic < chains.size(); iic++) {
            if(iic==0) {
                target = Target::sample1;
            }
            else if(iic==1){
                target = Target::sample2;
            }
            else{
                target = Target::none;
            }

            m = chains[iic];
            move mt;
            double move_n = urn(gen);

            //====================================================
            //receipt//
            //====================================================
            if (move_n <= pl) mt = move::perturb;
            if (move_n > pl && move_n <= bl) mt = move::birth;
            if (move_n > bl && move_n <= dl) mt = move::death;
            if (move_n > dl && move_n <= il) mt = move::iso_switch;
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

            if ((status != SamplerStatus::burn_in) && (m.beta == 1)) { // only collect samples from t0 chain
                h_n_inter(m.nodes.size());
                hll(m.logL);
                for (int i = 0; i < n_z_bins; i++) {
                    auto dz = (prior[paramType::depth].second - prior[paramType::depth].first) /
                              static_cast<double>(n_z_bins);
                    auto z_hist = static_cast<double>(i) * dz + (dz * 0.5);
//            auto z_hist = prior[paramType::depth].first;
                    //            auto z_hist = (x.bin(1).upper() - x.bin(1).lower())*0.5 + x.bin(1).lower();
                    auto sm_hist = m.getNode(z_hist).params[paramType::sigmaMean].getValue();
                    auto sr_hist = m.getNode(z_hist).params[paramType::sigmaRatio].getValue();
                    auto bs_hist = m.getNode(z_hist).params[paramType::beta].getValue();
                    if(target==Target::sample1) {
                        h_sigmaMean1(sm_hist, z_hist);
                    }else if(target==Target::sample2){
                        h_sigmaMean2(sm_hist, z_hist);
                    }
                    // only add entries to anisotropy images if anisotropy is present
                    if(not std::isnan(sr_hist)){
                        h_sigmaRatio(sr_hist, z_hist);
                    }
                    if(not std::isnan(bs_hist)){
                        h_betaStrike(bs_hist, z_hist);
                    }
//            std::cout << z_hist << "\t" << sm_hist << "\n";
                }
            }

//        return 0;

            if (((itern + 1) % 1000 == 0) && (iic==0)) {
                auto t1 = Clock::now();
                std::cout << "reached iter: " << itern + 1 << "in " <<
                          std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count() << " seconds. " <<
                          "log(L) = " << m.logL << "\n";

            }

        }


        /*===============================================================
         * CONVERGENCE TEST SECTION
        =================================================================*/
        if ((status == SamplerStatus::sampling) && (itern + 1) % n_iter_between_convergence_checks == 0) { // if i am sampling I check the convergence between the samples
            for (auto i = 0; i<n_z_bins;i++){
                std::vector<double> sample1(n_sigma_bins,0);
                std::vector<double> sample2(n_sigma_bins,0);
                for (auto j = 0; j< n_sigma_bins; j++){
                    sample1[j] = h_sigmaMean1.at(j,i);
                    sample2[j] = h_sigmaMean2.at(j,i);
                }
                auto test_res = cvt::chi2twoBins(sample1, sample2);
                if (std::get<cvt::chi2twoBinsResults::significance>(test_res) < significance){

                }
            }
        }

            /*===============================================================
             * PARALLEL TEMPERING SECTION
            =================================================================*/
        if ((itern + 1) % n_iter_in_pt == 0) { // propose exchange between chains
            parallel_tempering_swap(chains);

            if(status==SamplerStatus::burn_in) {    // if I am in burn-in, check if burn-in is over
                for (auto iic = 0; iic < chains.size(); iic++) {
                    if (isLogLhoodExpected(chains[iic].logL, el, sl, chains[iic].beta)) {
                        status = SamplerStatus::sampling;
                    } else {
                        status = SamplerStatus::burn_in;
                        break;
                    }
                }
                if(status!=SamplerStatus::burn_in) std::cout << "My principle states that burn-in is done after " << itern + 1<< " iterations.\n";

                if(!vm.count("Infer-burn-in")){
                    if(itern<=burn_in_n){
                        status = SamplerStatus::burn_in;}
                    else{
                        status = SamplerStatus::sampling;}
                }
            }


        }
        // ================================================================== //

    }
// print dataset statistics
    std::cout << "\n===================================\n" << m << "\n===================================\n";
    std::cout <<
              "ELog(L): " << mtobj::expectedLogL(d, cov) <<
              "\nLog[L(m*)]: " << mtobj::logL(chains[0], d, cov) <<
              "\nstd(ELogL): " << mtobj::stdLogL(d) <<
              "\nvar(ELogL): " << mtobj::varLogL(d) << "\n";

    // print histogram data
    std::ofstream histogram_file;
    histogram_file.open(base_filename+"hist1_gp_mean_data.res");
//    std::ofstream norm_file;
//    norm_file.open(base_filename+"norm_hist1_mean_data.res");
    int linecount = 0;
//    double integral = 0.;
    for(auto &&x : boost::histogram::indexed(h_sigmaMean1)){
        auto sm_hist = (x.bin(0).upper() - x.bin(0).lower())*0.5 + x.bin(0).lower();
        auto z_hist = (x.bin(1).upper() - x.bin(1).lower())*0.5 + x.bin(1).lower();
        histogram_file << sm_hist << " " << z_hist << " " << *x << "\n";
//        integral += (x.bin(0).upper() - x.bin(0).lower()) * (*x);
        linecount++;
        if(linecount%n_sigma_bins==0) {
            histogram_file << "\n";
//            norm_file << integral << "\n";
//            integral = 0.;
        }
    }
    histogram_file.close();

    histogram_file.open(base_filename+"hist2_gp_mean_data.res");
    linecount = 0;
    for(auto &&x : boost::histogram::indexed(h_sigmaMean2)){
        auto sm_hist = (x.bin(0).upper() - x.bin(0).lower())*0.5 + x.bin(0).lower();
        auto z_hist = (x.bin(1).upper() - x.bin(1).lower())*0.5 + x.bin(1).lower();
        histogram_file << sm_hist << " " << z_hist << " " << *x << "\n";
        linecount++;
        if(linecount%n_sigma_bins==0) histogram_file << "\n";
    }
    histogram_file.close();


    histogram_file.open(base_filename+"hist_gp_ratio_data.res");
    linecount = 0;
    for(auto &&x : boost::histogram::indexed(h_sigmaRatio)){
        auto sm_hist = (x.bin(0).upper() - x.bin(0).lower())*0.5 + x.bin(0).lower();
        auto z_hist = (x.bin(1).upper() - x.bin(1).lower())*0.5 + x.bin(1).lower();
        histogram_file << sm_hist << " " << z_hist << " " << *x << "\n";
        linecount++;
        if(linecount%n_sigma_bins==0) histogram_file << "\n";
    }
    histogram_file.close();

    histogram_file.open(base_filename+"hist_gp_strike_data.res");
    linecount = 0;
    for(auto &&x : boost::histogram::indexed(h_betaStrike)){
        auto sm_hist = (x.bin(0).upper() - x.bin(0).lower())*0.5 + x.bin(0).lower();
        auto z_hist = (x.bin(1).upper() - x.bin(1).lower())*0.5 + x.bin(1).lower();
        histogram_file << sm_hist << " " << z_hist << " " << *x << "\n";
        linecount++;
        if(linecount%n_sigma_bins==0) histogram_file << "\n";
    }
    histogram_file.close();

    std::ofstream hll_file(base_filename+"logL.res");
    std::ofstream hin_file(base_filename+"hInter.res");
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


//    ▀███▀▀▀███▀███▄   ▀███▀███▀▀▀██▄     ▀███▀▀▀██▄▀███▀▀▀██▄   ▄▄█▀▀██▄   ▄▄█▀▀▀█▄█▀███▀▀▀██▄       ██     ▀████▄     ▄███▀
//      ██    ▀█  ███▄    █   ██    ▀██▄     ██   ▀██▄ ██   ▀██▄▄██▀    ▀██▄██▀     ▀█  ██   ▀██▄     ▄██▄      ████    ████
//      ██   █    █ ███   █   ██     ▀██     ██   ▄██  ██   ▄██ ██▀      ▀███▀       ▀  ██   ▄██     ▄█▀██▄     █ ██   ▄█ ██
//      ██████    █  ▀██▄ █   ██      ██     ███████   ███████  ██        ███           ███████     ▄█  ▀██     █  ██  █▀ ██
//      ██   █  ▄ █   ▀██▄█   ██     ▄██     ██        ██  ██▄  ██▄      ▄███▄    ▀████ ██  ██▄     ████████    █  ██▄█▀  ██
//      ██     ▄█ █     ███   ██    ▄██▀     ██        ██   ▀██▄▀██▄    ▄██▀██▄     ██  ██   ▀██▄  █▀      ██   █  ▀██▀   ██
//    ▄█████████████▄    ██ ▄████████▀     ▄████▄    ▄████▄ ▄███▄ ▀▀████▀▀   ▀▀███████▄████▄ ▄███▄███▄   ▄████▄███▄ ▀▀  ▄████▄





/* parsing program options and config */

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
int generate_configuration_file(boost::program_options::variables_map& p_vm){
    try{
        auto conf_filename = p_vm["config"].as<std::string>();
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
        os << "n-iterations-between-pt-swaps=1000" << "\n";
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
        os << "base-filename=cg_model_1\n";
        return 0;}
    catch (...){
        return -1;
    }
}