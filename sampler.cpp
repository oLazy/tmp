//
// Created by Eric Mandolesi on 03/06/2020.
//
#include "objects.h"
#include <string>
#include <chrono>
#include <boost/timer/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
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

#define _RJMCMCFUNC

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
    const int n_iter_between_convergence_checks{vm["n-iter-between-convergence-checks"].as<int>()};
    const double significance{vm["significance"].as<double>()};
    const double min_convergence_ratio{vm["min-convergence-ratio"].as<double>()};
    const int subs1{vm["n1-subsample"].as<int>()};
    const int subs2{vm["n2-subsample"].as<int>()};
    int n_sample_collected{0};
    int n_sample_collected1{0};
    int n_sample_collected2{0};
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
    Buffer outbuffer{static_cast<unsigned long>(n_iter_in_pt)};
    unsigned long iBuffer{0};
    std::ofstream sample_output(base_filename+"_sample_out.bin");
    boost::archive::binary_oarchive sample_oa(sample_output);

    initPrior({0.,max_depth},
              {vm["prior-min-sigma-mean"].as<double>(),vm["prior-max-sigma-mean"].as<double>()},
              {vm["prior-min-sigma-ratio"].as<double>(),vm["prior-max-sigma-ratio"].as<double>()},
              {vm["prior-min-beta-strike"].as<double>(),vm["prior-max-beta-strike"].as<double>()});
    initProposal({0.,max_depth},
                 {vm["prior-min-sigma-mean"].as<double>(),vm["prior-max-sigma-mean"].as<double>()},
                 {vm["prior-min-sigma-ratio"].as<double>(),vm["prior-max-sigma-ratio"].as<double>()},
                 {vm["prior-min-beta-strike"].as<double>(),vm["prior-max-beta-strike"].as<double>()},
                 vm["proposal-scale"].as<double>());
//    Cov0 cov;
    Cov1 cov;
    Dataset d;
    std::ifstream is(base_filename+"_rep.dat");
    boost::archive::text_iarchive ia(is);
    ia >> d >> cov;
    is.close();
    // init model zero.
    model m;
    m.nodes.push_back({0,-3});
    m.nodes.push_back({3000,-3});
    model ml_model; // {MAXIMUM LIKELIHOOD MODEL}
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
    ml_model = m;
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
    auto h_sigmaMean = boost::histogram::make_histogram(sm_ax1, sm_ax2);
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
    auto h_sigmaRatio1 = boost::histogram::make_histogram(sr_ax1,sr_ax2);
    auto h_sigmaRatio2 = boost::histogram::make_histogram(sr_ax1,sr_ax2);
    // BETA STRIKE
    auto bs_ax1 = boost::histogram::axis::regular<>(n_sigma_bins,
                                                    prior[paramType::beta].first,
                                                    prior[paramType::beta].second, "beta");
    auto bs_ax2 = boost::histogram::axis::regular<>(n_z_bins,
                                                    prior[paramType::depth].first,
                                                    prior[paramType::depth].second, "depth");
    auto h_betaStrike = boost::histogram::make_histogram(bs_ax1,bs_ax2);
    auto h_betaStrike1 = boost::histogram::make_histogram(bs_ax1,bs_ax2);
    auto h_betaStrike2 = boost::histogram::make_histogram(bs_ax1,bs_ax2);

    auto logl_ax = boost::histogram::axis::regular<>(101,el-7*sl,el+7*sl,"logL");
    auto hll = boost::histogram::make_histogram(logl_ax);
    auto hll1 = boost::histogram::make_histogram(logl_ax);
    auto hll2 = boost::histogram::make_histogram(logl_ax);
    // anis probability===============================================================================================
    auto anisax = boost::histogram::axis::regular<>(n_z_bins,
                                                    prior[paramType::depth].first,
                                                    prior[paramType::depth].second, "depth");
    auto h_anis = boost::histogram::make_histogram(anisax);
    auto h_anis1 = boost::histogram::make_histogram(anisax);
    auto h_anis2 = boost::histogram::make_histogram(anisax);

    //================================================================================================================
    auto ni_ax = boost::histogram::axis::regular<>(max_interfaces,1,max_interfaces,"interfaces");
    auto h_n_inter = boost::histogram::make_histogram(ni_ax); // TODO: update this histogram to save k in spite of #interfaces
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


    // START THE ALGORITHM //

    for (auto itern=1; itern!=itern_max;itern++) {
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
                case move::perturb:
                    rjmcmc::perturb(iic, m, d, cov, chains, proposed, accepted, timer);
                    break;
                case move::birth:
                    rjmcmc::birth(iic, m, d, cov, chains, max_interfaces, proposed, accepted, timer);
                    break;
                case move::death:
                    rjmcmc::death(iic, m, d, cov, chains, proposed, accepted, timer);
                    break;
                case move::iso_switch:
                    rjmcmc::isoswap(iic, m, d, cov, chains, proposed, accepted, timer);
                    break;
            }
            if(m.beta == 1){
                if(m.logL > ml_model.logL) ml_model = m; // save maximum likelihood model
            }
//            // fill the histogram
//            ▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
//            ██░▄▄▄██▄██░██░████░█████▄██░▄▄█▄░▄█░▄▄
//            ██░▄▄███░▄█░██░████░▄▄░██░▄█▄▄▀██░██▄▄▀
//            ██░████▄▄▄█▄▄█▄▄███▄██▄█▄▄▄█▄▄▄██▄██▄▄▄
//            ▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀

            if(status != SamplerStatus::burn_in) {// only collect samples if the burn-in phase is over
                if (m.beta == 1) { // only collect samples from chains sampling at temperature t=1
                    n_sample_collected++;
                    iBuffer++;
                    h_n_inter(m.nodes.size());
                    hll(m.logL);
                    if ((target == Target::sample1) && (itern % subs1 == 0)) {
                        hll1(m.logL);
                        n_sample_collected1++;
                    }
                    if ((target == Target::sample2) && (itern % subs2 == 0)) {
                        hll2(m.logL);
                        n_sample_collected2++;
                    }
                    outbuffer.push_back(m);
                    if(iBuffer==outbuffer.capacity()){// buffer is full
                        for(auto outm=outbuffer.begin();outm!=outbuffer.end();outm++){
                            sample_oa << *outm;}
                        iBuffer=0;
                    }
                    for (int i = 0; i < n_z_bins; i++) {
                        auto dz = (prior[paramType::depth].second - prior[paramType::depth].first) /
                                  static_cast<double>(n_z_bins);
                        auto z_hist = static_cast<double>(i) * dz + (dz * 0.5);
                        auto sm_hist = m.getNode(z_hist).params[paramType::sigmaMean].getValue();
                        auto sr_hist = m.getNode(z_hist).params[paramType::sigmaRatio].getValue();
                        auto bs_hist = m.getNode(z_hist).params[paramType::beta].getValue();
                        // eventually here I can compute sigma high and sigma low to fill further histograms
                        h_sigmaMean(sm_hist, z_hist);
                        if(not std::isnan(sr_hist)){
                            h_sigmaRatio(sr_hist, z_hist);
                            h_anis(z_hist);
                        }
                        if(not std::isnan(bs_hist)){
                            h_betaStrike(bs_hist, z_hist);
                        }
                        if ((target == Target::sample1) && (itern % subs1 == 0)) {
                            h_sigmaMean1(sm_hist, z_hist);
                            if(not std::isnan(sr_hist)){
                                h_sigmaRatio1(sr_hist, z_hist);
                                h_anis1(z_hist);
                            }
                            if(not std::isnan(bs_hist)){
                                h_betaStrike1(bs_hist, z_hist);
                            }
                        }
                        if ((target == Target::sample2) && (itern % subs2 == 0)) {
                            h_sigmaMean2(sm_hist, z_hist);
                            if(not std::isnan(sr_hist)){
                                h_sigmaRatio2(sr_hist, z_hist);
                                h_anis2(z_hist);
                            }
                            if(not std::isnan(bs_hist)){
                                h_betaStrike2(bs_hist, z_hist);
                            }
                        }

                    }
                }
            }


            if ((itern % 1000 == 0) && (iic==0)) {
                auto t1 = Clock::now();
                std::cout << "reached iter: " << itern << " in " <<
                          std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count() << " seconds. " <<
                          "log(L) = " << m.logL << ". Best log(L) till now = " << ml_model.logL << "\n";
            }

        } // done  one iter for all chains


        /*===============================================================
         * CONVERGENCE TEST SECTION
        =================================================================*/
        if (    (status == SamplerStatus::sampling) &&
                ((itern % n_iter_between_convergence_checks) == 0) ) { // if i am sampling I check the convergence between the samples
            std::cout << "Testing for convergence ...\n";
            std::cout << "itern: " << itern << "\n";
            std::cout << "Sample1 size: " << n_sample_collected1 << "\n";
            std::cout << "Sample2 size: " << n_sample_collected2 << "\n";
            auto test_pass_ratio{0.};
            auto w=1./static_cast<double>(n_z_bins);
            for (auto i = 0; i < n_z_bins; i++) {
                std::vector<double> sample1(n_sigma_bins, 0);
                std::vector<double> sample2(n_sigma_bins, 0);
                for (auto j = 0; j < n_sigma_bins; j++) {
                    sample1[j] = h_sigmaMean1.at(j, i);
                    sample2[j] = h_sigmaMean2.at(j, i);
                }
                if(cvt::ks2test(sample1, sample2, n_sample_collected1, n_sample_collected2, significance)) test_pass_ratio+=w;
                if(test_pass_ratio>=(min_convergence_ratio)){
                    status = SamplerStatus::convergence;
                    break; // don't need other proofs, the two samples are in convergence
                }
            }
            std::cout << "KS test pass ratio: " << test_pass_ratio <<"\n";
            if (status == SamplerStatus::convergence) {
                std::cout << "Sample1 and Sample2 are in convergence after " << itern  << " iterations.\n";
                if(vm.count("Converge"))break;
            }
        }

        if(status==SamplerStatus::burn_in) {    // if I am in burn-in, check if burn-in is over
            for (auto iic = 0; iic < chains.size(); iic++) {
                if (isLogLhoodExpected(chains[iic].logL, el, sl, chains[iic].beta, 1)) {
                    status = SamplerStatus::sampling;
                } else {
                    status = SamplerStatus::burn_in;
                    break;
                }
            }
            if(!vm.count("Infer-burn-in")){
                if(itern<=burn_in_n){
                    status = SamplerStatus::burn_in;}
                else{
                    status = SamplerStatus::sampling;}
            }
            if(status==SamplerStatus::sampling) std::cout << "My principle states that burn-in is done after " << itern << " iterations.\n";
        }

        /*===============================================================
         * PARALLEL TEMPERING SECTION
        =================================================================*/
        if (itern % n_iter_in_pt == 0) { // propose exchange between chains
            parallel_tempering_swap(chains);}
        // ================================================================== //
    }
// print dataset statistics
    std::cout << "\n===================================\n" << m << "\n===================================\n";
    std::cout <<
              "ELog(L): " << mtobj::expectedLogL(d, cov) <<
              "\nLog[L(m*)]: " << mtobj::logL(chains[0], d, cov) <<
              "\nstd(ELogL): " << mtobj::stdLogL(d) <<
              "\nvar(ELogL): " << mtobj::varLogL(d) << "\n";

    if(status!=SamplerStatus::burn_in) {
        // print histogram data if burn-in is over
        gp_utils::d2hist2disk(h_sigmaMean, base_filename + "hist_gp_mean_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_sigmaMean, base_filename + "norm.res", n_sigma_bins, true);
        gp_utils::d2hist2disk(h_sigmaMean1, base_filename + "hist1_gp_mean_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_sigmaMean2, base_filename + "hist2_gp_mean_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_sigmaRatio, base_filename + "hist_gp_ratio_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_sigmaRatio1, base_filename + "hist1_gp_ratio_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_sigmaRatio2, base_filename + "hist2_gp_ratio_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_betaStrike, base_filename + "hist_gp_strike_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_betaStrike1, base_filename + "hist1_gp_strike_data.res", n_sigma_bins);
        gp_utils::d2hist2disk(h_betaStrike2, base_filename + "hist2_gp_strike_data.res", n_sigma_bins);

        std::ofstream hll_file(base_filename + "logL.res");
        std::ofstream hin_file(base_filename + "hInter.res");
        for (auto &&x : boost::histogram::indexed(hll)) {
            hll_file << 0.5 * (x.bin(0).lower() + x.bin(0).upper()) << " " << *x << "\n";
        }
        hll_file.close();
        for (auto &&x : boost::histogram::indexed(h_n_inter)) {
            hin_file << 0.5 * (x.bin(0).lower() + x.bin(0).upper()) << " " << *x << "\n";
        }
        hin_file.close();
    }
    // OUTPUT ML MODEL
    gp_utils::model2disk(ml_model, paramType::sigmaMean, prior, "ml_sigma_mean.res");
    gp_utils::model2disk(ml_model, paramType::sigmaRatio, prior, "ml_sigma_ratio.res");
    gp_utils::model2disk(ml_model, paramType::beta, prior, "ml_beta_strike.res");

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
    sample_output.close();

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
        return 0;}
    catch (...){
        return -1;
    }
}