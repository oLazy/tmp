//
// Created by Eric Mandolesi on 31/05/2020.
//
#include "objects.h"
#include "global.h"
#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>

// GLOBAL
// program configuration
boost::program_options::variables_map vm;

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm);
boost::program_options::options_description parse_config(boost::program_options::variables_map& p_vm);

int main(int argc, char* argv[]){
    /* add as config file
     * stdev
     * linear option
     * min T
     * max T
     * No. T
    */
    auto desc = parse_cmdline(argc, argv, vm);
    auto config_desc = parse_config(vm);
    if(vm.count("help")){
        std::cout << desc <<std::endl;
        std::cout << config_desc << std::endl;
        return 1;
    }
    auto conf_filename = vm["config"].as<std::string>();
    auto isBuildingConfig = vm["init_config"].as<bool>();
    if(isBuildingConfig){
        std::ofstream os;
        os.open(conf_filename, std::ios::trunc);
        os << "std_dev="<<0.03<<"\n";
        os << "min_T="<<-2<<"\n";
        os << "max_T="<<3<<"\n";
        os << "n_T="<<81<<"\n";
        os << "model-filename=cg_model_1.bin";
    }

    const double stdev=vm["std_dev"].as<double>();
    auto model_file = vm["model-filename"].as<std::string>();
    std::ifstream is;
    is.open(model_file);

    boost::archive::binary_iarchive ia(is);
    std::string out_f = model_file.replace(model_file.length()-4,4,"_rep.dat");
//    std::cout << out_f << "\n";
    std::string out_file_name{out_f};
//    boost::archive::text_iarchive ia(is);
    mtobj::model m;
    ia >> m;
    m.calc_params();
    std::vector<double> per, per2, rxx, rxy, ryx, ryy;
    std::vector<double> ixx, ixy, iyx, iyy;
    std::vector<MTTensor> z, z1;
    std::vector<double> nrxx, nrxy, nryx, nryy;
    std::vector<double> nixx, nixy, niyx, niyy;
    std::vector<double> errorbar;
    std::map<double, MTTensor> dataset; // period, tensor
    std::map<double, double> cov_0; // period, variance

    auto iTmax = vm["n_T"].as<int>();
    double Tmin, Tmax;
    if(vm["IS"].as<bool>()){
        Tmin = log10(vm["min_T"].as<double>());
        Tmax = log10(vm["max_T"].as<double>());
    }else{
        Tmin = (vm["min_T"].as<double>());
        Tmax = (vm["max_T"].as<double>());
    }
    auto deltaT = (Tmax-Tmin)/static_cast<double>(iTmax-1);
    auto ttmp = Tmin;
    for(int iT=0; iT<iTmax; iT++){
        per.push_back(pow(10.,ttmp));
        double T = pow(10.,ttmp);
        ttmp+=deltaT;
        auto z0 = m(T);
        z.push_back(z0);
        rxx.push_back(std::real(z0.xx));
        rxy.push_back(std::real(z0.xy));
        ryx.push_back(std::real(z0.yx));
        ryy.push_back(std::real(z0.yy));
        ixx.push_back(std::imag(z0.xx));
        ixy.push_back(std::imag(z0.xy));
        iyx.push_back(std::imag(z0.yx));
        iyy.push_back(std::imag(z0.yy));

        MTTensor sz{{rn(gen),rn(gen)},{rn(gen),rn(gen)},{rn(gen),rn(gen)},{rn(gen),rn(gen)}};
        auto ztmp = z0+(sz*z0.maxAbsImpedance()*stdev);
        errorbar.push_back(z0.maxAbsImpedance()*stdev);
        dataset[T] = ztmp;
        cov_0[T] = z0.maxAbsImpedance()*stdev;
        z1.push_back(ztmp);
        nrxx.push_back(std::real(ztmp.xx));
        nrxy.push_back(std::real(ztmp.xy));
        nryx.push_back(std::real(ztmp.yx));
        nryy.push_back(std::real(ztmp.yy));
        nixx.push_back(std::imag(ztmp.xx));
        nixy.push_back(std::imag(ztmp.xy));
        niyx.push_back(std::imag(ztmp.yx));
        niyy.push_back(std::imag(ztmp.yy));
    }
    // print dataset statistics
    std::cout <<
              "ELog(L): " << mtobj::expectedLogL(dataset, cov_0) <<
              "\nLog[L(m*)]: " << mtobj::logL(m, dataset, cov_0) <<
              "\nstd(ELogL): " << mtobj::stdLogL(dataset) <<
              "\nvar(ELogL): " << mtobj::varLogL(dataset) << "\n";

    // save tensors and periods
    std::ofstream os(out_file_name);
    boost::archive::text_oarchive oa(os);
    oa << dataset;
    oa << cov_0;
    // plot!!!
    Gnuplot gp;
    auto gp_filename = out_file_name.replace(out_file_name.length()-3,3,"pdf");
    gp << "set term pdfcairo size 4,3 font \"Times,9\"\n";
    gp << "set tics out\n";
    gp << "set format x \"10^{%T}\"\n";
    gp << "set out \"" << gp_filename << "\"\n";
    gp << "set xrange[" << pow(10.,Tmin)*0.5 << ":" << pow(10.,Tmax)*2. << "] reverse\n";
    gp << "set logscale x\n";
    gp << "set ytics rotate by 45\n";
    gp << "set pointsize 0.25\n";
    gp << "set ytics mirror\n";
    gp << "set multiplot layout 2,2\n";

    gp << "set margins 7, 1, 3, 1\n";
    gp << "set format x ''; unset xlabel\n";
    gp << "plot "<<
       "'-' with lines lc rgb '#F0BC42' title 're(z_{xx})', " << // with points pt 7
       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{xx})', " << //pt 6
       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'\n"; //pt 6
    gp.send1d(boost::make_tuple(per,rxx));
    gp.send1d(boost::make_tuple(per,ixx));
    gp.send1d(boost::make_tuple(per,nrxx,errorbar));
    gp.send1d(boost::make_tuple(per,nixx,errorbar));

    gp << "set margins 7, 1, 3, 1\n";
    gp << "plot "<<
       "'-' with lines lc rgb '#F0BC42' title 're(z_{xy})', " << // with points pt 7
       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{xy})', " << //pt 6
       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " <<
       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'"<<
       "\n";
    gp.send1d(boost::make_tuple(per,rxy));
    gp.send1d(boost::make_tuple(per,ixy));
    gp.send1d(boost::make_tuple(per,nrxy, errorbar));
    gp.send1d(boost::make_tuple(per,nixy, errorbar));

    gp << "set format x \"10^{%T}\"\n";
    gp << "set margins 7, 1, 3, 1\n";
    gp << "plot "<<
       "'-' with lines lc rgb '#F0BC42' title 're(z_{yx})', " << // with points pt 7
       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{yx})', " << //pt 6
       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'"<<
       "\n"; //pt 6
    gp.send1d(boost::make_tuple(per,ryx));
    gp.send1d(boost::make_tuple(per,iyx));
    gp.send1d(boost::make_tuple(per,nryx,errorbar));
    gp.send1d(boost::make_tuple(per,niyx,errorbar));

    gp << "set margins 7, 1, 3, 1\n";
    gp << "plot "<<
       "'-' with lines lc rgb '#F0BC42' title 're(z_{yy})', " << // with points pt 7
       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{yy})', " << //pt 6
       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'"<<
       "\n"; //pt 6
    gp.send1d(boost::make_tuple(per,ryy));
    gp.send1d(boost::make_tuple(per,iyy));
    gp.send1d(boost::make_tuple(per,nryy,errorbar));
    gp.send1d(boost::make_tuple(per,niyy,errorbar));
    gp << "unset multiplot\n";

    return 0;
}

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "display this message and exit.")
            ("IS", po::value<bool>()->default_value(false), "Read min and max periods in seconds.")
            ("config,C", po::value<std::string>()->default_value("MTrgen.cfg"), "Specify configuration file path")
            ("init_config,c", po::value<bool>()->default_value(false), "Create simple configuration file");
    po::store(po::parse_command_line(argc, argv, generic), p_vm);
    po::notify(p_vm);
    return generic;
}

boost::program_options::options_description parse_config(boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description config("Configuration");
    config.add_options()
            ("std_dev", po::value<double>()->default_value(0.03), "Error to be used to simulate "
                                                                  "data.\nError on real and imaginary parts are computed"
                                                                  " independently. Error on each component of the tensor"
                                                                  " element is computed as arg*nrand(0,1)*abs(Z_ij)"
                                                                  "with i, j chosen so that Z_ij is maximum.")
            ("min_T", po::value<double>(), "magnitude of minimum period for the simulated survey. If the "
                                           "--IS is specified, the value must be inserted in seconds.")
            ("max_T", po::value<double>(), "magnitude of maximum period for the simulated survey. If the "
                                            "--IS is specified, the value must be inserted in seconds.")
            ("n_T", po::value<int>()->default_value(81), "Number of periods for the simulated survey, equally log-spaced between min_T and max_T.")
            ("model-filename", po::value<std::string>()->default_value("cg_model_1.bin"), "Model file path.");
    auto isGeneratingConfig = p_vm["init_config"].as<bool>();
    if(not isGeneratingConfig and not p_vm.count("help")){
        po::store(po::parse_config_file(p_vm["config"].as<std::string>().c_str(), config),p_vm);
    }
    po::notify(p_vm);
    return config;
//    std::cout << config << std::endl;
}