//
// Created by Eric Mandolesi on 31/05/2020.
//

#include "objects.h"

#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include "global.h"

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm);
//! inputs: a model and a dataset
//! \param argc
//! \param argv
//! \return plots replica from the model and compares it with dataset. To be used to check forward modeling
int main(int argc, char* argv[]){

    auto desc = parse_cmdline(argc, argv, vm);
    if(vm.count("help")){
        std::cout << desc <<std::endl;
        return 0;
    }

    auto model_file = vm["model-filename"].as<std::string>();
    auto data_file = vm["data-filename"].as<std::string>();
    mtobj::model model = mtobj::io::load(model_file);
    mtobj::Dataset pred_data;
    mtobj::io::dataset meas_data;
    meas_data = mtobj::io::loadDataset(data_file);
    for (auto it : meas_data.d){
        double T = it.first;
        pred_data[T] = model(T);
    }
    std::string out_model_file = "model_plot.txt";
    std::string out_data_file = "data_plot.txt";
    std::string out_replica_file = "replica_plot.txt";
    std::ofstream os;
    os.open(out_replica_file);
    gp_utils::dataset2disk(pred_data,os,true);
    os.close();
    os.open(out_data_file);
    gp_utils::dataset2disk(meas_data,os,true);
    os.close();

//    mtobj::Dataset dataset;
//    mtobj::Cov1 cov;
//    std::cout << "Model file: " << model_file << "\nData file: " << data_file <<"\n";
//    {
//        std::ifstream is;
//        is.open(model_file);
//        boost::archive::binary_iarchive ia(is);
//        ia >> m;
//        is.close();
//    }
//    {
//        std::ifstream is;
//        is.open(data_file);
//        boost::archive::binary_iarchive ia(is);
//        ia >> dataset >> cov;
//        is.close();
//    }
//    m.calc_params();
//    std::vector<double> rxxd, rxyd, ryxd, ryyd;
//    std::vector<double> ixxd, ixyd, iyxd, iyyd;
//
//    std::vector<double> rxx0, rxy0, ryx0, ryy0;
//    std::vector<double> ixx0, ixy0, iyx0, iyy0;
////    std::vector<double> nrxx, nrxy, nryx, nryy;
//    std::vector<double> exx, exy, eyx, eyy;
//    std::vector<double>  per;
//
//    for (auto d:dataset){
//        double T=d.first;
//        auto z0 = m(T);
//        auto zd = d.second;
//        std::cout << "T: " << T << "\n" << z0 << "\n" << zd << "\n";
//        per.push_back(T);
//        rxxd.push_back(std::real(zd.xx));
//        rxyd.push_back(std::real(zd.xy));
//        ryxd.push_back(std::real(zd.yx));
//        ryyd.push_back(std::real(zd.yy));
//        ixxd.push_back(std::imag(zd.xx));
//        ixyd.push_back(std::imag(zd.xy));
//        iyxd.push_back(std::imag(zd.yx));
//        iyyd.push_back(std::imag(zd.yy));
//
//        rxx0.push_back(std::real(z0.xx));
//        rxy0.push_back(std::real(z0.xy));
//        ryx0.push_back(std::real(z0.yx));
//        ryy0.push_back(std::real(z0.yy));
//        ixx0.push_back(std::imag(z0.xx));
//        ixy0.push_back(std::imag(z0.xy));
//        iyx0.push_back(std::imag(z0.yx));
//        iyy0.push_back(std::imag(z0.yy));
//
//        exx.push_back(std::sqrt(std::real(cov[T].xx)));
//        exy.push_back(std::sqrt(std::real(cov[T].xy)));
//        eyx.push_back(std::sqrt(std::real(cov[T].yx)));
//        eyy.push_back(std::sqrt(std::real(cov[T].yy)));
//    }
//
//    //plot
//    Gnuplot gp;
//    auto gp_filename = data_file.replace(data_file.length()-3,3,"pdf");
//    gp << "set term pdfcairo size 4,3 font \"Times,9\"\n";
//    gp << "set tics out\n";
//    gp << "set format x \"10^{%T}\"\n";
//    double Tmin = std::log10(*std::min_element(per.begin(), per.end()));
//    double Tmax = std::log10(*std::max_element(per.begin(), per.end()));
//    gp << "set out \"" << gp_filename << "\"\n";
//    gp << "set xrange[" << pow(10.,Tmin)*0.5 << ":" << pow(10.,Tmax)*2. << "]" <<"\n"; // reverse\n";
//    gp << "set logscale x\n";
//    gp << "set ytics rotate by 45\n";
//    gp << "set pointsize 0.25\n";
//    gp << "set ytics mirror\n";
//    gp << "set multiplot layout 2,2\n";
//
//    gp << "set margins 7, 1, 3, 1\n";
//    gp << "set format x ''; unset xlabel\n";
//    gp << "plot "<<
//       "'-' with lines lc rgb '#F0BC42' title 're(z_{xx})', " << // with points pt 7
//       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{xx})', " << //pt 6
//       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
//       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'\n"; //pt 6
//    gp.send1d(boost::make_tuple(per,rxx0));
//    gp.send1d(boost::make_tuple(per,ixx0));
//    gp.send1d(boost::make_tuple(per,rxxd,exx));
//    gp.send1d(boost::make_tuple(per,ixxd,exx));
//
//    gp << "set margins 7, 1, 3, 1\n";
//    gp << "plot "<<
//       "'-' with lines lc rgb '#F0BC42' title 're(z_{xy})', " << // with points pt 7
//       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{xy})', " << //pt 6
//       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " <<
//       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'"<<
//       "\n";
//    gp.send1d(boost::make_tuple(per,rxy0));
//    gp.send1d(boost::make_tuple(per,ixy0));
//    gp.send1d(boost::make_tuple(per,rxyd, exy));
//    gp.send1d(boost::make_tuple(per,ixyd, exy));
//
//    gp << "set format x \"10^{%T}\"\n";
//    gp << "set margins 7, 1, 3, 1\n";
//    gp << "plot "<<
//       "'-' with lines lc rgb '#F0BC42' title 're(z_{yx})', " << // with points pt 7
//       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{yx})', " << //pt 6
//       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
//       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'"<<
//       "\n"; //pt 6
//    gp.send1d(boost::make_tuple(per,ryx0));
//    gp.send1d(boost::make_tuple(per,iyx0));
//    gp.send1d(boost::make_tuple(per,ryxd,eyx));
//    gp.send1d(boost::make_tuple(per,iyxd,eyx));
//
//    gp << "set margins 7, 1, 3, 1\n";
//    gp << "plot "<<
//       "'-' with lines lc rgb '#F0BC42' title 're(z_{yy})', " << // with points pt 7
//       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{yy})', " << //pt 6
//       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
//       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'"<<
//       "\n"; //pt 6
//    gp.send1d(boost::make_tuple(per,ryy0));
//    gp.send1d(boost::make_tuple(per,iyy0));
//    gp.send1d(boost::make_tuple(per,ryyd,eyy));
//    gp.send1d(boost::make_tuple(per,iyyd,eyy));
//    gp << "unset multiplot\n";
//
//    //summary
//    std::cout <<
//              "ELog(L): " << mtobj::expectedLogL(dataset, cov) <<
//              "\nLog[L(m*)]: " << mtobj::logL(m, dataset, cov) <<
//              "\nstd(ELogL): " << mtobj::stdLogL(dataset) <<
//              "\nvar(ELogL): " << mtobj::varLogL(dataset) << "\n";
    return 0;
}

boost::program_options::options_description parse_cmdline(int argc, char *argv[], boost::program_options::variables_map& p_vm){
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "display this message and exit.")
            ("model-filename,m", po::value<std::string>(), "Model file path.")
            ("data-filename,d", po::value<std::string>(), "Data file path. Cov expected as tensor");
    po::store(po::parse_command_line(argc, argv, generic), p_vm);
    po::notify(p_vm);
    return generic;
}