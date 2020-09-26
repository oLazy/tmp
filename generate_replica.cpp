//
// Created by Eric Mandolesi on 31/05/2020.
//
#include "objects.h"
#include "global.h"
#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>

int main(int argc, char* argv[]){
    const double stdev{0.03};
    std::ifstream is(argv[1]);
    std::string model_filename{argv[1]};
    std::string out_f = model_filename.replace(model_filename.length()-4,4,"_rep.dat");
//    std::cout << out_f << "\n";
    std::string out_file_name{out_f};
    boost::archive::text_iarchive ia(is);
    mtobj::model m;
    ia >> m;
    m.calc_params();
    std::vector<double> per, rxx, rxy, ryx, ryy;
    std::vector<double> ixx, ixy, iyx, iyy;
    std::vector<MTTensor> z, z1;
    std::vector<double> nrxx, nrxy, nryx, nryy;
    std::vector<double> nixx, nixy, niyx, niyy;
    std::vector<double> errorbar;
    std::map<double, MTTensor> dataset; // period, tensor
    std::map<double, double> cov_0; // period, variance

    for (int iper = -30; iper< 51; ++iper) {
        double T = pow(10., 0.1 * double(iper));
        per.push_back(T);
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
//        std::cout << "T===============================================\n";
//        std::cout << "std:" << z0.maxAbsImpedance() <<"\n";
//        std::cout << sz << "\n" << z0 << "\n" << (sz*z0.maxAbsImpedance()) << "\n";
//        std::cout << "================================================\n";
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
//    gp << "set yrange ["<<min_rxx<<":"<< max_rxx <<"] reverse\n";
//    gp << "set y2range ["<<min_ixx<<":"<< max_ixx <<"] reverse\n";
    gp << "set logscale x\n";
    gp << "set ytics rotate by 45\n";
    gp << "set pointsize 0.25\n";
    gp << "set ytics mirror\n";
    gp << "set multiplot layout 2,2\n";

    gp << "set margins 7, 1, 3, 1\n";
    gp << "set format x ''; unset xlabel\n";
    gp << "plot "<<
       "'-' with lines lc rgb '#F0BC42' title 're(z_{yx})', " << // with points pt 7
       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{yx})', " << //pt 6
       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'\n"; //pt 6
//    gp << "plot '-' with lines lc rgb '#F0BC42' title 're(z_{xx})', " << // with points pt 7
//       "'-' axis x1y2 with lines lc rgb '#8E1F2F' title 'im(z_{xx})'\n"; //pt 6
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
       "'-' with lines lc rgb '#F0BC42' title 're(z_{yx})', " << // with points pt 7
       "'-' with lines lc rgb '#8E1F2F' title 'im(z_{yx})', " << //pt 6
       "'-' with errorbars pt 6 lc rgb '#F0BC42' title '-', " << // with points pt 7
       "'-' with errorbars pt 6 lc rgb '#8E1F2F' title '-'"<<
       "\n"; //pt 6
//    gp << "plot '-' with lines lc rgb '#F0BC42' title 're(z_{yy})', " << // with points pt 7
//       "'-' axis x1y2 with lines lc rgb '#8E1F2F' title 'im(z_{yy})'\n"; //pt 6
    gp.send1d(boost::make_tuple(per,ryy));
    gp.send1d(boost::make_tuple(per,iyy));
    gp.send1d(boost::make_tuple(per,nryy,errorbar));
    gp.send1d(boost::make_tuple(per,niyy,errorbar));
//    gp << "set margins 1, 1, 1, 1\n";

    gp << "unset multiplot\n";


    return 0;
}