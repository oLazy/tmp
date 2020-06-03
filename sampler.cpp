//
// Created by Eric Mandolesi on 03/06/2020.
//
#include "objects.h"
#include "global.h"
#include <string>
using namespace mtobj;
int main(int argn, char* argv[]) {
    if (argn > 2) {
        std::cout << "Error! this program requires a single argument.\n";
        return 1;
    }
    std::string base_filename{argv[1]};
    // init
    initPrior({0.,150000.}, {-5,2}, {-3,0}, {-90,90});
    initProposal({0.,150000.}, {-5,2}, {-3,0}, {-90,90},100.);
    Cov0 cov;
    Dataset d;
    std::ifstream is(base_filename+"_rep.dat");
    boost::archive::text_iarchive ia(is);
    ia >> d >> cov;
    is.close();
    // init model zero. it is a conductive crust lying over a resistive basement
    model m;
    m.nodes.push_back({0,-1});
    m.nodes.push_back({35000,-3});
    if(!m.isInPrior()){
        std::cerr << "model not in prior\n";
        return 17;
    }
    m.calc_params();
    // print dataset statistics
    std::cout <<
              "ELog(L): " << mtobj::expectedLogL(d, cov) <<
              "\nLog[L(m0)]: " << mtobj::logL(m, d, cov) <<
              "\nstd(ELogL): " << mtobj::stdLogL(d) <<
              "\nvar(ELogL): " << mtobj::varLogL(d) << "\n";

    for (auto itern=0; itern<300000;itern++) {
        for (int n = 0; n < m.nodes.size(); n++) {
            for (int pt = paramType::begin; pt != paramType::end; pt++) {
                if (m.nodes[n].params[pt].isActive()) {
                    auto m1 = perturb(m, n, static_cast<paramType>(pt));
                    if (m1.isInPrior()) {
                        m1.calc_params();
                        auto u = urn(gen);
                        auto l0 = logL(m, d, cov);
                        auto l1 = logL(m1, d, cov);
                        if (u < exp(l1 - l0)) {
                            m = m1;
                        }
                    }
                }
            }
        }
    }
    // print dataset statistics
    std::cout << "\n===================================\n" << m << "\n===================================\n";
    std::cout <<
              "ELog(L): " << mtobj::expectedLogL(d, cov) <<
              "\nLog[L(m*)]: " << mtobj::logL(m, d, cov) <<
              "\nstd(ELogL): " << mtobj::stdLogL(d) <<
              "\nvar(ELogL): " << mtobj::varLogL(d) << "\n";
    return 0;
}