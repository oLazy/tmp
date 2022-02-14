#define BOOST_TEST_MODULE MyTest
#include <boost/test/unit_test.hpp>
#include "../objects.h"
#include "../global.h"

BOOST_AUTO_TEST_SUITE(test_mt_model)
BOOST_AUTO_TEST_CASE( cov0_vs_cov1 )
{
    //Generate a model
    mtobj::model m;
    m.nodes.push_back({0,-1});
    m.nodes.push_back({2000,-3,-2,45.});
    m.nodes.push_back({6000,-1,-3,-15.});
    m.calc_params();
    //Generate Dataset
    mtobj::Dataset dataset;
    for (auto i=0; i<31; i++){
        double increment = 5./30.;
        double per = pow(10,-2+static_cast<double>(i)*increment);
        dataset[per] = m(per);
    }
    //Generate two equivalent covariances
    mtobj::Cov0 cov0;
    mtobj::Cov1 cov1;
    for (auto it : dataset){
        auto T = it.first;
        
    }
}
BOOST_AUTO_TEST_SUITE_END()
