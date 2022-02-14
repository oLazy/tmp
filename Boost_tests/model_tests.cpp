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
            auto var = it.second.maxAbsImpedance() * 0.03;
            cov0[T] = var;
            cov1[T] = {{var,0},{var,0},{var,0},{var,0}};
        }
        BOOST_CHECK_EQUAL(mtobj::logL(m,dataset,cov0),mtobj::logL(m,dataset,cov1));
    }
    BOOST_AUTO_TEST_CASE( init_cov1_from_cov0 )
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
        for (auto it : dataset){
            auto T = it.first;
            auto var = it.second.maxAbsImpedance() * 0.03;
            cov0[T] = var;
        }
        mtobj::Cov1 cov1{mtobj::initFrom(cov0)};
        BOOST_CHECK_EQUAL(mtobj::logL(m,dataset,cov0),mtobj::logL(m,dataset,cov1));
        BOOST_CHECK_EQUAL(mtobj::expectedLogL(dataset,cov0),mtobj::expectedLogL(dataset,cov1));
    }

    BOOST_AUTO_TEST_CASE( calc_kappa ) {
        //Generate a model
        mtobj::model m;
        m.nodes.push_back({0, -1});
        m.nodes.push_back({2000, -3, -2, 45.});
        m.nodes.push_back({6000, -1, -3, -15.});
        m.sort_nodes();
        m.calc_params();
        BOOST_CHECK_EQUAL(m.bin2dec(), 3);
        BOOST_CHECK_EQUAL(m.calcK(), 10);
        m.nodes[0] = {0, -1, -1, 0.};
        m.sort_nodes();
        m.calc_params();
        BOOST_CHECK_EQUAL(m.bin2dec(), 7);
        BOOST_CHECK_EQUAL(m.calcK(), 14);
        mtobj::model m1;
        m1.nodes.push_back({0, -1});
        m1.nodes.push_back({2000, -3});
        BOOST_CHECK_EQUAL(m1.bin2dec(), 0);
        BOOST_CHECK_EQUAL(m1.calcK(), 3);
        m.nodes.push_back({1, -1});
        m.nodes.push_back({2, -3, -2, 45.});
        m.nodes.push_back({6, -1, -3, -15.});
        m.nodes.push_back({10, -1});
        m.nodes.push_back({20, -3, -2, 45.});
        m.nodes.push_back({60, -1, -3, -15.});
        m.sort_nodes();
        m.calc_params();
        BOOST_CHECK_EQUAL(m.bin2dec(), 367);
        BOOST_CHECK_EQUAL(m.calcK(), 511+367);
        mtobj::model m2;
        m2.nodes.push_back({0, -1});
        BOOST_CHECK_EQUAL(m2.bin2dec(), 0);
        BOOST_CHECK_EQUAL(m2.calcK(), 1);
    }
BOOST_AUTO_TEST_SUITE_END()
