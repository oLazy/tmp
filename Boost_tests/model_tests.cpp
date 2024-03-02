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
    BOOST_AUTO_TEST_CASE( MTTensor_determinant ) {
        { // small local test
            //Generate a tensor
            MTTensor z{{1,  -2},
                       {3,  -4},
                       {-5, 2},
                       {-1, 2}};
            MTComplex det_z_pek = z.xx * z.yy - z.xy * z.yx;
            MTComplex det_z_eric = z.det();
            BOOST_CHECK_EQUAL(det_z_eric, det_z_pek);
        }
        { // improved coverage test
            mtobj::model m;
            m.nodes.push_back({0, -1});
            m.nodes.push_back({2000, -3, -2, 45.});
            m.nodes.push_back({6000, -1, -3, -15.});
            m.calc_params();
            //Generate and check Dataset
            mtobj::Dataset dataset;
            for (auto i = 0; i < 31; i++) {
                double increment = 5. / 30.;
                double per = pow(10, -2 + static_cast<double>(i) * increment);
                auto const z = m(per);
                dataset[per] = z;
                MTComplex det_z_pek = z.xx * z.yy - z.xy * z.yx;
                MTComplex det_z_eric = z.det();
                BOOST_CHECK_EQUAL(det_z_eric, det_z_pek);
            }
            gp_utils::dataset2disk(dataset,std::cout, true);
        }
    }
    BOOST_AUTO_TEST_CASE( test_mt_model_fw ){
        mtobj::model m;
        m.nodes.push_back({0, -1});
//            m.nodes.push_back({2000, -3, -2, 45.});
//            m.nodes.push_back({6000, -1, -3, -15.});
        m.calc_params();
        //Generate and check Dataset
        double per;
        mtobj::Dataset dataset;
        for (auto i = 0; i < 31; i++) {
            double increment = 5. / 30.;
            per = pow(10, -2 + static_cast<double>(i) * increment);
            double xx = std::real(m(per).xx);
            double mxx = std::real(m.compute_impedance_at_the_basement_top(per, &(*m.nodes.crbegin())).xx);
            double xy = std::real(m(per).xy);
            double mxy = std::real(m.compute_impedance_at_the_basement_top(per, &(*m.nodes.crbegin())).xy);
            double yx = std::real(m(per).yx);
            double myx = std::real(m.compute_impedance_at_the_basement_top(per, &(*m.nodes.crbegin())).yx);
            double yy = std::real(m(per).yy);
            double myy = std::real(m.compute_impedance_at_the_basement_top(per, &(*m.nodes.crbegin())).yy);
            BOOST_CHECK_EQUAL(xx,mxx);
            BOOST_CHECK_EQUAL(xy,mxy);
            BOOST_CHECK_EQUAL(yx,myx);
            BOOST_CHECK_EQUAL(yy,myy);
        }
        m.nodes.push_back({2000, -3, -2, 45.});
        m.calc_params();
        double xy = std::real(m(per).xy);
        double mxy = std::real(m.myOperator(per).xy);
        BOOST_CHECK_EQUAL(xy,mxy);
    }
    BOOST_AUTO_TEST_CASE( test_print_element ){
        gp_utils::print_element("pippo", 15, std::cout);
        gp_utils::print_element("pluto", 15, std::cout);
        gp_utils::print_element("topolino", 15, std::cout);
        std::cout << std::endl;
        gp_utils::print_element("pluto", 15, std::cout);
        gp_utils::print_element("topolino", 15, std::cout);
        gp_utils::print_element("pippo", 15, std::cout);
        std::cout << std::endl;
        gp_utils::print_element("topolino", 15, std::cout);
        gp_utils::print_element("pluto", 15, std::cout);
        gp_utils::print_element("pippo", 15, std::cout);

        std::cout << std::setprecision(2) << std::scientific;
        gp_utils::print_element(22.2, 15, std::cout);
        gp_utils::print_element(2.22, 15, std::cout);
        gp_utils::print_element(2.938567983475e-34, 15, std::cout);
        gp_utils::print_element(0.000000000042525324532, 15, std::cout);
        std::cout << std::endl;

    }

BOOST_AUTO_TEST_SUITE_END()
