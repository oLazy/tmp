#define BOOST_TEST_MODULE MyTest
#include <boost/test/unit_test.hpp>
#include "../objects.h"
#include "../global.h"

BOOST_AUTO_TEST_SUITE(test_mt_model)
BOOST_AUTO_TEST_CASE( cov0_vs_cov1 )
{
    BOOST_CHECK_EQUAL(1,1);
}
    BOOST_AUTO_TEST_CASE( cov1_vs_cov1 )
    {
        BOOST_CHECK_EQUAL(1,1);
    }
BOOST_AUTO_TEST_SUITE_END()
