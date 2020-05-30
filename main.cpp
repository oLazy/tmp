#include "objects.h"
#include "global.h"
// declare global variables
using namespace mtobj;

int main() {
    MTTensor z{{1,1},{2,2},{3,3},{4,4}};
    std::cout << z << std::endl;
    // initialize prior and proposal sd
    initPrior({0.,410000.}, {-5,2}, {-3,0}, {-90,90});
    initProposal({0.,410000.}, {-5,2}, {-3,0}, {-90,90},100.);
    std::cout << proposal[0] << std::endl;
    //    Parameter p(paramType::sigmaMean, true, 3.1);
//    std::ofstream os("aTestNode.txt");

//    node n(0,0.7);
//
//    std::ifstream is("aTestNode.txt");
////    boost::archive::text_oarchive oa(os);
//    boost::archive::text_iarchive ia(is);
//
////    oa << n;
//    node q;
//    ia >> q;
//    std::cout << q << std::endl;
//    if (q.getType()==paramType::sigmaMean){
//        std::cout << "evviva!\n";
//    }

//    node p{q};

//
//    // modify all active parameters
//    for (int type = paramType::begin; type!=paramType::end; type++) {
//        if (p.params[type].isActive()) {
//            p.params[type].setValue(p.params[type].getValue() + 1.0);
//            // do metropolis hastings
//        }
//    }


//    for (int pt=paramType::begin; pt!=paramType::end; pt++){
//        std::cout << "q: " << q.params[pt] <<"\n";
//        std::cout << "p: " << p.params[pt] <<"\n";
//    }
    model m;
    m.nodes.push_back(node(0,1));
    m.nodes.push_back(node(100,-2));
    m.nodes.push_back(node(300,-3, -2, -45.2));
    m.nodes.push_back(node(6000,-4));

    if (m.isInPrior()) m.calc_params();
    std::cout << m(0.77);

    int node_id = 1; // the node I am going to perturb
    model m2 = perturb(m,node_id,paramType::depth);
    m2.calc_params();
    std::cout << m << std::endl;
    std::cout << std::endl;
    std::cout << m2 << std::endl;

//    std::cout << m.isValid() << "\n";
//    oa << p;
    return 0;
}
