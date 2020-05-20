#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <map>
#include <boost/serialization/map.hpp>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <fstream>
#include <cmath>
#include <numeric>
#include <exception>

typedef std::pair<double, double> limits ;
enum paramType{begin=0, depth=begin, sigmaMean, sigmaRatio, beta, end};
std::map<int,std::pair<double,double> > prior;
void initPrior(limits depth, limits mean, limits ratio, limits beta){
    prior[paramType::depth] = depth;
    prior[paramType::sigmaMean] = mean;
    prior[paramType::sigmaRatio] = ratio;
    prior[paramType::beta] = beta;
}
class Parameter{
    paramType type;
    bool active{true};
    double value;
    friend boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version){
        ar & type;
        ar & active;
        ar & value;
    }
public:
    Parameter() = default;
    Parameter(paramType type, bool isActive, double value) : type(type), active(isActive), value(value) {}

    paramType getType() const {
        return type;
    }

    bool isActive() const {
        return active;
    }

    double getValue() const {
        return value;
    }

    void setType(paramType type) {
        Parameter::type = type;
    }

    void setIsActive(bool isActive) {
        Parameter::active = isActive;
    }

    void setValue(double value) {
        Parameter::value = value;
    }

    friend std::ostream &operator<<(std::ostream &os, const Parameter &parameter) {
        os << "type: " << parameter.type << " active: " << parameter.active << " value: " << parameter.value;
        return os;
    }
};

class ModelError : std::logic_error{
public:
    ModelError(const std::string &string) : logic_error(string) {message=string;}

private:
    std::string message;

    const char * what () const throw ()
    {char *str;
        strcpy(static_cast<char *>(str), "ModelError: ");
        strcat(static_cast<char *>(str), message.c_str());
        return (str);
    }
};
struct node{
    std::map<int,Parameter> params;
public:
    node() = default;
    node(double depth, double sigmaMean){
        if(depth > 0){
            params[paramType::depth] = Parameter(paramType::depth,
                                                 true,
                                                 depth);
        }else if (depth==0){
            params[paramType::depth] = Parameter(paramType::depth,
                                                 false, // ground node cannot be moved
                                                 depth);
        }else{
            throw std::logic_error("negative depth in node construction (iso).");
        }
        params[paramType::sigmaMean] = Parameter(paramType::sigmaMean,
                                                 true,
                                                 sigmaMean);
        params[paramType::sigmaRatio] = Parameter(paramType::sigmaRatio,
                                                  false,
                                                  nan(" "));
        params[paramType::beta] = Parameter(paramType::beta,
                                            false,
                                            nan(" "));
    }
    node(double depth, double sigmaMean, double sigmaRatio, double beta){
        if (depth > 0) {
            params[paramType::depth] = Parameter(paramType::depth,
                                                 true,
                                                 depth);
        } else if (depth == 0){
            params[paramType::depth] = Parameter(paramType::depth,
                                                 false,
                                                 depth);
        } else {
            throw std::logic_error("negative depth in node construction (anis).");
        }
        params[paramType::sigmaMean] = Parameter(paramType::sigmaMean,
                                                 true,
                                                 sigmaMean);
        params[paramType::sigmaRatio] = Parameter(paramType::sigmaRatio,
                                                  true,
                                                  sigmaRatio);
        params[paramType::beta] = Parameter(paramType::beta,
                                            true,
                                            beta);
    }

//    friend std::ostream &operator<<(std::ostream &os, const node &node) {
//        os << "params: \n";
//        for (int i = paramType::begin; i != paramType::end; i++){
//            os << params[i] << "\n";
//        }
//        return os;
//    }

private:
    friend boost::serialization::access;
    template<typename Archive>
    void serialize(Archive &ar, const unsigned int version){
        ar & params;
    }

};

struct model{
    std::vector<node> nodes;
    node getNode(int i) const {return nodes[i];}
    node getNode(double depth) {
        for(int i=1; i<nodes.size();i++){
            if(nodes[i].params[paramType::depth].getValue()>depth &&
               nodes[i-1].params[paramType::depth].getValue()<depth){
                return nodes[i];
            }
        }
        throw ModelError("no value available at specified depth.");
    }
    bool isValid(){
        std::vector<double> z;
        for(auto n : nodes){
            z.push_back(n.params[paramType::depth].getValue());
        }
        std::vector<double> diff(z.size());
        std::adjacent_difference(z.begin(),z.end(),diff.begin());
        std::cout << "diff elements:\n";
        for (auto d:diff){
            std::cout << d << "\n";
        }
        if(std::any_of(diff.begin()++,diff.end(),[](double x){return x<0;})){
            for (auto d=diff.begin()++;d!=diff.end();d++) {
                std::cout << "ifany=true loop\n";
                std::cout << *d << "\n";
            }
            return false;
        }
        return true;
    }

    bool isInPrior() {

        auto params_in_prior = [](node x){
            for (int p = paramType::begin; p != paramType::end; p++) {
                if (x.params[p].getValue() < prior[p].first || x.params[p].getValue() > prior[p].second) {
                    return false;
                }
                return true;
            }
        };
        return std::all_of(nodes.begin(),nodes.end(),params_in_prior);
    }

};

int main() {
    initPrior({0.,410000.}, {-5,2}, {0,-3}, {-90,90});
    std::cout << "Hello, World!" << std::endl;
//    Parameter p(paramType::sigmaMean, true, 3.1);
//    std::ofstream os("aTestNode.txt");

//    node n(0,0.7);

    std::ifstream is("aTestNode.txt");
//    boost::archive::text_oarchive oa(os);
    boost::archive::text_iarchive ia(is);

//    oa << n;
    node q;
    ia >> q;
//    std::cout << q << std::endl;
//    if (q.getType()==paramType::sigmaMean){
//        std::cout << "evviva!\n";
//    }

    node p{q};


    // modify all active parameters
    for (int type = paramType::begin; type!=paramType::end; type++) {
        if (p.params[type].isActive()) {
            p.params[type].setValue(p.params[type].getValue() + 1.0);
            // do metropolis hastings
        }
    }
    for (int pt=paramType::begin; pt!=paramType::end; pt++){
        std::cout << "q: " << q.params[pt] <<"\n";
        std::cout << "p: " << p.params[pt] <<"\n";
    }
    model m;
    m.nodes.push_back(node(0,1));
    m.nodes.push_back(node(1,2));
    m.nodes.push_back(node(3,3));
    m.nodes.push_back(node(6,4));

    std::cout << m.isValid() << "\n";
//    oa << p;
    return 0;
}
