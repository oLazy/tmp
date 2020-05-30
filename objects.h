//
// Created by Eric Mandolesi on 30/05/2020.
//

#ifndef MT1DANISMODELPARAMS_OBJECTS_H
#define MT1DANISMODELPARAMS_OBJECTS_H
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <map>
#include <boost/serialization/map.hpp>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/random.hpp>
#include <fstream>
#include <cmath>
#include <numeric>
#include <exception>

extern int seed;
extern boost::random::mt19937 gen;
extern boost::random::normal_distribution<double> rn;
namespace mtobj {
    typedef std::pair<double, double> limits;
    enum paramType {
        begin = 0, depth = begin, sigmaMean, sigmaRatio, beta, end
    };
    std::map<int, std::string> param_map{{0, "depth"},
                                         {1, "mean"},
                                         {2, "ratio"},
                                         {3, "beta"}};
    std::map<int, std::pair<double, double> > prior;

    void initPrior(limits depth, limits mean, limits ratio, limits beta) {
        prior[paramType::depth] = depth;
        prior[paramType::sigmaMean] = mean;
        prior[paramType::sigmaRatio] = ratio;
        prior[paramType::beta] = beta;
    }

    std::map<int, double> proposal; //stores the sd for perturbations
    void initProposal(limits depth, limits mean, limits ratio, limits beta, double scale) {
        proposal[paramType::depth] = (depth.second - depth.first) / scale;
        proposal[paramType::sigmaMean] = (mean.second - mean.first) / scale;
        proposal[paramType::sigmaRatio] = (ratio.second - ratio.first) / scale;
        proposal[paramType::beta] = (beta.second - beta.first) / scale;
    }

    class Parameter {
        paramType type;
        bool active{true};
        double value;
        friend boost::serialization::access;

        template<typename Archive>
        void serialize(Archive &ar, const unsigned int version) {
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

    class ModelError : std::logic_error {
    public:
        ModelError(const std::string &string) : logic_error(string) { message = string; }

    private:
        std::string message;

        const char *what() const throw() {
            char *str;
            strcpy(static_cast<char *>(str), "ModelError: ");
            strcat(static_cast<char *>(str), message.c_str());
            return (str);
        }
    };

    struct node {
        std::map<int, Parameter> params;
    public:
        node() = default;

        node(double depth, double sigmaMean) {
            if (depth > 0) {
                params[paramType::depth] = Parameter(paramType::depth,
                                                     true,
                                                     depth);
            } else if (depth == 0) {
                params[paramType::depth] = Parameter(paramType::depth,
                                                     false, // ground node cannot be moved
                                                     depth);
            } else {
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

        node(double depth, double sigmaMean, double sigmaRatio, double beta) {
            if (depth > 0) {
                params[paramType::depth] = Parameter(paramType::depth,
                                                     true,
                                                     depth);
            } else if (depth == 0) {
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
        void serialize(Archive &ar, const unsigned int version) {
            ar & params;
        }

    };

    struct model {
        model() = default;

        model(const model &m) {
            for (const auto &n:m.nodes) {
                this->nodes.push_back(n);
            }
        } // implement copy constructor, it must copy the nodes but not other fields.
        model &operator=(const model &m) {
            model res;
            for (const auto &n:m.nodes) {
                this->nodes.push_back(n);
            }
            return *this;
        } // implement copy assignment, it must copy the nodes but not other fields.
        std::vector<node> nodes;
        std::vector<double> _h, _al, _at, _blt;

        node operator[](int i) const { return nodes[i]; }

        node &operator[](int i) { return nodes[i]; }

        node getNode(double depth) {
            for (int i = 0; i < nodes.size() - 1; i++) {
                if (nodes[i + 1].params[paramType::depth].getValue() > depth &&
                    nodes[i].params[paramType::depth].getValue() < depth) {
                    return nodes[i];
                }
            }
            throw ModelError("no value available at specified depth.");
        }

        bool isValid() {
            std::vector<double> z;
            for (auto n : nodes) {
                z.push_back(n.params[paramType::depth].getValue());
            }
            std::vector<double> diff(z.size());
            std::adjacent_difference(z.begin(), z.end(), diff.begin());
            std::cout << "diff elements:\n";
            for (auto d:diff) {
                std::cout << d << "\n";
            }
            if (std::any_of(diff.begin()++, diff.end(), [](double x) { return x < 0; })) {
                for (auto d = diff.begin()++; d != diff.end(); d++) {
                    std::cout << "if any=true loop\n";
                    std::cout << *d << "\n";
                }
                return false;
            }
            return true;
        }

        bool isInPrior() {

            auto params_in_prior = [](node x) {
                for (int p = paramType::begin; p != paramType::end; p++) {
                    if (x.params[p].isActive()) {
                        if (x.params[p].getValue() < prior[p].first || x.params[p].getValue() > prior[p].second) {
                            return false;
                        }
                    } // only check for active parameter, the inactive parameters are set to nan.
                }
                return true;
            };
            return std::all_of(nodes.begin(), nodes.end(), params_in_prior);
        }

        void calc_params() {
            for (int i = 0; i < nodes.size(); i++) {
                if (i != nodes.size() - 1) {
                    auto h0 = nodes[i].params[paramType::depth].getValue();
                    auto h1 = nodes[i + 1].params[paramType::depth].getValue();
                    auto thick = h1 - h0;
                    _h.push_back(thick);
                } else {
                    _h.push_back(0);
                }
                auto sm = pow(10., nodes[i].params[paramType::sigmaMean].getValue()); // 10^param, sampling in log space
                if (nodes[i].params[paramType::sigmaRatio].isActive()) {
                    auto sr = pow(10.,
                                  nodes[i].params[paramType::sigmaRatio].getValue()); // 10^param, sampling in log space
                    auto beta = nodes[i].params[paramType::beta].getValue() * M_PI /
                                180.; // angle stored in deg, used in rads
                    _al.push_back(2 * sm * sr / (1. + sr));
                    _at.push_back(2 * sm / (1. + sr));
                    _blt.push_back(beta);
                } else {
                    _al.push_back(sm);
                    _at.push_back(sm);
                    _blt.push_back(0);
                }
            }
        }

        friend std::ostream &operator<<(std::ostream &os, const model &model) {
            int node_id = 0;
            for (auto n : model.nodes) {
                for (int pt = paramType::begin; pt != paramType::end; pt++) {
                    os << "[" << node_id << "]" << param_map[pt] << ":" << n.params[pt] << "\n";
                }
                node_id++;
            }
            return os;
        }
    };

    model perturb(model const &m0, int node_id, paramType pt) {
        model mp = m0;
        if (mp[node_id].params[pt].isActive()) {
            auto ran = rn(gen);
            std::cerr << "ran = " << ran << "\n"; // log the random number
            auto p_i = mp[node_id].params[pt].getValue() + ran * proposal[pt];
            mp[node_id].params[pt].setValue(p_i);
        }
        return mp;
    }
}
#endif //MT1DANISMODELPARAMS_OBJECTS_H
