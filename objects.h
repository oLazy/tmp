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
#include <algorithm>
#include <boost/histogram.hpp>
#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#ifdef _OMP
#include <omp.h>
#endif


#include "MTTensor.h"


extern int seed;
extern boost::random::mt19937 gen;
extern boost::random::normal_distribution<double> rn;
extern boost::random::uniform_real_distribution<double> urn;

namespace mtobj {
    /*include code from the old project */
    //TODO adjust the calc_params method to ensure physical parameter units are correct

    typedef std::complex<double> MTComplex;
    typedef std::complex<double> dcomp;
    typedef std::pair<double, double> limits;

    typedef std::map<double, double> Cov0;
    typedef std::map<double, MTTensor> Dataset;
    MTComplex static constexpr ic{0, 1};
    double static constexpr pi{M_PI};

    class ModelOutOfPrior : public std::exception {
        virtual const char *what() const throw() {
            return "Model out prior";
        }
    };


    MTTensor rot_z(MTTensor const &za, double beta_rad) {
        if (std::isnan(beta_rad)) beta_rad = 0.;
        MTTensor result;
        double co2 = cos(2. * beta_rad);
        double si2 = sin(2. * beta_rad);

        dcomp sum1 = za.xx + za.yy;
        dcomp sum2 = za.xy + za.yx;

        dcomp dif1 = za.xx - za.yy;
        dcomp dif2 = za.xy - za.yx;

        result.xx = 0.5 * (sum1 + dif1 * co2 + sum2 * si2);
        result.xy = 0.5 * (dif2 + sum2 * co2 - dif1 * si2);
        result.yx = 0.5 * (-dif2 + sum2 * co2 - dif1 * si2);
        result.yy = 0.5 * (sum1 - dif1 * co2 - sum2 * si2);
        return result;
    }

    [[nodiscard]] inline dcomp dfp(dcomp const &x) noexcept {
        return 1.0 + std::exp(-2.0 * x);
    }

    [[nodiscard]] inline dcomp dfm(dcomp const &x) noexcept {
        return 1.0 - std::exp(-2.0 * x);
    }

    enum paramType {
        begin = 0, depth = begin, sigmaMean, sigmaRatio, beta, end
    };
    std::map<int, std::string> param_map{{0, "depth"},
                                         {1, "mean"},
                                         {2, "ratio"},
                                         {3, "beta"}};
    typedef std::map<int, std::pair<double, double> > Prior;
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

    bool isLogLhoodExpected(double logL, double ElogL, double logLStd, double beta, int nStd=3){
        auto semi_interval = (double)nStd*pow(beta,-1)*logLStd;
        if( ElogL-semi_interval <= logL)return true; // only check for lower bound
        return false;
    }

    std::tuple<double, double, double, double> balanceReceiptWeights(boost::program_options::variables_map& vm_p){
        auto pp = vm_p["perturb"].as<double>();
        auto pb = vm_p["birth"].as<double>();
        auto pd = vm_p["death"].as<double>();
        auto ps = vm_p["iso-switch"].as<double>();
        auto sum{pp+pb+pd+ps};
        return std::make_tuple(pp/sum,pb/sum,pd/sum,ps/sum);
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

        void setType(paramType _type) {
            Parameter::type = _type;
        }

        void setIsActive(bool isActive) {
            Parameter::active = isActive;
        }

        void setValue(double _value) {
            Parameter::value = _value;
        }

        friend std::ostream &operator<<(std::ostream &os, const Parameter &parameter) {
            os << "type: " << parameter.type << " active: " << parameter.active << " value: " << parameter.value;
            return os;
        }
    };

    class ModelError : public std::logic_error {
    public:
        explicit ModelError(const std::string &string) : logic_error(string) { message = string; }

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

    private:
        friend boost::serialization::access;

        template<typename Archive>
        void serialize(Archive &ar, const unsigned int version) {
            ar & params;
        }

    };

    node generateIsoNode(Prior const &p) {

        double const mind = p.at(paramType::depth).first;
        double const maxd = p.at(paramType::depth).second;
        double depth = boost::random::uniform_real_distribution<double>(mind, maxd)(gen);

        double const mins = p.at(paramType::sigmaMean).first;
        double const maxs = p.at(paramType::sigmaMean).second;
        double sm = boost::random::uniform_real_distribution<double>(mins, maxs)(gen);
        node res(depth, sm);
        return res;
    }

    node generateAnisNode(Prior const &p) {
        node res;
        for (int pt = paramType::begin; pt != paramType::end; pt++) {
            double const min = p.at(pt).first;
            double const max = p.at(pt).second;
            double value = boost::random::uniform_real_distribution<double>(min, max)(gen);
            res.params[pt].setType(static_cast<paramType>(pt));
            res.params[pt].setIsActive(true);
            res.params[pt].setValue(value);
        }
        return res;
    }

    struct model {
        model() = default;

        model(const model &m) {
            for (const auto &n:m.nodes) {
                this->nodes.push_back(n);
            }
            this->logL = m.logL;
            this->beta = m.beta;
        } // implement copy constructor, it must copy the nodes but not other fields.

        model &operator=(const model &m) = default;

        // fields ::
        std::vector<node> nodes;
        std::vector<double> _h, _al, _at, _blt;
        double logL{0};
        double beta{1.};

        node operator[](int i) const { return nodes[i]; }

        node &operator[](int i) { return nodes[i]; }

        node getNode(double depth) {
            auto nnode = nodes.size() - 1;
            auto deep_z = nodes[nnode].params[paramType::depth].getValue();
            for (int i = 0; i < nnode; i++) {
                auto z1 = nodes[i + 1].params[paramType::depth].getValue();
                auto z0 = nodes[i].params[paramType::depth].getValue();
                if ((z1 > depth) && (z0 <= depth)) {
                    return nodes[i];
                }
            }
            return nodes[nnode];

            throw ModelError("no value available at specified depth.");
        }

        void setLogL(double logL) {
            model::logL = logL;
        }

        void setBeta(double beta) {
            model::beta = beta;
        }

        bool isValid() {
            std::vector<double> z;
            for (auto n : nodes) {
                z.push_back(n.params[paramType::depth].getValue());
            }
            std::vector<double> diff(z.size());
            std::adjacent_difference(z.begin(), z.end(), diff.begin());

            if (std::any_of(diff.begin()++, diff.end(), [](double x) { return x < 0; })) {
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

        void sort_nodes() {
            auto comp = [](const node &lhs, const node &rhs) {
                return lhs.params.at(paramType::depth).getValue() < rhs.params.at(paramType::depth).getValue();
            };
            std::sort(nodes.begin(), nodes.end(), comp);
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
                    _at.push_back(2 * sm * sr / (1. + sr));
                    _al.push_back(2 * sm / (1. + sr));
                    _blt.push_back(beta);
                } else {
                    _al.push_back(sm);
                    _at.push_back(sm);
                    _blt.push_back(0);
                }
            }
        }

        MTTensor operator()(const double &x) const noexcept { // pek algorithm, my implementation. TESTED

            dcomp k0{(1.0 - ic) * 2. * pi * pow(10., -3.) / sqrt(10. * x)};
            // compute the impedance on the top of the homogeneous
            // basement in the direction of its strike
            auto number_of_layers = nodes.size() - 1;
            auto number_of_interfaces = nodes.size();

            auto i_layer = static_cast<unsigned long>(number_of_layers);

            double a1 = _al[i_layer];
            double a2 = _at[i_layer];
            double bs = _blt[i_layer];
            double a1is = 1. / sqrt(a1);
            double a2is = 1. / sqrt(a2);
            MTTensor z{{0, 0},
                       {0, 0},
                       {0, 0},
                       {0, 0}};
            MTTensor z_rot;
            MTTensor z_bot;
            z_rot.xx = 0.;
            z_rot.xy = k0 * a1is;
            z_rot.yx = -k0 * a2is;
            z_rot.yy = 0.;
            /*
             c> If no more layers are present in the model, rotate the
             c> impedance into the original coordinate system and return
             */
            if (number_of_interfaces == 1) {

                z = rot_z(z_rot, -bs);
                return z;
            }

            double bs_ref = bs;


            for (unsigned long layer = (number_of_layers - 1); layer != static_cast<unsigned long>(-1); --layer) {
                i_layer = static_cast<unsigned long>(layer);
                /* in the old code this was in km.
                 * now the physical parameterization is computed on a different place, so the unit
                 * transformation is not needed here */
//                double hd = 1000. * _h[i_layer];
                double hd = _h[i_layer];
                a1 = _al[i_layer];
                a2 = _at[i_layer];
                bs = _blt[i_layer];
                /*
                 c
                 c> If the strike direction differs from that of the previous
                 c> layer, rotate the impedance into the coordinate system of
                 c> the current anisotropy strike
                 c
                 */
                dcomp dt_z_bot = z_rot.xx * z_rot.yy - z_rot.xy * z_rot.yx;
                if (bs != bs_ref && a1 != a2) {
                    z_bot = rot_z(z_rot, bs - bs_ref);

                } else {
                    z_bot = z_rot;
                    bs = bs_ref;
                }

                dcomp k1 = k0 * sqrt(a1);
                dcomp k2 = k0 * sqrt(a2);
                a1is = 1. / sqrt(a1);
                a2is = 1. / sqrt(a2);
                dcomp dz1 = k0 * a1is;
                dcomp dz2 = k0 * a2is;
                dcomp ag1 = k1 * hd;
                dcomp ag2 = k2 * hd;

                /*
                 c
                 c> Propagate the impedance tensor from the bottom to the top
                 c> of the current layer
                 c
                 */

                dcomp z_denominator = dt_z_bot * dfm(ag1) * dfm(ag2) / (dz1 * dz2) +
                                      z_bot.xy * dfm(ag1) * dfp(ag2) / dz1 -
                                      z_bot.yx * dfp(ag1) * dfm(ag2) / dz2 +
                                      dfp(ag1) * dfp(ag2);

                z_rot.xx = 4. * z_bot.xx * std::exp(-ag1 - ag2) / z_denominator;
                z_rot.xy = (z_bot.xy * dfp(ag1) * dfp(ag2) -
                            z_bot.yx * dfm(ag1) * dfm(ag2) * dz1 / dz2 +
                            dt_z_bot * dfp(ag1) * dfm(ag2) / dz2 +
                            dfm(ag1) * dfp(ag2) * dz1) / z_denominator;
                z_rot.yx = (z_bot.yx * dfp(ag1) * dfp(ag2) -
                            z_bot.xy * dfm(ag1) * dfm(ag2) * dz2 / dz1 -
                            dt_z_bot * dfm(ag1) * dfp(ag2) / dz1 -
                            dfp(ag1) * dfm(ag2) * dz2) / z_denominator;
                z_rot.yy = 4. * z_bot.yy * std::exp(-ag1 - ag2) / z_denominator;

                bs_ref = bs;

            }

            if (bs_ref != 0.0) {
                z = rot_z(z_rot, -bs_ref);
            } else {
                z = z_rot;
            }
            return z;


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

    private:
        friend boost::serialization::access;

        template<typename Archive>
        void serialize(Archive &ar, const unsigned int version) {
            ar & nodes;
        }
    };

    model perturb(model const &m0, int node_id, paramType pt) {
        model mp = m0;
        if (mp[node_id].params[pt].isActive()) {
            auto ran = rn(gen);
//            std::cerr << "ran = " << ran << "\n"; // log the random number
            auto p_i = mp[node_id].params[pt].getValue() + ran * proposal[pt];
            mp[node_id].params[pt].setValue(p_i);
        }
        return mp;
    }

    model death(model const &m0) {
        if (m0.nodes.size() < 2) {
            throw ModelOutOfPrior();
        }
        boost::random::uniform_int_distribution<int> r_node(1, m0.nodes.size() - 1);
        auto nn = r_node(gen);
        model res;
        int i{0};
        for (auto n:m0.nodes) {
            if (i != nn) res.nodes.push_back(n); // copy only "live" nodes
            i++;
        }
        res.setBeta(m0.beta);
        return res;
    }

    enum class birthType {
        iso, anis, any
    };

    model birth(model const &m0, birthType t) {
        node newNode;
        switch (t) {
            case birthType::iso: {
                newNode = generateIsoNode(prior);
            }
                break;
            case birthType::anis: {
                newNode = generateAnisNode(prior);
            }
                break;
            case birthType::any: {
                double rand_num = urn(gen);
                if (rand_num <= 0.5) {
                    newNode = generateIsoNode(prior);
                } else {
                    newNode = generateAnisNode(prior);
                }
            }
        }
        model m1 = m0;
        m1.nodes.push_back(newNode);
        auto comp = [](node const &a, node const &b) {
            return (
                    a.params.at(paramType::depth).getValue() <
                    b.params.at(paramType::depth).getValue()
            );
        };
//        std::sort(m1.nodes.begin(), m1.nodes.end(), comp); // TODO TEST THIS ACTUALLY RETURNS m1 SORTED PROPERLY
        m1.sort_nodes();
        return m1;
    }

    model iso_switch(model const &m0, int node_id) {
        model mp = m0;
        if (mp[node_id].params[paramType::sigmaRatio].isActive()) {
            mp[node_id].params[paramType::sigmaRatio].setValue(nan(""));
            mp[node_id].params[paramType::beta].setValue(nan(""));
            mp[node_id].params[paramType::sigmaRatio].setIsActive(false);
            mp[node_id].params[paramType::beta].setIsActive(false);
        } else {
            // generate ratio
            double const min = prior.at(paramType::sigmaRatio).first;
            double const max = prior.at(paramType::sigmaRatio).second;
            double value = boost::random::uniform_real_distribution<double>(min, max)(gen);
            mp[node_id].params[paramType::sigmaRatio].setIsActive(true);
            mp[node_id].params[paramType::sigmaRatio].setValue(value);
            // generate beta
            double const minb = prior.at(paramType::beta).first;
            double const maxb = prior.at(paramType::beta).second;
            value = boost::random::uniform_real_distribution<double>(minb, maxb)(gen);
            mp[node_id].params[paramType::beta].setIsActive(true);
            mp[node_id].params[paramType::beta].setValue(value);

        }
        return mp;
    }

    void parallel_tempering_swap(std::vector<model> &chains){
        auto n_chains = chains.size();
        boost::random::uniform_int_distribution<int> chain_picker(0,n_chains-1);
        for (auto i = 0; i < n_chains * n_chains; i++){ // propose n^2 switches
            int j = chain_picker(gen);
            int k = chain_picker(gen);
            if(j!=k){
                auto bj = chains[j].beta;
                auto bk = chains[k].beta;
                if(urn(gen) <  pow(exp(chains[j].logL - chains[k].logL),(bk-bj))){
                    std::cerr << "swapping chain no. " << j << " with chain no. " << k <<"\n";
                    chains[j].setBeta(bk);
                    chains[k].setBeta(bj);
                    std::swap(chains[j],chains[k]); // i-th temperature remains associated to the i-th chain
                }
            }
        }
    }


//    model isoJump(model const &m0, int node_id){
//
//    }
/*** ASSUMPTIONS:
 * the error statistics is normal, all the impedances are affected by the same error
 * which is dominated by the noise in the electric channel.
 * @param m : model
 * @param d : data
 * @param cov : covariance
 * @return : log(likelihood)
 */
    double logL(model const &m,
                std::map<double, MTTensor> const &d,
                std::map<double, double> const &cov) {

        static const double f{8}; // 8=4 real and 4 imag parts
        auto N = d.size() * f;
        auto c = -N * 0.5 * log(2 * M_PI);
        double sum_log_sigma{0};
        double sum_res{0};
#ifdef _OMP
        #pragma omp parallel for reduction (+:sum_log_sigma, sum_res)
        for (int i = 0; i < d.size(); i++) {
            auto it = d.begin();
            std::advance(it, i);
            auto dat = *it;
#else
        for (auto dat: d){
#endif
            double const T = dat.first;
            MTTensor const z_meas = dat.second;
            MTTensor z_pred = m(T);
            auto const sigma = cov.at(T);
            sum_log_sigma += f * log(sigma);
            sum_res += pow((std::real(z_meas.xx) - std::real(z_pred.xx)) / sigma, 2);
            sum_res += pow((std::real(z_meas.xy) - std::real(z_pred.xy)) / sigma, 2);
            sum_res += pow((std::real(z_meas.yx) - std::real(z_pred.yx)) / sigma, 2);
            sum_res += pow((std::real(z_meas.yy) - std::real(z_pred.yy)) / sigma, 2);
            sum_res += pow((std::imag(z_meas.xx) - std::imag(z_pred.xx)) / sigma, 2);
            sum_res += pow((std::imag(z_meas.xy) - std::imag(z_pred.xy)) / sigma, 2);
            sum_res += pow((std::imag(z_meas.yx) - std::imag(z_pred.yx)) / sigma, 2);
            sum_res += pow((std::imag(z_meas.yy) - std::imag(z_pred.yy)) / sigma, 2);
        } // end of omp parallel region. Here sum_res and sum_sigma should be reduced
        return c - sum_log_sigma - 0.5 * sum_res;
    }


    double expectedLogL(
            std::map<double, MTTensor> const &d,
            std::map<double, double> const &cov) {
        static const double f{8}; // 8=4 real and 4 imag parts
        auto N = d.size() * f;
        auto c = -N * 0.5 * log(2 * M_PI);
        double sum_log_sigma{0};
        for (auto dat: d) {
            double const T = dat.first;
            MTTensor const z_meas = dat.second;
            auto const sigma = cov.at(T);
            sum_log_sigma += f * log(sigma);
        }
        return c - sum_log_sigma - 0.5 * N;
    }

    double varLogL(std::map<double, MTTensor> const &d) {
        static const double f{8}; // 8=4 real and 4 imag parts
        auto N = d.size() * f;
        return 2 * N;
    }

    double stdLogL(std::map<double, MTTensor> const &d) {
        static const double f{8}; // 8=4 real and 4 imag parts
        auto N = d.size() * f;
        return sqrt(2 * N);
    }

    void calc_beta(int n_temp, double max_temp, model &m, std::vector<model> &chains, bool printSchedule=false) { // make tempering schedule
        double length = log10(max_temp);
        double delta = length / static_cast<double>(n_temp - 1);
        for (int i = -1; i < n_temp; i++) {
            double temperature{0};
            if (i == -1) {
                temperature = 1.;
            } else {
                temperature = pow(10, static_cast<double>(i) * delta);
            }
            double beta = pow(temperature, -1);
            if(printSchedule) {
                std::cout << "t[" << i << "]: " << temperature << "\n";
            }
            m.setBeta(beta);
            chains.push_back(m);
        }
    }

}

namespace cvt{ // convergence tools
        /* here I will put the implementation of chi2 test for 2 binned datasets as described in
        @book{nr1985,
         title={Numerical Recipes: Example Book: Fortran},
         author={Press, William H and Vetterling, William T and Teukolsky, Saul A and Flannery, Brian P},
         year={1985},
         publisher={Cambridge Univ. P.},
         pages={471-472}}
         */
}
#endif //MT1DANISMODELPARAMS_OBJECTS_H
