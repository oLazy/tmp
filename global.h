//
// Created by Eric Mandolesi on 30/05/2020.
//

#ifndef MT1DANISMODELPARAMS_GLOBAL_H
#define MT1DANISMODELPARAMS_GLOBAL_H
int seed{23};
boost::random::mt19937 gen{static_cast<std::uint32_t>(seed)};
boost::random::normal_distribution<double> rn;

#endif //MT1DANISMODELPARAMS_GLOBAL_H
