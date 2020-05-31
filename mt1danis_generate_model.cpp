//
// Created by Eric Mandolesi on 30/05/2020.
//
#include "objects.h"
#include "global.h"
#include <sstream>
#include <string>
#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>

int main(int argn, char* argv[]){
    if (argn > 2){
        std::cout << "Error! this program requires a single argument.\n";
        return 1;
    }
    std::string out_file_name{"model.dat"}; // TODO make this custom
    std::ifstream is(argv[1]);
    std::string line;
    int ln = 0;
    mtobj::model m;
    while (std::getline(is,line))
    {
        std::istringstream iss(line);
        double z, h, l, b;
        if((iss >> z >> h)){

            if((iss >> l >> b)){
                //process
                m.nodes.push_back(mtobj::node(z, log10(0.5*(h+l)), log10(l/h), b));
            }else{
                //process
                m.nodes.push_back(mtobj::node(z, log10(h)));}
        }
        else{
            std::cerr << "error in input file, line " << ln << ". Please check.\n";
            return -ln;
        }

    }
    // check node 0 has depth 0
    if(m.nodes[0].params[mtobj::paramType::depth].getValue() == 0){
        std::ofstream os(out_file_name);
        boost::archive::text_oarchive oa(os);
        oa << m ;
    }else{
        return 2;
    }
    m.calc_params();
    // plot section: here the model is plot both in sampling space and in physical space
    double tmpd=m.nodes[m.nodes.size()-1].params[mtobj::paramType::depth].getValue()+1000;
    double tmpsm=m.nodes[m.nodes.size()-1].params[mtobj::paramType::sigmaMean].getValue();
    std::vector<double> x, y, x1, x2, x3;
//    x.push_back(m.nodes[0].params[mtobj::paramType::sigmaMean].getValue());
//    y.push_back(0);
//    x1.push_back(m._al[0]*pow(10,-6));

    for(int i=0; i<m.nodes.size();i++){
        x.push_back(m.nodes[i].params[mtobj::paramType::sigmaMean].getValue());
        x1.push_back(m._al[i]);
        x2.push_back(m._at[i]);
        x3.push_back(m._blt[i]);
        y.push_back(m.nodes[i].params[mtobj::paramType::depth].getValue());
        x.push_back(m.nodes[i].params[mtobj::paramType::sigmaMean].getValue());
        x1.push_back(m._al[i]);
        x2.push_back(m._at[i]);
        x3.push_back(m._blt[i]);
        if (i!=m.nodes.size()-1){
            y.push_back(m.nodes[i+1].params[mtobj::paramType::depth].getValue());}
        else{
            y.push_back(tmpd);
        }
        std::cout << "i= " << i <<"\n";

    }
//    x1.push_back(m._al[m.nodes.size()]);
    std::cout << x1.size() << " " << y.size() << "\n";
    for (int i=0; i<x1.size(); i++){
        std::cout << x1[i] << "\t" << y[i] << "\n";
    }
//    y.push_back(tmpd);
//    x.push_back(tmpsm);

    Gnuplot gp;
    auto gp_filename = out_file_name.replace(out_file_name.length()-3,3,"png");
    gp << "set term pngcairo\n";
    gp << "set tics out\n";
    gp << "set grid\n";
    gp << "set out \"" << gp_filename << "\"\n";
    gp << "set logscale x\n";
    gp << "set xrange [0.00001:10]\nset yrange ["<<tmpd<<":0] reverse\n";
    gp << "plot '-' with lines lw 2 title '{/Symbol s}_{hi}', '-' with lines title '{/Symbol s}_{lo}'\n";
    gp.send1d(boost::make_tuple(x1,y));
    gp.send1d(boost::make_tuple(x2,y));


    return 0;
}