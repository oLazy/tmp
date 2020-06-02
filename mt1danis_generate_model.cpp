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
    double tmpd=m.nodes[m.nodes.size()-1].params[mtobj::paramType::depth].getValue()*1.05;
    double tmpsm=m.nodes[m.nodes.size()-1].params[mtobj::paramType::sigmaMean].getValue();
    std::vector<double> x, y, x1, x2;
    std::vector<double> xb, yb, dxb, dyb;
//    x.push_back(m.nodes[0].params[mtobj::paramType::sigmaMean].getValue());
//    y.push_back(0);
//    x1.push_back(m._al[0]*pow(10,-6));

    for(int i=0; i<m.nodes.size();i++){
        x.push_back(m.nodes[i].params[mtobj::paramType::sigmaMean].getValue());
        x1.push_back(m._al[i]);
        x2.push_back(m._at[i]);
        if(m.nodes[i].params[mtobj::paramType::beta].isActive()) {
            xb.push_back(m.nodes[i].params[mtobj::paramType::beta].getValue());
            yb.push_back(m.nodes[i].params[mtobj::paramType::depth].getValue());
            dxb.push_back(0);
            if(i!=m.nodes.size()-1){
//                dyb.push_back(m.nodes[i+1].params[mtobj::paramType::depth].getValue());
                dyb.push_back(m._h[i]);
            }else{
                dyb.push_back(tmpd-m.nodes[i].params[mtobj::paramType::depth].getValue());
            }
        }
        y.push_back(m.nodes[i].params[mtobj::paramType::depth].getValue());
        x.push_back(m.nodes[i].params[mtobj::paramType::sigmaMean].getValue());
        x1.push_back(m._al[i]);
        x2.push_back(m._at[i]);
//        x3.push_back(m._blt[i]);
        if (i!=m.nodes.size()-1){
            y.push_back(m.nodes[i+1].params[mtobj::paramType::depth].getValue());}
        else{
            y.push_back(tmpd);
        }
//        std::cout << "i= " << i <<"\n";

    }
//    x1.push_back(m._al[m.nodes.size()]);
//    std::cout << x1.size() << " " << y.size() << "\n";
//    for (int i=0; i<x1.size(); i++){
//        std::cout << x1[i] << "\t" << y[i] << "\n";
//    }
//    y.push_back(tmpd);
//    x.push_back(tmpsm);

    Gnuplot gp;
    auto gp_filename = out_file_name.replace(out_file_name.length()-3,3,"pdf");
    gp << "set term pdfcairo size 1.75,2 font \"Times,9\"\n";
    gp << "set tics out\n";
    gp << "set format x \"10^{%T}\"\n";
    gp << "set grid\n";
    gp << "set out \"" << gp_filename << "\"\n";
    gp << "set logscale x\n";
    gp << "set xrange [0.000005:15]\n";
    gp << "set yrange ["<<tmpd<<":0] reverse\n";
    gp << "set xtics nomirror\n";
    gp << "set ytics rotate by 45\n";
    gp << "set x2range [-90:90]\n";
    gp << "set x2tics -90,30\n";
    gp << "set xlabel \"Conductivity (S/m)\"\n";
    gp << "set x2label \"Strike ({/Symbol \\260})\"\n";
    gp << "set ylabel \"Depth (m)\"\n";
    gp << "plot '-' with lines lw 2.5 lc rgb '#F0BC42' title '{/Symbol s}_{hi}',"<<
       "'-' with lines lw 1 lc rgb '#8E1F2F' title '{/Symbol s}_{lo}'," <<
       "'-' with vectors nohead lc rgb '#000000' axis x2y1 title '{/Symbol b}_{strike}'\n";
    gp.send1d(boost::make_tuple(x1,y));
    gp.send1d(boost::make_tuple(x2,y));
    gp.send1d(boost::make_tuple(xb,yb,dxb,dyb));


    return 0;
}