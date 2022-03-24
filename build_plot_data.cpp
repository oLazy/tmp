//
// Created by eric on 28/02/2022.
//
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/histogram.hpp>
#include "objects.h"
#include "global.h"
#include <boost/archive/binary_iarchive.hpp>
int main(int argc, char* argv[]){
    using namespace boost::filesystem;
    using namespace boost::histogram;
    auto pwd = current_path();
    auto fileName = std::string(argv[1]);
    path out (current_path().append("plot_data"));
    if(!is_regular_file(argv[1]))std::cout << argv[1] << "file not found.\n";
    if (!is_directory(out)) {
        std::cout << "making plot_data folder...\n";
        create_directories(out);
    };
    // create output folder
    std::cout << out << "\n";
    if (is_directory(out))std::cout << out << " exists!\n";
    const unsigned n_z_bins = 4096;
    const unsigned n_sm_bins = 512;
    auto z = axis::regular<>(n_z_bins,0,300000,"depth");
    auto sm = axis::regular<>(n_sm_bins, -8, 2, "sigma_mean");
    auto sr = axis::regular<>(n_sm_bins, -6, 0, "sigma_ratio");
    auto bs = axis::regular<>(180, -90, 90, "angle");
    auto logl = axis::regular<>(21,920.297-(2*920.297),920.297+(2*920.297),"logL");

    std::string extension = boost::filesystem::extension(fileName);
    if (extension != ".bin") throw std::runtime_error("only binary files accepted.\n");
    std::ifstream is(fileName);
    boost::archive::binary_iarchive ia(is);
    mtobj::model model;
    bool done{false};
    unsigned m_read{0};
    auto anis_histogram = make_histogram(z);
    auto interf_histogram = make_histogram(z);
    auto sigma_mean_histogram = make_histogram(sm, z);
    auto sigma_ratio_histogram = make_histogram(sr, z);
    auto beta_strike_histogram = make_histogram(bs, z);
    auto logl_histogram = make_histogram(logl);
    while (!done) {
        try {
            ia >> model;
            m_read++;
            // fill one-d histogram
            for (auto i = 0; i<n_z_bins; i++){
                auto depth = z.bin(i).center();
                if(model.getNode(depth).params[mtobj::paramType::beta].isActive()){
                    anis_histogram(depth);
                    sigma_ratio_histogram(model.getNode(depth).params[mtobj::paramType::sigmaRatio].getValue(), depth);
                    beta_strike_histogram(model.getNode(depth).params[mtobj::paramType::beta].getValue(), depth);
                }
                auto this_node_z = model.getNode(depth).params[mtobj::paramType::depth].getValue();
                if(z.bin(i).lower() <= this_node_z && this_node_z < z.bin(i).upper()){
                    interf_histogram(depth);
                }
                sigma_mean_histogram(model.getNode(depth).params[mtobj::paramType::sigmaMean].getValue(), depth);

            }

            logl_histogram(model.logL);
        }
        catch (boost::archive::archive_exception &e){
            std::cout << "catch exception " << e.what() <<"\n";
            done=true;
            std::cout << m_read << " models read from file " << fileName << "\n";
            break;
        }
        catch (...) {
            std::cout << "Unrecognised problem.\n";
            throw std::bad_exception();
        }
    }

    gp_utils::d1hist2disk(anis_histogram, out.append("anis_prob.dat").string(), m_read);
    gp_utils::d1hist2disk(logl_histogram, out.parent_path().append("logl_hist.dat").string(), false);
    gp_utils::d1hist2disk(interf_histogram, out.parent_path().append("interf_prob.dat").string(), m_read);
    gp_utils::d2hist2disk(sigma_mean_histogram, out.parent_path().append("sigma_mean.dat").string(), n_sm_bins, true);
    gp_utils::d2hist2disk(sigma_ratio_histogram, out.parent_path().append("sigma_ratio.dat").string(), n_sm_bins, true);
    gp_utils::d2hist2disk(beta_strike_histogram, out.parent_path().append("beta_strike.dat").string(), n_sm_bins, true);
    return 0;
}