#include <iostream>
#include <boost/program_options.hpp>
#include "src/launcher/Launcher.h"

namespace po = boost::program_options;

int main(int argc, char ** argv)
{
    po::options_description desc ("Options");

    std::string taskFile;
    std::string outputDir;

    desc.add_options ()
        ("help, h", "Show help")
        ("task", po::value<std::string>(&taskFile), "File with task")
        ("output", po::value<std::string>(&outputDir), "Output directory")
        ;

    po::variables_map vm;
    try {
        po::parsed_options parsed = po::command_line_parser(argc, argv).
                                        options(desc).allow_unregistered().run();
        po::store (parsed, vm);
        po::notify (vm);

        std::cout << taskFile << std::endl;
        std::cout << outputDir << std::endl;
    }
    catch (std::exception & exception) {
        std::cerr << "Error during parsing command line: "
                  << exception.what() << std::endl;
        return 1;
    }

    return 0;
}