#include <iostream>
#include <boost/program_options.hpp>
#include "src/launcher/Launcher.h"

namespace po = boost::program_options;

int main(int argc, char ** argv)
{
    po::options_description desc ("Options");

    std::string taskFile;
    std::string outputDir;

    auto taskFileOption = po::value<std::string>(&taskFile);
    taskFileOption->value_name("task")->required();

    auto outputDirOption = po::value<std::string>(&outputDir);
    outputDirOption->value_name("output")->default_value(".");

    desc.add_options ()
        ("help,H"   ,                  "Show help")
        ("task,T"   , taskFileOption , "File with task")
        ("output,O" , outputDirOption, "Output directory")
        ;

    try {
        po::variables_map vm;
        po::store (po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc << std::endl;
            return 0;
        }

        po::notify (vm);

        std::cout << taskFile << std::endl;
        std::cout << outputDir << std::endl;
    }
    catch (std::exception & exception) {
        std::cerr << "Error during parsing command line: " << std::endl
                  << exception.what() << std::endl;
        return 1;
    }

    return 0;
}