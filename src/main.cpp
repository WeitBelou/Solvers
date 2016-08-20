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
    outputDirOption->value_name("output")->default_value("./");

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
    }
    catch (std::exception & exc) {
        std::cerr << "Error during parsing command line: " << std::endl
                  << exc.what() << std::endl;
        return 1;
    }

    Launcher launcher;
    try {
        launcher.set_output_dir(outputDir);
        launcher.set_task(taskFile);
    }
    catch (std::exception & exc) {
        std::cerr << "Error during initializing launcher: " << std::endl
                  << exc.what() << std::endl;
        return 1;
    }

    try {
        launcher.run();
    }
    catch (std::exception & exc) {
        std::cerr << "Runtime error: " << std::endl
                  << exc.what() << std::endl;
    }

    return 0;
}