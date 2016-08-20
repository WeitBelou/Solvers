//
// Created by ivan on 08.08.16.
//
#include <iostream>

#include "Launcher.h"
#include "src/libsolvers/ElasticitySolver/DummyProblem.hpp"

Launcher::Launcher ()
{

}

void Launcher::set_task (const std::string &taskName)
{
    this->taskName = taskName;
    boost::system::error_code err;
    if (!bfs::exists (this->taskName, err))
    {
        throw bfs::filesystem_error (taskName, err);
    }
}

void Launcher::set_output_dir (const std::string &outputDirName)
{
    outputDir = outputDirName;

    if (!bfs::exists (outputDir))
    {
        std::cout << "Output direction does not exist. Trying to create it."
                  << std::endl;
        bfs::create_directories (outputDir);
    }
}

void Launcher::run ()
{
    std::cout << "Running task: " << taskName.filename ()
              << " from: " << bfs::canonical (taskName).remove_filename ()
              << std::endl;

    DummyProblem::run_pipe_task();

    std::cout << "Output will be writen in: " << bfs::canonical (outputDir)
              << std::endl;

    std::cout << "Task completed." << std::endl;
}

