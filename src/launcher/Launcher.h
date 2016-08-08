//
// Created by ivan on 08.08.16.
//

#ifndef ELASTICITY_LAUNCHER_H
#define ELASTICITY_LAUNCHER_H

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

class Launcher
{
public:
    Launcher ();

    void set_task (const std::string &taskName);
    void set_output_dir (const std::string & outputDirName);

    void run ();
private:
    bfs::path outputDir;
    bfs::path taskName;
};


#endif //ELASTICITY_LAUNCHER_H
