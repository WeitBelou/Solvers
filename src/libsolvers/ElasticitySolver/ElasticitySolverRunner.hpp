//
// Created by ivan on 20.08.16.
//

#ifndef SOLVERS_ELASTICITYSOLVERRUNNER_HPP
#define SOLVERS_ELASTICITYSOLVERRUNNER_HPP

#include "ElasticitySolver.hpp"
#include "src/launcher/util/Parameters.hpp"
#include "src/libsolvers/global.hpp"

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

namespace PipeTask {

void run_pipe_task(const Parameters::ElasticitySolverParameters &par);

void write_pipe_grid(const std::string &file_name);
void read_triangulation(dealii::Triangulation<DIM> &triangulation,
                        std::string file_name);

}

#endif //SOLVERS_ELASTICITYSOLVERRUNNER_HPP
