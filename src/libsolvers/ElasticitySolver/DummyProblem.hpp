//
// Created by ivan on 20.08.16.
//

#ifndef SOLVERS_DUMMYPROBLEM_HPP
#define SOLVERS_DUMMYPROBLEM_HPP

#include "ElasticitySolver.hpp"
#include "Parameters.hpp"
#include "global.hpp"

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
namespace PipeTask {

using namespace dealii;

void run_pipe_task(const Parameters::All &par);

void write_pipe_grid(const std::string &file_name);
void read_triangulation(Triangulation<DIM> &triangulation,
                        std::string file_name);

}

#endif //SOLVERS_DUMMYPROBLEM_HPP
