//
// Created by ivan on 20.08.16.
//

#ifndef SOLVERS_DUMMYPROBLEM_HPP
#define SOLVERS_DUMMYPROBLEM_HPP

#include "ElasticitySolver.hpp"

namespace DummyProblem
{
const int dim = 3;

using namespace dealii;

void run_pipe_task();

void set_triangulation(Triangulation<dim> &triangulation);

class BodyForce: public Function<dim>
{
public:
    BodyForce();
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &value_list) const override;

};

}

#endif //SOLVERS_DUMMYPROBLEM_HPP
