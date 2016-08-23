//
// Created by ivan on 20.08.16.
//

#ifndef SOLVERS_DUMMYPROBLEM_HPP
#define SOLVERS_DUMMYPROBLEM_HPP

#include "ElasticitySolver.hpp"

#include "global.hpp"
namespace DummyProblem
{

using namespace dealii;

void run_pipe_task();

void set_triangulation(Triangulation<DIM> &triangulation);

class BodyForce: public Function<DIM>
{
public:
    BodyForce();
    virtual void
    vector_value(const Point<DIM> &p, Vector<double> &values) const override;

    virtual void
    vector_value_list(const std::vector<Point<DIM>> &points,
                      std::vector<Vector<double>> &value_list) const override;

};

}

#endif //SOLVERS_DUMMYPROBLEM_HPP
