//
// Created by ivan on 22.08.16.
//

#ifndef SOLVERS_BOUNDARYCONDITIONS_HPP
#define SOLVERS_BOUNDARYCONDITIONS_HPP

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>

#include "global.hpp"

namespace ElasticityEquation
{
//begin namespace ElasticityEquation
using namespace dealii;

class DirichletBoundary;

class FunctionTimeBoundaryConditions: public Subscriptor
{
public:
    FunctionTimeBoundaryConditions(const std::map<types::boundary_id,
                                                      std::vector<std::string>> &boundary_functions_map,
                                       const double timestep);

    void reinit(double present_time);
    void update(double present_timestep);

    std::map<types::global_dof_index, double> interpolate(const DoFHandler<DIM> &dof_handler);

private:
    std::map<types::boundary_id, DirichletBoundary> boundary_functions_map;
};

class DirichletBoundary: public Function<DIM>
{
public:
    DirichletBoundary(const std::vector<std::string> &function,
                 const ComponentMask &mask = ComponentMask(),
                 const double timestep = 1.0);
    DirichletBoundary(const DirichletBoundary &other);

    ComponentMask get_mask();
    virtual void vector_value(const Point<DIM> &p, Vector<double> &values) const override;

    virtual ~DirichletBoundary();

private:
    const ComponentMask mask;
    const std::vector<std::string> function;

    double present_timestep;
};

//end namespace ElasticityEquation
}


#endif //SOLVERS_BOUNDARYCONDITIONS_HPP
