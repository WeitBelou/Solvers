//
// Created by ivan on 22.08.16.
//

#ifndef SOLVERS_BOUNDARYCONDITIONS_HPP
#define SOLVERS_BOUNDARYCONDITIONS_HPP

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>

#include "src/libsolvers/global.hpp"

namespace BoundaryConditions
{
using namespace dealii;

class DirichletBoundary;
class BoundaryMaskGroup;

class FunctionBoundaryConditions: public Subscriptor
{
public:
    FunctionBoundaryConditions(const std::map<types::boundary_id, std::string> &boundary_functions_map,
                                   const std::map<types::boundary_id,
                                                  std::string> &bondary_functions_mask,
                                   const BoundaryMaskGroup &mask_group,
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
    DirichletBoundary(const std::string &function,
                      const ComponentMask &mask,
                      const double timestep);
    DirichletBoundary(const DirichletBoundary &other);

    virtual void vector_value(const Point<DIM> &p, Vector<double> &values) const override;
    ComponentMask get_mask() const;

    virtual ~DirichletBoundary();

private:
    const std::string function;
    const ComponentMask mask;
    double present_timestep;
};

class BoundaryMaskGroup
{
public:
    BoundaryMaskGroup(const ComponentMask &x_mask,
                      const ComponentMask &y_mask,
                      const ComponentMask &z_mask);

    ComponentMask get_mask_from_string(const std::string &mask) const;
private:
    const ComponentMask x_mask;
    const ComponentMask y_mask;
    const ComponentMask z_mask;
};

}


#endif //SOLVERS_BOUNDARYCONDITIONS_HPP
