//
// Created by ivan on 22.08.16.
//

#ifndef SOLVERS_BOUNDARYCONDITIONS_HPP
#define SOLVERS_BOUNDARYCONDITIONS_HPP

#include <deal.II/base/function_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/dofs/function_map.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/component_mask.h>

#include "global.hpp"

namespace BoundaryConditions
{
//begin namespace BoundaryConditions
using namespace dealii;

enum BoundaryType
{
    Dirichlet,
    Neumann
};

class BaseBoundary;
class DirichletBoundary;
class NeumannBoundary;

class FunctionTimeBoundaryConditions : public Subscriptor
{
public:
    FunctionTimeBoundaryConditions(const std::map<types::boundary_id,
                                                  BaseBoundary *> &boundary_functions_map);

    void reinit (double present_time);
    void update (double present_timestep);

    std::map<types::global_dof_index, double> interpolate(const DoFHandler<DIM> & dof_handler);

    BaseBoundary *function(types::boundary_id id);
private:
    std::map<types::boundary_id, BaseBoundary *> boundary_functions_map;
};

class BaseBoundary : public Function<DIM>
{
public:
    BaseBoundary(const ComponentMask & mask);

    ComponentMask get_mask();

    virtual ~BaseBoundary();

    virtual BoundaryType type () const = 0;

private:
    const ComponentMask mask;
};

class DirichletBoundary : public BaseBoundary
{
public:
    DirichletBoundary (const ComponentMask & mask = ComponentMask());
    virtual ~DirichletBoundary();

    virtual BoundaryType type () const override;
};

class NeumannBoundary : public BaseBoundary
{
public:
    NeumannBoundary (const ComponentMask & mask = ComponentMask());
    virtual ~NeumannBoundary();

    virtual BoundaryType type () const override;
};

class IncrementalBoundaryValues: public DirichletBoundary
{
public:
    IncrementalBoundaryValues(const double velocity, const ComponentMask & mask = ComponentMask());

    virtual void advance_time(double present_timestep) override;

    virtual void
    vector_value(const Point<DIM> &p,
                 Vector<double> &values) const override;
    virtual void
    vector_value_list(const std::vector<Point<DIM>> &points,
                      std::vector<Vector<double>> &value_list) const override;
private:
    const double velocity;
    double present_timestep;
};

class ZeroFunctionBoundaryValues : public DirichletBoundary
{
public:
    ZeroFunctionBoundaryValues (const ComponentMask & mask = ComponentMask());

    virtual void
    vector_value(const Point<DIM> &,
                 Vector<double> &values) const override;
};

//end namespace BoundaryConditions
}


#endif //SOLVERS_BOUNDARYCONDITIONS_HPP
