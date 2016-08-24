//
// Created by ivan on 22.08.16.
//
#include <algorithm>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>

#include "BoundaryConditions.hpp"

using namespace BoundaryConditions;

FunctionTimeBoundaryConditions::FunctionTimeBoundaryConditions(const std::map<types::boundary_id,
                                                               BaseBoundary *> &boundary_values_map)
    : Subscriptor()
{
    this->boundary_functions_map = boundary_values_map;

    for (auto item : this->boundary_functions_map)
    {
        item.second->set_time(0.0);
    }
}

void FunctionTimeBoundaryConditions::reinit(double present_time)
{
    for (auto item : this->boundary_functions_map)
    {
        item.second->set_time(present_time);
    }
}

void FunctionTimeBoundaryConditions::update(double present_timestep)
{
    for (auto item : this->boundary_functions_map)
    {
        item.second->advance_time(present_timestep);
    }
}

std::map<types::global_dof_index, double>
FunctionTimeBoundaryConditions::interpolate(const DoFHandler<DIM> &dof_handler)
{
    std::map<types::global_dof_index, double> temp;

    for (auto item : boundary_functions_map)
    {
        VectorTools::interpolate_boundary_values (dof_handler,
                                                  item.first,
                                                  *item.second,
                                                  temp,
                                                  item.second->get_mask());
    }

    return temp;
}

BaseBoundary *FunctionTimeBoundaryConditions::function(types::boundary_id id)
{
    return boundary_functions_map.at(id);
}

//
BaseBoundary::BaseBoundary(const ComponentMask & mask)
    :
    Function<DIM>(DIM),
    mask(mask)
{

}

ComponentMask BaseBoundary::get_mask()
{
    return mask;
}

BaseBoundary::~BaseBoundary()
{

}

//
IncrementalBoundaryValues::IncrementalBoundaryValues(const Point<DIM> &velocity, const ComponentMask &mask)
    :
    BaseBoundary(mask),
    velocity(velocity)
{
    present_timestep = 1.0;
}

void IncrementalBoundaryValues::advance_time(double present_timestep)
{
    this->present_timestep = present_timestep;
    BaseBoundary::advance_time(present_timestep);
}
void IncrementalBoundaryValues::vector_value(const Point<DIM> &/*p*/, Vector<double> &values) const
{
    Assert(values.size() == DIM, ExcDimensionMismatch(values.size(), DIM));
    for (size_t i = 0; i < DIM; ++i)
    {
        values(i) = velocity(i) * present_timestep;
    }
}
void IncrementalBoundaryValues::vector_value_list(const std::vector<Point<DIM>> &points,
                                                  std::vector<Vector<double>> &value_list) const
{
    const unsigned int n_points = points.size();
    Assert (value_list.size() == n_points,
            ExcDimensionMismatch(value_list.size(), n_points));
    for (unsigned int p = 0; p < n_points; ++p)
    {
        IncrementalBoundaryValues::vector_value(points[p],
                                                value_list[p]);
    }
}

//

ZeroFunctionBoundaryValues::ZeroFunctionBoundaryValues(const ComponentMask &mask)
    :
    BaseBoundary(mask)
{

}

void ZeroFunctionBoundaryValues::vector_value(const Point<DIM> &/*p*/, Vector<double> &values) const
{
    Assert (values.size() == DIM, ExcDimensionMismatch(values.size(), DIM));

    values = 0;
}
