//
// Created by ivan on 22.08.16.
//
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "BoundaryConditions.hpp"

using namespace ElasticityEquation;

FunctionTimeBoundaryConditions::FunctionTimeBoundaryConditions(const std::map<types::boundary_id,
                                                                              std::vector<std::string>> &boundary_values_map,
                                                               const double timestep)
    : Subscriptor()
{
    for (auto item : boundary_values_map)
    {
        this->boundary_functions_map.insert(std::make_pair(item.first, BaseBoundary(item.second,
                                                                                    ComponentMask(),
                                                                                    timestep)));
    }
}

void FunctionTimeBoundaryConditions::reinit(double present_time)
{
    for (auto &item : this->boundary_functions_map)
    {
        item.second.set_time(present_time);
    }
}

void FunctionTimeBoundaryConditions::update(double present_timestep)
{
    for (auto &item : this->boundary_functions_map)
    {
        item.second.advance_time(present_timestep);
    }
}

std::map<types::global_dof_index, double>
FunctionTimeBoundaryConditions::interpolate(const DoFHandler<DIM> &dof_handler)
{
    std::map<types::global_dof_index, double> temp;

    for (auto item : boundary_functions_map)
    {
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 item.first,
                                                 item.second,
                                                 temp,
                                                 item.second.get_mask());
    }

    return temp;
}

//
BaseBoundary::BaseBoundary(const std::vector<std::string> &function,
                           const ComponentMask &mask,
                           const double timestep)
    :
    Function<DIM>(DIM),
    mask(mask),
    function(function),
    present_timestep(timestep)
{

}

BaseBoundary::BaseBoundary(const BaseBoundary &other)
    :
    Function<DIM>(DIM),
    mask(other.mask),
    function(other.function),
    present_timestep(other.present_timestep)
{

}

ComponentMask BaseBoundary::get_mask()
{
    return mask;
}

BaseBoundary::~BaseBoundary()
{

}

void BaseBoundary::vector_value(const Point<DIM> &p, Vector<double> &values) const
{
    FunctionParser<DIM> function_parser(DIM, get_time());
    const std::string vars("x, y, z, t");
    const std::map<std::string, double> constants = {std::make_pair("dt", present_timestep)};
    function_parser.initialize(vars, function, constants, true);

    function_parser.vector_value(p, values);
}