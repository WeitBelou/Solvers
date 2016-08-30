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
        this->boundary_functions_map.insert(std::make_pair(item.first, DirichletBoundary(item.second,
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
                                                 temp);
    }

    return temp;
}

//
DirichletBoundary::DirichletBoundary(const std::vector<std::string> &function,
                                     const double timestep)
    :
    Function<DIM>(DIM),
    function(function),
    present_timestep(timestep)
{

}

DirichletBoundary::DirichletBoundary(const DirichletBoundary &other)
    :
    Function<DIM>(DIM),
    function(other.function),
    present_timestep(other.present_timestep)
{

}

DirichletBoundary::~DirichletBoundary()
{

}

void DirichletBoundary::vector_value(const Point<DIM> &p, Vector<double> &values) const
{
    FunctionParser<DIM> function_parser(DIM, get_time());
    const std::string vars("x, y, z, t");
    const std::map<std::string, double> constants = {std::make_pair("dt", present_timestep)};
    function_parser.initialize(vars, function, constants, true);

    function_parser.vector_value(p, values);
}
