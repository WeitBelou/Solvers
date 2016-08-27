#include "Parameters.hpp"

using namespace Parameters;


void FiniteElementSystem::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Finite element system");
    {
        prm.declare_entry("Polynomial degree", "1",
                          Patterns::Integer(0),
                          "Displacement system polynomial order");
        prm.declare_entry("Quadrature order", "2",
                          Patterns::Integer(0),
                          "Gauss quadrature order");
    }
    prm.leave_subsection();
}

void FiniteElementSystem::parse_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Finite element system");
    {
        polynomial_degree = prm.get_integer("Polynomial degree");
        quadrature_degree = prm.get_integer("Quadrature order");
    }
    prm.leave_subsection();
}

void Geometry::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Geometry");
    {
        prm.declare_entry("Path to grid", "",
                          Patterns::FileName(),
                          "Path to file that contains grid");
    }
    prm.leave_subsection();
}

void Geometry::parse_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Geometry");
    {
        path_to_grid = prm.get("Path to grid");
    }
    prm.leave_subsection();
}

void Material::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Material");
    {
        prm.declare_entry("Lambda", "9.695e10",
                          Patterns::Double(0),
                          "Lambda Lame paramter");
        prm.declare_entry("Mu", "7.617e10",
                          Patterns::Double(0),
                          "Mu Lame paramter");
    }
    prm.leave_subsection();
}

void Material::parse_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Material");
    {
        lambda = prm.get_double("Lambda");
        mu = prm.get_double("Mu");
    }
    prm.leave_subsection();
}

void Time::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Time");
    {
        prm.declare_entry("Timestep", "1.0",
                          Patterns::Double(0.01),
                          "Timestep");
        prm.declare_entry("End time", "10.0",
                          Patterns::Double(1.0, 15.0),
                          "End time");
    }
    prm.leave_subsection();
}

void Time::parse_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Time");
    {
        timestep = prm.get_double("Timestep");
        end_time = prm.get_double("End time");
    }
    prm.leave_subsection();
}

All::All(const std::string &input_file)
{
    ParameterHandler prm;
    declare_parameters(prm);
    prm.read_input(input_file);
    parse_parameters(prm);
}

void All::declare_parameters(ParameterHandler &prm)
{
    FiniteElementSystem::declare_parameters(prm);
    Geometry::declare_parameters(prm);
    Material::declare_parameters(prm);
    Time::declare_parameters(prm);
}

void All::parse_parameters(ParameterHandler &prm)
{
    FiniteElementSystem::parse_parameters(prm);
    Geometry::parse_parameters(prm);
    Material::parse_parameters(prm);
    Time::parse_parameters(prm);
}
