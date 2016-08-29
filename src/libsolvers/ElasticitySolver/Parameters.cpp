#include "Parameters.hpp"

using namespace ElasticityEquation;


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

void BoundaryConditions::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Boundary conditions");
    {
        prm.declare_entry("Dirichlet boundary conditions",
                          "", Patterns::Map(Patterns::Integer(),
                                            Patterns::Anything()));
    }
    prm.leave_subsection();
}

void BoundaryConditions::parse_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Boundary conditions");
    {
        const std::string function_map = prm.get("Dirichlet boundary conditions");

        const std::vector<std::string> function_id_vector = Utilities::split_string_list(function_map, ',');

        for (auto pair: function_id_vector)
        {
            const std::vector<std::string> function_pair = Utilities::split_string_list(pair, ':');

            const types::boundary_id index = static_cast<types::boundary_id>(Utilities::string_to_int(function_pair[0]));

            const std::vector<std::string> function_components = Utilities::split_string_list(function_pair[1], ';');

            boundary_functions.insert(std::make_pair(index, function_components));
        }
    }
    prm.leave_subsection();
}

Parameters::Parameters(const std::string &input_file)
{
    ParameterHandler prm;
    declare_parameters(prm);
    prm.read_input(input_file);
    parse_parameters(prm);
}

void Parameters::declare_parameters(ParameterHandler &prm)
{
    FiniteElementSystem::declare_parameters(prm);
    Geometry::declare_parameters(prm);
    Material::declare_parameters(prm);
    Time::declare_parameters(prm);
    BoundaryConditions::declare_parameters(prm);
}

void Parameters::parse_parameters(ParameterHandler &prm)
{
    FiniteElementSystem::parse_parameters(prm);
    Geometry::parse_parameters(prm);
    Material::parse_parameters(prm);
    Time::parse_parameters(prm);
    BoundaryConditions::parse_parameters(prm);
}
