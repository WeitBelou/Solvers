#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <deal.II/base/parameter_handler.h>

#include "BoundaryConditions.hpp"

namespace ElasticityEquation
{
using namespace dealii;

class FiniteElementSystem
{
public:
    size_t polynomial_degree;
    size_t quadrature_degree;

    static void declare_parameters(ParameterHandler &prm);
    void parse_parameters(ParameterHandler &prm);
};

class Geometry
{
public:
    std::string path_to_grid;

    static void declare_parameters(ParameterHandler &prm);
    void parse_parameters(ParameterHandler &prm);
};

struct BoundaryConditions
{
    std::map<types::boundary_id, std::vector<std::string>> boundary_functions;

    static void declare_parameters(ParameterHandler &prm);
    void parse_parameters(ParameterHandler &prm);
};

class Material
{
public:
    double lambda;
    double mu;

    static void declare_parameters(ParameterHandler &prm);
    void parse_parameters(ParameterHandler &prm);
};

class Time
{
public:
    double timestep;
    double end_time;

    static void declare_parameters(ParameterHandler &prm);
    void parse_parameters(ParameterHandler &prm);
};

class Parameters: public FiniteElementSystem,
                  public Geometry,
                  public Material,
                  public Time,
                  public BoundaryConditions
{
public:
    Parameters(const std::string &input_file);

    static void declare_parameters(ParameterHandler &prm);
    void parse_parameters(ParameterHandler &prm);
};

}

#endif // PARAMETERS_HPP
