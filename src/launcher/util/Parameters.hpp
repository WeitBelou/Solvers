#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <deal.II/base/parameter_handler.h>

namespace Parameters
{

class FiniteElementSystem
{
public:
    size_t polynomial_degree;
    size_t quadrature_degree;

    static void declare_parameters(dealii::ParameterHandler &prm);
    void parse_parameters(dealii::ParameterHandler &prm);
};

class Geometry
{
public:
    std::string path_to_grid;

    static void declare_parameters(dealii::ParameterHandler &prm);
    void parse_parameters(dealii::ParameterHandler &prm);
};

struct BoundaryConditions
{
    std::map<dealii::types::boundary_id, std::string> boundary_functions;
    std::map<dealii::types::boundary_id, std::string> boundary_conditions_mask;

    static void declare_parameters(dealii::ParameterHandler &prm);
    void parse_parameters(dealii::ParameterHandler &prm);
};

class Material
{
public:
    double lambda;
    double mu;

    static void declare_parameters(dealii::ParameterHandler &prm);
    void parse_parameters(dealii::ParameterHandler &prm);
};

class Time
{
public:
    double timestep;
    double end_time;

    static void declare_parameters(dealii::ParameterHandler &prm);
    void parse_parameters(dealii::ParameterHandler &prm);
};

class ElasticitySolverParameters: public FiniteElementSystem,
                  public Geometry,
                  public Material,
                  public Time,
                  public BoundaryConditions
{
public:
    ElasticitySolverParameters(const std::string &input_file);

    static void declare_parameters(dealii::ParameterHandler &prm);
    void parse_parameters(dealii::ParameterHandler &prm);
};

}

#endif // PARAMETERS_HPP
