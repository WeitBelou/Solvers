//
// Created by ivan on 10.08.16.
//

#ifndef SOLVERS_ELASTICITYSOLVER_HPP
#define SOLVERS_ELASTICITYSOLVER_HPP

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/base/multithread_info.h>

#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace ElasticitySolver
{
//begin namespace ElasticitySolver
using namespace dealii;

static const int dim = 3;

struct PointHistory
{
    SymmetricTensor<2, dim> old_stress;
};

class IncrementalBoundaryValues: public Function<dim>
{
public:
    IncrementalBoundaryValues();

    void reinit(double present_time, double present_timestep);

    virtual void
    vector_value(const Point<dim> &p,
                 Vector<double> &values) const;
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &value_list) const;
private:
    const double velocity;
    double present_time;
    double present_timestep;
};


inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values, const size_t shape_func, const size_t q_point)
{
    SymmetricTensor<2, dim> strain;

    for (size_t i = 0; i < dim; ++i)
    {
        strain[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];
    }

    for (size_t i = 0; i < dim; ++i)
    {
        for (size_t j = i + 1; j < dim; ++j)
        {
            strain[i][j] = (fe_values.shape_grad_component(shape_func, q_point, i)[j]
                            + fe_values.shape_grad_component(shape_func, q_point, j)[i]) / 2;
        }
    }

    return strain;
}

inline SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda, const double mu)
{
    auto d = [](size_t i, size_t j)
    {
        return (i == j) ? (1.0) : (0.0);
    };

    SymmetricTensor<4, dim> stress_strain;
    for (size_t i = 0; i < dim; ++i)
    {
        for (size_t j = 0; j < dim; ++j)
        {
            for (size_t k = 0; k < dim; ++k)
            {
                for (size_t l = 0; l < dim; ++l)
                {
                    stress_strain[i][j][k][l] = mu * (d(i, k) * d(j, l) + d(i, l) * d(j, k))
                                                + lambda * d(i, j) * d(k, l);
                }
            }
        }
    }
    return stress_strain;
}

inline SymmetricTensor<2, dim> get_strain(const std::vector<Tensor<1, dim>> &grad)
{
    SymmetricTensor<2, dim> strain;

    Assert(grad.size() == dim, ExcDimensionMismatch(grad.size(), dim));

    for (size_t i = 0; i < dim; ++i)
    {
        strain[i][i] = grad[i][i];
    }

    for (size_t i = 0; i < dim; ++i)
    {
        for (size_t j = i + 1; j < dim; ++j)
        {
            strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
        }
    }

    return strain;
}

inline Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>> &grad_u)
{
    const Point<3> curl(grad_u[2][1] - grad_u[1][2],
                        grad_u[0][2] - grad_u[2][0],
                        grad_u[1][0] - grad_u[0][1]);

    const double angle = std::atan(curl.norm());

    if (angle < 1e-9)
    {
        static const double rotation[3][3] = {{1, 0, 0},
                                              {0, 1, 0},
                                              {0, 0, 1}};
        static const Tensor<2, 3> rot(rotation);
        return rot;
    }

    const double s = std::sin(angle);
    const double c = std::cos(angle);
    const double t = 1 - c;

    const Point<3> axis = curl / curl.norm();

    const double rotation[3][3] =
        {{
             t * axis[0] * axis[0] + c,
             t * axis[0] * axis[1] + s * axis[2],
             t * axis[0] * axis[2] - s * axis[1]
         },
         {
             t * axis[0] * axis[1] - s * axis[2],
             t * axis[0] * axis[1] + c,
             t * axis[1] * axis[2] + s * axis[0]
         },
         {
             t * axis[0] * axis[2] + s * axis[1],
             t * axis[1] * axis[1] - s * axis[0],
             t * axis[2] * axis[2] + c
         }
        };

    return Tensor<2, 3, double>(rotation);
}

class TopLevel
{
public:
    TopLevel(Triangulation <dim> &triangulation,
                 const FESystem <dim> &fe,
                 const QGauss <dim> &quadrature,
                 const Function <dim> &body_force,
                 IncrementalBoundaryValues &incremental_boundary_values);
    ~TopLevel();
    void run();

private:
    const SmartPointer<Triangulation<dim>> triangulation;
    const SmartPointer<const FESystem<dim>> fe;
    const SmartPointer<const QGauss<dim>> quadrature;

    DoFHandler<dim> dof_handler;

    const SmartPointer<IncrementalBoundaryValues> incremental_boundary_values;
    const SmartPointer<const Function<dim>> body_force;

    Vector<double> incremental_displacement;

    std::vector<PointHistory> quadrature_point_history;
    static const SymmetricTensor<4, dim> stress_strain_tensor;

    double present_time;
    double present_timestep;
    double end_time;

    size_t timestep_no;
    void refine_grid();
    void setup_quadrature_point_history();
    void update_quadrature_point_history();

    void do_initial_timestep();
    void do_timestep();

    void solve_timestep();
    void move_mesh();
    void output_results() const;

    class LinearSystem
    {
    public:
        LinearSystem(const DoFHandler<dim> &dof_handler);

        void solve(Vector<double> &solution) const;

        ConstraintMatrix hanging_node_constraints;
        SparsityPattern sparsity_pattern;
        SparseMatrix<double> matrix;
        Vector<double> rhs;
    };

    struct AssemblyScratchData
    {
        AssemblyScratchData(const FiniteElement <dim> &fe,
                                    const Quadrature <dim> &quadrature);
        AssemblyScratchData(const AssemblyScratchData &scratch);
        FEValues<dim> fe_values;
    };

    struct AssemblyCopyData
    {
        FullMatrix<double> cell_matrix;
        Vector<double> cell_rhs;

        std::vector<types::global_dof_index> local_dofs_indices;
    };

    void assemble_linear_system(LinearSystem &linear_system);

    void local_assemble_system(const typename DoFHandler<dim>::active_cell_iterator &cell,
                               AssemblyScratchData &scratch_data,
                               AssemblyCopyData &copy_data) const;

    void copy_local_to_global(const AssemblyCopyData &copy_data,
                              LinearSystem &linear_system) const;
};

//end namespace ElasticitySolver
}


#endif //SOLVERS_ELASTICITYSOLVER_HPP
