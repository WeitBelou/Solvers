//
// Created by ivan on 20.08.16.
//

#include "DummyProblem.hpp"

using namespace DummyProblem;

void ::DummyProblem::run_pipe_task()
{
    Triangulation<dim> triangulation;
    const double inner_radius = 0.8,
        outer_radius = 1;
    GridGenerator::cylinder_shell(triangulation,
                                  3, inner_radius, outer_radius);
    for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
    {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {
            if (cell->face(f)->at_boundary())
            {
                const Point<dim> face_center = cell->face(f)->center();
                if (face_center[2] == 0)
                {
                    cell->face(f)->set_boundary_id(0);
                }
                else if (face_center[2] == 3)
                {
                    cell->face(f)->set_boundary_id(1);
                }
                else if (std::sqrt(face_center[0] * face_center[0] +
                                   face_center[1] * face_center[1])
                         <
                         (inner_radius + outer_radius) / 2)
                {
                    cell->face(f)->set_boundary_id(2);
                }
                else
                {
                    cell->face(f)->set_boundary_id(3);
                }
            }
        }
    }
    static const CylindricalManifold<dim> cylindrical_manifold(2);
    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold(0, cylindrical_manifold);

    FESystem<dim> fe(FE_Q<dim>(1), dim);
    QGauss<dim> quadrature(2);
    ElasticitySolver::IncrementalBoundaryValues incremental_boundary_values;
    BodyForce body_force;

    ElasticitySolver::TopLevel top_level(triangulation, fe, quadrature, body_force, incremental_boundary_values);
    top_level.run();
}

BodyForce::BodyForce() : Function<dim>(dim)
{

}

void BodyForce::vector_value(const Point<dim> &/*p*/, Vector<double> &values) const
{
    Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));

    const double g = 9.81;
    const double rho = 7700;

    values = 0;
    values(dim - 1) = -rho * g;
}
void BodyForce::vector_value_list(const std::vector<Point<dim>> &points, std::vector<Vector<double>> &value_list) const
{
    const size_t n_points = points.size();

    Assert(value_list.size() == n_points, ExcDimensionMismatch(value_list.size(), n_points));

    for (size_t p = 0; p < n_points; ++p)
    {
        BodyForce::vector_value(points[p], value_list[p]);
    }
}
