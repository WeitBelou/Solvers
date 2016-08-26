//
// Created by ivan on 20.08.16.
//

#include "DummyProblem.hpp"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_q.h>

using namespace DummyProblem;
namespace bc = BoundaryConditions;

void ::DummyProblem::run_pipe_task()
{
    Triangulation<DIM> triangulation;
    set_triangulation(triangulation);

    FESystem<DIM> fe(FE_Q<DIM>(1), DIM);
    QGauss<DIM> quadrature(2);

    BodyForce body_force;

    FEValuesExtractors::Scalar z_component(DIM - 1);
    ComponentMask z_mask = fe.component_mask(z_component);
    bc::IncrementalBoundaryValues inc_bv({0, 0, -0.1}, z_mask);
    bc::ZeroFunctionBoundaryValues zero_bv;
    bc::FunctionTimeBoundaryConditions boundary_conditions( {
        std::make_pair(0, &zero_bv),
        std::make_pair(1, &inc_bv)
    });

    ElasticitySolver::TopLevel top_level(triangulation, fe, quadrature,
                                         body_force, boundary_conditions);
    top_level.run();
}

void DummyProblem::set_triangulation(Triangulation<DIM> &triangulation)
{
    const double inner_radius = 0.8,
                 outer_radius = 1;
    GridGenerator::cylinder_shell(triangulation,
                                  3, inner_radius, outer_radius);
    for (typename Triangulation<DIM>::active_cell_iterator
            cell = triangulation.begin_active();
            cell != triangulation.end(); ++cell) {
        for (unsigned int f = 0; f < GeometryInfo<DIM>::faces_per_cell; ++f) {
            if (cell->face(f)->at_boundary()) {
                const Point<DIM> face_center = cell->face(f)->center();
                if (face_center[2] == 0) {
                    cell->face(f)->set_boundary_id(0);
                } else if (face_center[2] == 3) {
                    cell->face(f)->set_boundary_id(1);
                } else if (sqrt(face_center[0] * face_center[0] +
                                face_center[1] * face_center[1])
                           <
                           (inner_radius + outer_radius) / 2) {
                    cell->face(f)->set_boundary_id(2);
                } else {
                    cell->face(f)->set_boundary_id(3);
                }
            }
        }
    }
    static const CylindricalManifold<DIM> cylindrical_manifold(2);
    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold(0, cylindrical_manifold);
    triangulation.refine_global(1);
}

BodyForce::BodyForce() : Function<DIM>(DIM)
{

}

void BodyForce::vector_value(const Point<DIM> &/*p*/, Vector<double> &values) const
{
    Assert(values.size() == DIM, ExcDimensionMismatch(values.size(), DIM));

    const double g = 9.81;
    const double rho = 7700;

    values = 0;
    values(DIM - 1) = -rho * g;
}
void BodyForce::vector_value_list(const std::vector<Point<DIM>> &points,
                                  std::vector<Vector<double>> &value_list) const
{
    const size_t n_points = points.size();

    Assert(value_list.size() == n_points, ExcDimensionMismatch(value_list.size(), n_points));

    for (size_t p = 0; p < n_points; ++p) {
        BodyForce::vector_value(points[p], value_list[p]);
    }
}
