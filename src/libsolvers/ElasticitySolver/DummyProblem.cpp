//
// Created by ivan on 20.08.16.
//

#include "DummyProblem.hpp"
#include "BodyForce.hpp"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/fe/fe_q.h>

using namespace PipeTask;
namespace es = ElasticityEquation;

void ::PipeTask::run_pipe_task(const ElasticityEquation::Parameters &par)
{
    Triangulation<DIM> triangulation;
    read_triangulation(triangulation, par.path_to_grid);

    FESystem<DIM> fe(FE_Q<DIM>(par.polynomial_degree), DIM);
    QGauss<DIM> quadrature(par.quadrature_degree);

    es::GravityForce body_force;

    FEValuesExtractors::Scalar z_component(DIM - 1);
    ComponentMask z_mask = fe.component_mask(z_component);
    es::FunctionTimeBoundaryConditions boundary_conditions(par.boundary_functions, par.timestep);

    es::ElasticitySolver top_level(triangulation, fe, quadrature,
                                           body_force, boundary_conditions);
    top_level.run(par.timestep, par.end_time);
}

void PipeTask::write_pipe_grid(const std::string &file_name)
{
    Triangulation<DIM> triangulation;

    const double inner_radius = 0.8,
        outer_radius = 1;
    GridGenerator::cylinder_shell(triangulation,
                                  3, inner_radius, outer_radius);
    for (typename Triangulation<DIM>::active_cell_iterator
             cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
    {
        for (unsigned int f = 0; f < GeometryInfo<DIM>::faces_per_cell; ++f)
        {
            if (cell->face(f)->at_boundary())
            {
                const Point<DIM> face_center = cell->face(f)->center();
                if (face_center[2] == 0)
                {
                    cell->face(f)->set_boundary_id(0);
                }
                else if (face_center[2] == 3)
                {
                    cell->face(f)->set_boundary_id(1);
                }
                else if (sqrt(face_center[0] * face_center[0] +
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
    static const CylindricalManifold<DIM> cylindrical_manifold(2);
    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold(0, cylindrical_manifold);
    triangulation.refine_global(1);

    GridOut grid_out;

    std::ofstream out(file_name);

    GridOutFlags::Msh msh_flags(true, true);
    grid_out.set_flags(msh_flags);

    grid_out.write_msh(triangulation, out);
}

void PipeTask::read_triangulation(Triangulation<DIM> &triangulation,
                                  std::string file_name)
{
    GridIn<DIM> grid_in;
    grid_in.attach_triangulation(triangulation);

    std::ifstream in(file_name);

    grid_in.read_msh(in);
}
