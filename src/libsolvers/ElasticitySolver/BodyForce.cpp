#include "BodyForce.hpp"

using namespace BodyForce;

BodyForce::GravityForce::GravityForce(const double rho,
                                      const Point<DIM> &g) : Function<DIM>(DIM),
    rho(rho), g(g)
{

}

void BodyForce::GravityForce::vector_value(const Point<DIM> &/*p*/, Vector<double> &values) const
{
    Assert(values.size() == DIM, ExcDimensionMismatch(values.size(), DIM));

    values = 0;
    for (size_t i = 0; i < DIM; ++i)
    {
        values(i) = rho * g(i);
    }
}
void BodyForce::GravityForce::vector_value_list(const std::vector<Point<DIM>> &points,
                                                std::vector<Vector<double>> &value_list) const
{
    const size_t n_points = points.size();

    Assert(value_list.size() == n_points, ExcDimensionMismatch(value_list.size(), n_points));

    for (size_t p = 0; p < n_points; ++p) {
        GravityForce::vector_value(points[p], value_list[p]);
    }
}
