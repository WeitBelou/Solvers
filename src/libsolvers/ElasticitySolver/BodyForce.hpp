#ifndef BODYFORCE_HPP
#define BODYFORCE_HPP

#include "global.hpp"

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

namespace BodyForce {
using namespace dealii;

class GravityForce: public Function<DIM>
{
public:
    GravityForce(const double rho = 7700,
                 const Point<DIM> & g = Point<DIM>(0.0, 0.0, -9.81));
    virtual void
    vector_value(const Point<DIM> &p, Vector<double> &values) const override;

    virtual void
    vector_value_list(const std::vector<Point<DIM>> &points,
                      std::vector<Vector<double>> &value_list) const override;
private:
    const double rho;
    const Point<DIM> g;
};

}

#endif // BODYFORCE_HPP
