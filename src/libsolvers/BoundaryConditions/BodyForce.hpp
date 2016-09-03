#ifndef BODYFORCE_HPP
#define BODYFORCE_HPP

#include "src/libsolvers/global.hpp"

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

namespace BodyForce
{

class GravityForce: public dealii::Function<DIM>
{
public:
    GravityForce(const double rho = 7700,
                 const dealii::Point<DIM> &g = dealii::Point<DIM>(0.0, 0.0, -9.81));
    virtual void
    vector_value(const dealii::Point<DIM> &p, dealii::Vector<double> &values) const override;

    virtual void
    vector_value_list(const std::vector<dealii::Point<DIM>> &points,
                      std::vector<dealii::Vector<double>> &value_list) const override;
private:
    const double rho;
    const dealii::Point<DIM> g;
};

}

#endif // BODYFORCE_HPP
