//
// Created by ivan on 10.08.16.
//

#ifndef SOLVERS_EVALUATIONBASE_HPP
#define SOLVERS_EVALUATIONBASE_HPP

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>

namespace Postprocessor
{
//begin namespace Postprocessor
using namespace dealii;

template <int dim>
class EvaluationBase
{
public:
    virtual ~EvaluationBase ()
    {}

    virtual void operator() (const DoFHandler<dim> &doFHandler,
                             const Vector<double> &solution) const = 0;
};



//end namespace Postprocessor
}


#endif //SOLVERS_EVALUATIONBASE_HPP
