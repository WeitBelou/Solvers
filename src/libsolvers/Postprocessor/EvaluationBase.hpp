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

template <int dim>
class OutputResults: public EvaluationBase<dim>
{
public:
    OutputResults (const std::string &outputBaseName,
                   const DataOutBase::OutputFormat outputFormat,
                   const std::string &suffix = "");

    virtual void operator() (const DoFHandler<dim> &doFHandler,
                             const Vector<double> &solution) const override;
private:
    const std::string outputBaseName;
    const std::string suffix;

    const DataOutBase::OutputFormat outputFormat;
};

template <int dim>
OutputResults<dim>::OutputResults (const std::string &outputBaseName,
                                   const DataOutBase::OutputFormat outputFormat,
                                   const std::string &suffix)
    :
    outputBaseName (outputBaseName),
    suffix (suffix),
    outputFormat (outputFormat)
{

}

template <int dim>
void OutputResults<dim>::operator() (const DoFHandler<dim> &doFHandler,
                                     const Vector<double> &solution) const
{
    DataOut <dim> dataOut;
    dataOut.attach_dof_handler (doFHandler);

    dataOut.add_data_vector (solution, "solution");
    dataOut.build_patches ();

    std::string filename = outputBaseName
                           + suffix
                           + dataOut.default_suffix (outputFormat);

    std::ofstream out (filename);
    dataOut.write (out, outputFormat);
}

//end namespace Postprocessor
}


#endif //SOLVERS_EVALUATIONBASE_HPP
