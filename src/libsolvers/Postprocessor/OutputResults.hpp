//
// Created by ivan on 11.08.16.
//

#ifndef SOLVERS_OUTPUTRESULTS_HPP
#define SOLVERS_OUTPUTRESULTS_HPP

#include "EvaluationBase.hpp"
namespace Postprocessor
{
//brgin namespace Postprocessor
template <int dim>
class OutputResults: public EvaluationBase<dim>
{
public:
    OutputResults (const std::string &outputDir,
                   const std::string &outputBaseName,
                   const DataOutBase::OutputFormat outputFormat,
                   const std::string &suffix = "");

    virtual void operator() (const DoFHandler<dim> &doFHandler,
                             const Vector<double> &solution) const override;
private:
    const std::string outputDir;
    const std::string outputBaseName;
    const std::string suffix;

    const DataOutBase::OutputFormat outputFormat;
};

template <int dim>
OutputResults<dim>::OutputResults (const std::string &outputDir,
                                   const std::string &outputBaseName,
                                   const DataOutBase::OutputFormat outputFormat,
                                   const std::string &suffix)
    :
    outputDir (outputDir),
    outputBaseName (outputBaseName),
    suffix (suffix),
    outputFormat (outputFormat)
{

}

template <int dim>
void OutputResults<dim>::operator() (const DoFHandler<dim> &doFHandler,
                                     const Vector<double> &solution) const
{
    DataOut<dim> dataOut;
    dataOut.attach_dof_handler (doFHandler);

    dataOut.add_data_vector (solution, "solution");
    dataOut.build_patches ();

    std::string filename = outputDir
                           + outputBaseName
                           + suffix
                           + dataOut.default_suffix (outputFormat);

    std::ofstream out (filename);
    dataOut.write (out, outputFormat);
}

//end namespace Posprocessor
}


#endif //SOLVERS_OUTPUTRESULTS_HPP
