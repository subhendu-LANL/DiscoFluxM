// Subhendu Chakraborty, Abigail Hunter, Darby J. Luscher
// Funding: LDRD project XX9A, XXNA
// contacts: 
// Subhendu Chakraborty(schakraborty@lanl.gov) 
// Abigail Hunter(ahunter@lanl.gov)
// Darby J. Luscher(djl@lanl.gov)
// Unit system: mm-MPa-s 
/*----------
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for 
Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC 
for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the 
U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, 
irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to 
permit others to do so.
------------*/
#include "ArrayMassLumpedTimeDerivative.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"

registerMooseObject("MooseApp", ArrayMassLumpedTimeDerivative);

InputParameters
ArrayMassLumpedTimeDerivative::validParams()
{
  InputParameters params = ArrayTimeKernel::validParams();
  params.addClassDescription("Array time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<MaterialPropertyName>("time_derivative_coefficient",
                                        "The name of the time derivative coefficient. "
                                        "Can be scalar, vector, or matrix material property.");
  return params;
}

ArrayMassLumpedTimeDerivative::ArrayMassLumpedTimeDerivative(const InputParameters & parameters)
  : ArrayTimeKernel(parameters),
    _has_coefficient(isParamValid("time_derivative_coefficient")),
    _coeff(_has_coefficient && hasMaterialProperty<Real>("time_derivative_coefficient")
               ? &getMaterialProperty<Real>("time_derivative_coefficient")
               : nullptr),
    _coeff_array(_has_coefficient &&
                         hasMaterialProperty<RealEigenVector>("time_derivative_coefficient")
                     ? &getMaterialProperty<RealEigenVector>("time_derivative_coefficient")
                     : nullptr),
    _coeff_2d_array(_has_coefficient &&
                            hasMaterialProperty<RealEigenMatrix>("time_derivative_coefficient")
                        ? &getMaterialProperty<RealEigenMatrix>("time_derivative_coefficient")
                        : nullptr),
  _u_dot_nodal(_var.dofValuesDot())
{
  if (!_coeff && !_coeff_array && !_coeff_2d_array && _has_coefficient)
  {
    MaterialPropertyName mat = getParam<MaterialPropertyName>("time_derivative_coefficient");
    mooseError("Property " + mat + " is of unsupported type for ArrayMassLumpedTimeDerivative");
  }
}

void
ArrayMassLumpedTimeDerivative::computeQpResidual(RealEigenVector & residual)
{
  if (!_has_coefficient)
    residual =  _u_dot_nodal[_i]*_test[_i][_qp];
  else if (_coeff)
    residual = (*_coeff)[_qp] * _test[_i][_qp] * _u_dot_nodal[_i];
  else if (_coeff_array)
  {
    mooseAssert((*_coeff_array)[_qp].size() == _var.count(),
                "time_derivative_coefficient size is inconsistent with the number of components "
                "in array variable");
    residual.noalias() = (*_coeff_array)[_qp].asDiagonal() * _test[_i][_qp] * _u_dot_nodal[_i];
  }
  else
  {
    mooseAssert((*_coeff_2d_array)[_qp].cols() == _var.count(),
                "time_derivative_coefficient size is inconsistent with the number of components "
                "in array variable");
    mooseAssert((*_coeff_2d_array)[_qp].rows() == _var.count(),
                "time_derivative_coefficient size is inconsistent with the number of components "
                "in array variable");
    residual.noalias() = (*_coeff_2d_array)[_qp] * _test[_i][_qp] * _u_dot_nodal[_i];
  }
}

RealEigenVector
ArrayMassLumpedTimeDerivative::computeQpJacobian()
{
  Real tmp = _test[_i][_qp] * _du_dot_du[_qp];
  if (!_has_coefficient)
    return RealEigenVector::Constant(_var.count(), tmp);
  else if (_coeff)
    return RealEigenVector::Constant(_var.count(), tmp * (*_coeff)[_qp]);
  else if (_coeff_array)
    return tmp * (*_coeff_array)[_qp];
  else
    return tmp * (*_coeff_2d_array)[_qp].diagonal();
}


RealEigenMatrix
ArrayMassLumpedTimeDerivative::computeQpOffDiagJacobian(const MooseVariableFEBase & jvar)
{
  if (jvar.number() == _var.number() && _coeff_2d_array)
    return _phi[_j][_qp] * _test[_i][_qp] * _du_dot_du[_qp] * (*_coeff_2d_array)[_qp];
  else
    return ArrayKernel::computeQpOffDiagJacobian(jvar);
}


void
ArrayMassLumpedTimeDerivative::computeOffDiagJacobian(const unsigned int jvar_num)
{
  RealEigenMatrix work_matrix;
  const auto & jvar = getVariable(jvar_num);

  bool same_var = (jvar_num == _var.number());

  prepareMatrixTag(_assembly, _var.number(), jvar_num);

  auto phi_size = jvar.dofIndices().size();

  precalculateOffDiagJacobian(jvar_num);
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    initQpOffDiagJacobian(jvar);
    for (_i = 0; _i < _test.size(); _i++)
      {
        work_matrix = computeQpOffDiagJacobian(jvar) * _JxW[_qp] * _coord[_qp];
        _assembly.saveFullLocalArrayJacobian(
            _local_ke, _i, _test.size(), _i, phi_size, _var.number(), jvar_num, work_matrix);
      }
  }

  accumulateTaggedLocalMatrix();
}
