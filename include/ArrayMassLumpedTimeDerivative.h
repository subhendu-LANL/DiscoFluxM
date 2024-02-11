// Subhendu Chakraborty, Abigail Hunter, Darby J. Luscher
// Funding: LDRD project XX9A, XXNA
// contacts: 
// Subhendu Chakraborty(schakraborty@lanl.gov) 
// Abigail Hunter(ahunter@lanl.gov)
// Darby J. Luscher(djl@lanl.gov)
// Unit system: mm-MPa-s 

#pragma once

#include "ArrayTimeKernel.h"

class ArrayMassLumpedTimeDerivative : public ArrayTimeKernel
{
public:
  static InputParameters validParams();

  ArrayMassLumpedTimeDerivative(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual RealEigenVector computeQpJacobian() override;
  virtual RealEigenMatrix computeQpOffDiagJacobian(const MooseVariableFEBase & jvar) override;

  /// Computes full Jacobian of jvar and the array variable this kernel operates on
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

  /// whether or not the coefficient property is provided
  const bool _has_coefficient;
  /// scalar time derivative coefficient
  const MaterialProperty<Real> * const _coeff;
  /// array time derivative coefficient
  const MaterialProperty<RealEigenVector> * const _coeff_array;
  /// matrix time derivative coefficient
  const MaterialProperty<RealEigenMatrix> * const _coeff_2d_array;

  const ArrayVariableValue & _u_dot_nodal;
};
