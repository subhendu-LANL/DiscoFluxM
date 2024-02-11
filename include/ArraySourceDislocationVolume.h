// Subhendu Chakraborty, Abigail Hunter, Darby J. Luscher
// Funding: LDRD project XX9A, XXNA
// contacts: 
// Subhendu Chakraborty(schakraborty@lanl.gov) 
// Abigail Hunter(ahunter@lanl.gov)
// Darby J. Luscher(djl@lanl.gov)
// UNit system: mm-MPa-s 
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

#pragma once

#include "ArrayKernel.h"

class ArraySourceDislocationVolume : public ArrayKernel
{
public:
  static InputParameters validParams();

  ArraySourceDislocationVolume(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual RealEigenVector computeQpJacobian() override;
  void computeDisloVelocity();
  
  const Real _dislo_density_factor_CDT;
  const Real _C_multi, _C_trap, _C_m_ann, _C_im_ann, _burgers_vector_mag;

  const MaterialProperty<std::vector<Real>> & _dislo_velocity_CP_edge;
  const MaterialProperty<std::vector<Real>> & _dislo_velocity_CP_screw;

  // Dislocation character
  const enum class DislocationCharacter { edge, screw } _dislocationcharacter;  
  // Dislocation sign
  const enum class DislocationSign { positive, negative } _dislocationsign;
  
  const MaterialProperty<std::vector<Real>> & _dislocation_mobile_old;
  const MaterialProperty<std::vector<Real>> & _slip_increment;
  const MaterialProperty<std::vector<Real>> & _slip_increment_old;
  
  Real _dt;

  RealEigenVector dislo_velocity;
  
};
