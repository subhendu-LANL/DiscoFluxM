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
#include "ArrayAdvectionCDT.h"

registerMooseObject("TensorMechanicsApp", ArrayAdvectionCDT);

InputParameters
ArrayAdvectionCDT::validParams()
{
  InputParameters params = ArrayKernel::validParams();
  params.addClassDescription("Continuum transport of dislocations(array variable) is modeled using advection model.");
  params.addRequiredParam<MaterialPropertyName>("advection_coefficient",
                                                "The name of the advdectivity, "
                                                "can be scalar, vector, or matrix.");
  MooseEnum upwinding_type("none full", "full");
  params.addParam<MooseEnum>("upwinding_type",
                             upwinding_type,
                             "Type of upwinding used.  None: Typically results in overshoots and "
                             "undershoots, but numerical diffusion is minimized.  Full: Overshoots "
                             "and undershoots are avoided, but numerical diffusion is large");
  MooseEnum dislocation_character("edge screw");
  params.addRequiredParam<MooseEnum>("dislocation_character", dislocation_character, "Character of dislocation");
  MooseEnum dislocation_sign("positive negative");
  params.addRequiredParam<MooseEnum>("dislocation_sign", dislocation_sign, "Sign of dislocation");
  params.addClassDescription(
      "The array Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
      "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  return params;
}

ArrayAdvectionCDT::ArrayAdvectionCDT(const InputParameters & parameters)
  : ArrayKernel(parameters),
	_dislo_velocity_CP_edge(getMaterialProperty<std::vector<Real>>("dislo_velocity_edge")),
	_dislo_velocity_CP_screw(getMaterialProperty<std::vector<Real>>("dislo_velocity_screw")),
	
  _upwinding(getParam<MooseEnum>("upwinding_type").getEnum<UpwindingType>()),
	_u_nodal(_var.dofValues()),
  _upwind_node(0),
  _dtotal_mass_out(0),
  work_vector(_count),

	_dislocationcharacter(getParam<MooseEnum>("dislocation_character").getEnum<DislocationCharacter>()),
	_dislocationsign(getParam<MooseEnum>("dislocation_sign").getEnum<DislocationSign>()), 
	
	_slip_direction_edge(getMaterialProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_direction_screw(getMaterialProperty<std::vector<RealVectorValue>>("slip_direction_screw")),
	_slip_plane_normalboth(getMaterialProperty<std::vector<RealVectorValue>>("slip_plane_normalboth"))
	
{
}


RealEigenVector
ArrayAdvectionCDT::negativeVelocityGradTestQp()
{    
	RealVectorValue _slip_direction_rotated;
  _Dislo_Velocity.resize(_var.count(),LIBMESH_DIM);
  _Dislo_Velocity.setZero();
  
  for (unsigned int i = 0; i < _var.count(); ++i) 
	{
    _slip_direction_rotated = _slip_direction_edge[_qp][i]; // already rotated

		switch (_dislocationsign)
		{
			case DislocationSign::negative:
			_slip_direction_rotated *= (-1);
			break;
		}
		for (unsigned int j = 0; j < LIBMESH_DIM; ++j) 
        _Dislo_Velocity(i,j) = _dislo_velocity_CP_edge[_qp][i] * _slip_direction_rotated(j);
  }
	
	return (-1) * _Dislo_Velocity * _array_grad_test[_i][_qp]; 
}

void
ArrayAdvectionCDT::computeQpResidual(RealEigenVector & residual)
{
  residual = (_u[_qp].asDiagonal() * negativeVelocityGradTestQp());
}

RealEigenVector
ArrayAdvectionCDT::computeQpJacobian()
{	
  RealEigenVector jac ;
  jac.resize(_count);
  jac.setZero();

	jac = _phi[_j][_qp] * negativeVelocityGradTestQp();

	return jac;
}

RealEigenMatrix
ArrayAdvectionCDT::computeQpOffDiagJacobian(const MooseVariableFEBase & jvar)
{
    return ArrayKernel::computeQpOffDiagJacobian(jvar);
}


void
ArrayAdvectionCDT::computeResidual()
{
  switch (_upwinding)
  {
    case UpwindingType::none:
      ArrayKernel::computeResidual();
      break;
    case UpwindingType::full:
      fullUpwind(JacRes::CALCULATE_RESIDUAL);
      break;
  }
}

void
ArrayAdvectionCDT::computeOffDiagJacobian(const unsigned int jvar_num)
{
  switch (_upwinding)
  {
    case UpwindingType::none:
      ArrayKernel::computeOffDiagJacobian(jvar_num);
      break;
    case UpwindingType::full:
      jvar_num02 = jvar_num;
      fullUpwind(JacRes::CALCULATE_JACOBIAN);
      break;
  }
}

void
ArrayAdvectionCDT::computeJacobian() // not needed, will be removed
{
  switch (_upwinding)
  {
    case UpwindingType::none:
      ArrayKernel::computeJacobian();
      break;
    case UpwindingType::full:
      ArrayKernel::computeJacobian();
      break;
  }
}

// Stabilization of the advection kernel 
void
ArrayAdvectionCDT::fullUpwind(JacRes res_or_jac)
{

  const unsigned int num_nodes = _test.size();

  prepareVectorTag(_assembly, _var.number());

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    prepareMatrixTag(_assembly, _var.number(), _var.number());

  _upwind_node.resize(_count,std::vector<bool>(num_nodes));
  RealEigenMatrix work_vector_res;
  work_vector_res.resize(_count,num_nodes); work_vector_res.setZero();
  work_vector.setZero(_count);
  RealEigenVector SpeedQp(_count); SpeedQp.setZero();
  RealEigenVector SpeedQpJac(_count); SpeedQpJac.setZero();
  

for (_i = 0; _i < num_nodes; _i++)
  {

   for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
      SpeedQp =negativeVelocityGradTestQp(); 
      for (unsigned int i_var=0; i_var<_count; i_var++)
      work_vector_res(i_var,_i) += _JxW[_qp] * _coord[_qp] * SpeedQp[i_var];
  }

  for (unsigned int i_var=0; i_var<_count; i_var++)
  _upwind_node[i_var][_i] = (work_vector_res(i_var,_i) >= 0.0);
  }

  std::vector<Real> total_mass_out;
  std::vector<Real> total_in;
  total_mass_out.assign(_count, 0.0);
  total_in.assign(_count, 0.0);

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    for (unsigned int i_var=0; i_var<_count; i_var++)

for (unsigned int i_var = 0; i_var < _count; i_var++)
{
  for (unsigned int n = 0; n < num_nodes; n++)
  {
    if (_upwind_node[i_var][n])
    {
      work_vector_res(i_var,n) *= _u_nodal[n][i_var];
      total_mass_out[i_var] += work_vector_res(i_var,n);
    }
    else 
      total_in[i_var] -= work_vector_res(i_var,n);
  }
}

for (unsigned int i_var = 0; i_var < _count; i_var++)
{
  if (std::abs(total_mass_out[i_var]) > 1.0e-20) //only if  velocity is non-zero
  {
for (unsigned int n = 0; n < num_nodes; n++)
  {
    if (!_upwind_node[i_var][n])
    {
      work_vector_res(i_var,n) *= total_mass_out[i_var] / total_in[i_var];
    }
  }
  }
}

  if (res_or_jac == JacRes::CALCULATE_RESIDUAL)
  {
    for (_i = 0; _i < _test.size(); _i++)
    {
      work_vector.setZero(_count);
      work_vector02.setZero(_count);

    for (unsigned int i_var=0; i_var<_count; i_var++)
      work_vector(i_var) = work_vector_res(i_var,_i);
    _assembly.saveLocalArrayResidual(_local_re, _i, _test.size(), work_vector);

    }
  accumulateTaggedLocalResidual();
  }
  
  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
  {
  ArrayKernel::computeOffDiagJacobian(jvar_num02);
  }

}

