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
#include "CrystalPlasticityDislocationEdgeScrew.h"
#include "RankTwoTensor.h"

#include "SystemBase.h"
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

#include "Assembly.h" 
#include "MooseMesh.h"

#include "libmesh/node.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/remote_elem.h"

registerMooseObject("TensorMechanicsApp", CrystalPlasticityDislocationEdgeScrew);

InputParameters
CrystalPlasticityDislocationEdgeScrew::validParams()
{
  InputParameters params = CrystalPlasticityStressUpdateBase::validParams();
  
  params.addClassDescription("Calculates the updated stress based on DiscoFlux crystal plasticity material model.");
  
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");
  params.addParam<Real>("lattice_friction", 10, "initial lattice friction strength of the material in MPa");  
  params.addParam<Real>("burgers_vector_mag",1.0e-07,"Magnitude of the Burgers vector in mm");
  params.addParam<Real>("dislo_max_velocity",1.0e-02,"Maximum dislocation velocity (phonon drag) in mm/s");
  params.addParam<Real>("dislo_density_initial",1.0e+8,"Initial immobile dislocation density");
  params.addParam<Real>("dislo_density_factor_CDT",1.0,"factor to convert the dislocation density from CDT to CP, this is for scaling of solution variable");
  
  params.addParam<Real>("C_multi", 8.96e-05, "parameter for dislocation multiplication");
  params.addParam<Real>("C_trap", 3.0e-03, "parameter for dislocation trapping");
  params.addParam<Real>("C_m_ann", 0.5, "parameter for dislocation mobile annihilation");
  params.addParam<Real>("C_im_ann", 0.5, "parameter for dislocation immobile annihilation");
  params.addParam<Real>("Coeff_hardening", 0.5, "parameter to control the material hardening");
  
  params.addCoupledVar("DD_EdgePositive", 1.0e+5, "Coupled dislocation density, EdgePositive");
  params.addCoupledVar("DD_EdgeNegative", 1.0e+5, "Coupled dislocation density, EdgeNegative");
  params.addCoupledVar("DD_ScrewPositive", 1.0e+5, "Coupled dislocation density, ScrewPositive");
  params.addCoupledVar("DD_ScrewNegative", 1.0e+5, "Coupled dislocation density, ScrewNegative");

  return params;
}

CrystalPlasticityDislocationEdgeScrew::CrystalPlasticityDislocationEdgeScrew(
    const InputParameters & parameters)
  : CrystalPlasticityStressUpdateBase(parameters),
    _xm(getParam<Real>("xm")),
	_lattice_friction(getParam<Real>("lattice_friction")),
	
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_dislo_max_velocity(getParam<Real>("dislo_max_velocity")),
	_dislo_density_initial(getParam<Real>("dislo_density_initial")),
	_dislo_density_factor_CDT(getParam<Real>("dislo_density_factor_CDT")),
	_C_multi(getParam<Real>("C_multi")),
	_C_trap(getParam<Real>("C_trap")),
	_C_m_ann(getParam<Real>("C_m_ann")),
	_C_im_ann(getParam<Real>("C_im_ann")),
	_Coeff_hardening(getParam<Real>("Coeff_hardening")),
	
	//_Rho_EdgePositive_01(coupledValue("Rho_EdgePositive_01")),
	_DD_EdgePositive(coupledArrayValue("DD_EdgePositive")), 
	_DD_EdgeNegative(coupledArrayValue("DD_EdgeNegative")), 
	_DD_ScrewPositive(coupledArrayValue("DD_ScrewPositive")), 
	_DD_ScrewNegative(coupledArrayValue("DD_ScrewNegative")),
	
	_DD_EdgePositive_Grad(coupledArrayGradient("DD_EdgePositive")),
	_DD_EdgeNegative_Grad(coupledArrayGradient("DD_EdgeNegative")),
	_DD_ScrewPositive_Grad(coupledArrayGradient("DD_ScrewPositive")),
	_DD_ScrewNegative_Grad(coupledArrayGradient("DD_ScrewNegative")),
	
	_slip_direction_edge(declareProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_direction_screw(declareProperty<std::vector<RealVectorValue>>("slip_direction_screw")),
	_slip_plane_normalboth(declareProperty<std::vector<RealVectorValue>>("slip_plane_normalboth")),
	_crysrot(getMaterialPropertyOld<RankTwoTensor>(_base_name + "crysrot")),
	
	_dislocation_immobile(declareProperty<std::vector<Real>>(_base_name + "dislocation_immobile")),
	_dislocation_immobile_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "dislocation_immobile")),
	_dislocation_immobile_total(declareProperty<Real>(_base_name + "dislocation_immobile_total")),
	
	_dislo_velocity_edge(declareProperty<std::vector<Real>>("dislo_velocity_edge")),
	_dislo_velocity_screw(declareProperty<std::vector<Real>>("dislo_velocity_screw")),
		
	// DDC related variables
	_kappa_grad(_number_slip_systems, 0.00),
	_tau_b(_number_slip_systems, 0.00),

  // resize local caching vectors used for substepping
  _previous_substep_slip_resistance(_number_slip_systems, 0.00),
  _previous_substep_dislocation_mobile(_number_slip_systems, 0.00),
  _previous_substep_dislocation_immobile(_number_slip_systems, 0.00),
  
  _slip_resistance_before_update(_number_slip_systems, 0.00),
  _dislocation_mobile_before_update(_number_slip_systems, 0.00),
  _dislocation_immobile_before_update(_number_slip_systems, 0.00),

  // resize vectors used in the consititutive slip hardening
  _hb(_number_slip_systems, 0.00),
  _slip_resistance_increment(_number_slip_systems, 0.00),
  _dislocation_mobile_increment(_number_slip_systems, 0.00),
  _dislocation_immobile_increment(_number_slip_systems, 0.00),
  _dv_dtau(_number_slip_systems, 0.00),
  _L_bar(_number_slip_systems, 0.00),

  // resize local variables realted to gerDisoVelocity**
  t_wait(_number_slip_systems,0.00), // dislocation waiting time at barrier
  t_run(_number_slip_systems,0.00), // running time for dislocations between barrier
  vel_run(_number_slip_systems,0.00),
  dislocation_density(_number_slip_systems,0.00),
  tau_b(_number_slip_systems,0.00), // back stress, not evolving for time being.
  xi0(_number_slip_systems,0.00),
  tau_eff(_number_slip_systems,0.00),
  tau_effAbs(_number_slip_systems,0.00),
  tau_effSign(_number_slip_systems,0.00),
  slip_r(_number_slip_systems,0.00)
{
  
}

void
CrystalPlasticityDislocationEdgeScrew::initQpStatefulProperties()
{
  CrystalPlasticityStressUpdateBase::initQpStatefulProperties();
    
  _slip_direction_edge[_qp].resize(_number_slip_systems);
  _slip_direction_screw[_qp].resize(_number_slip_systems);
  _slip_plane_normalboth[_qp].resize(_number_slip_systems);


  _dislocation_immobile[_qp].resize(_number_slip_systems);
  _dislo_velocity_edge[_qp].resize(_number_slip_systems);
  _dislo_velocity_screw[_qp].resize(_number_slip_systems);
  
  _dislocation_immobile_total[_qp] = 0.00;
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {  
    _slip_direction_edge[_qp][i].zero();
	_slip_direction_screw[_qp][i].zero();
	_slip_plane_normalboth[_qp][i].zero();

    _slip_resistance[_qp][i] = _lattice_friction + _Coeff_hardening*mu*_burgers_vector_mag*std::sqrt(1.0*4*2 * _dislo_density_initial + 0.1*11*4*2*_dislo_density_initial); // approximate
    _slip_increment[_qp][i] = 0.0;
	
	_dislocation_immobile[_qp][i] = 4.0 * _dislo_density_initial;
	_dislocation_immobile_total[_qp] += _dislocation_immobile[_qp][i];
	_dislo_velocity_edge[_qp][i] = 0.00;
	_dislo_velocity_screw[_qp][i] = 0.00;
	
	// DDC related stuff
	//_kappa_grad[_qp][i].zero();
  }
  
}

void
CrystalPlasticityDislocationEdgeScrew::storeDislocationMobilityInformation()
{
  for (const auto i : make_range(_number_slip_systems))
  {
	_slip_direction_edge[_qp][i] = _slip_direction[i];
	_slip_direction_edge[_qp][i] /= _slip_direction_edge[_qp][i].norm();
	_slip_direction_edge[_qp][i] = _crysrot[_qp] * _slip_direction_edge[_qp][i]; 
	
	_slip_direction_screw[_qp][i] = _slip_direction[i].cross(_slip_plane_normal[i]);
	_slip_direction_screw[_qp][i] /= _slip_direction_screw[_qp][i].norm();
	_slip_direction_screw[_qp][i] = _crysrot[_qp] * _slip_direction_screw[_qp][i]; 
	
	_slip_plane_normalboth[_qp][i] = _slip_plane_normal[i];
	_slip_plane_normalboth[_qp][i] /= _slip_plane_normalboth[_qp][i].norm();
	_slip_plane_normalboth[_qp][i] = _crysrot[_qp] * _slip_plane_normalboth[_qp][i];
  }

}

void
CrystalPlasticityDislocationEdgeScrew::setInitialConstitutiveVariableValues()
{
  storeDislocationMobilityInformation();
  _slip_resistance[_qp] = _slip_resistance_old[_qp];
  _previous_substep_slip_resistance = _slip_resistance_old[_qp];
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
  _previous_substep_dislocation_mobile[i] = (_DD_EdgePositive[_qp][i] + _DD_EdgeNegative[_qp][i] + _DD_ScrewPositive[_qp][i] + _DD_ScrewNegative[_qp][i]) * _dislo_density_factor_CDT;
  }
	_dislocation_immobile[_qp] = _dislocation_immobile_old[_qp];
    _previous_substep_dislocation_immobile = _dislocation_immobile_old[_qp];
}

void
CrystalPlasticityDislocationEdgeScrew::setSubstepConstitutiveVariableValues()
{
  _slip_resistance[_qp] = _previous_substep_slip_resistance;
  _dislocation_immobile[_qp] = _previous_substep_dislocation_immobile;
}

bool
CrystalPlasticityDislocationEdgeScrew::calculateSlipRate()
{	
	// compute dislocation velocity according to DiscoFlux material model
    DDCUpdate();
    getDisloVelocity();

    for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
		_slip_increment[_qp][i] = (_DD_EdgePositive[_qp][i] + _DD_EdgeNegative[_qp][i]) * _dislo_density_factor_CDT * _burgers_vector_mag * _dislo_velocity_edge[_qp][i] +
								( _DD_ScrewPositive[_qp][i] + _DD_ScrewNegative[_qp][i]) * _dislo_density_factor_CDT * _burgers_vector_mag * _dislo_velocity_screw[_qp][i];

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
		return false;
    }
  }
  return true;
}

void
CrystalPlasticityDislocationEdgeScrew::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0;
    else
	{
	    dslip_dtau[i] = (_DD_EdgePositive[_qp][i] + _DD_EdgeNegative[_qp][i]) * _dislo_density_factor_CDT * _burgers_vector_mag * _dv_dtau[i];
	}	
 }
}

void
CrystalPlasticityDislocationEdgeScrew::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];
  _dislocation_immobile_before_update = _dislocation_immobile[_qp];
}

void
CrystalPlasticityDislocationEdgeScrew::calculateStateVariableEvolutionRateComponent()
{ 
  // calculate dislocation density increment
  getDDIncrements();
}

bool
CrystalPlasticityDislocationEdgeScrew::updateStateVariables()
{
	Real Hij, eff_dislocation_density = 0.00;

     _dislocation_immobile_total[_qp] = 0.0;
     for (unsigned int i = 0; i < _number_slip_systems; ++i)
      {  
	    _dislocation_immobile[_qp][i] = _previous_substep_dislocation_immobile[i] + _dislocation_immobile_increment[i] * _substep_dt;
	  
	    _dislocation_immobile_total[_qp] += _dislocation_immobile[_qp][i];
      }

	  for (unsigned int i = 0; i < _number_slip_systems; ++i)
	  {
		  eff_dislocation_density = 0.00;
		  for (unsigned int j = 0; j < _number_slip_systems; ++j)
		  {
			if(i==j ? Hij=1.0: Hij=0.1) 
		    eff_dislocation_density += Hij*(_previous_substep_dislocation_mobile[j] + _dislocation_immobile[_qp][j]); 
		  }
	  _slip_resistance[_qp][i] = _lattice_friction + _Coeff_hardening*mu*_burgers_vector_mag*std::sqrt(eff_dislocation_density); 
	  }

  return true;
}

bool
CrystalPlasticityDislocationEdgeScrew::areConstitutiveStateVariablesConverged()
{
	bool flagSlipResistanceConverged, flagDislocationDensityConverged;
	Real _DislocationDensity_tol = 0.1;

// Check convergence of _slip_resistance
   flagSlipResistanceConverged = isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                              _slip_resistance_before_update,
                                              _previous_substep_slip_resistance,
                                              _resistance_tol);

// Check convergence of _dislocation_immobile
   flagDislocationDensityConverged = isConstitutiveStateVariableConverged(_dislocation_immobile[_qp],
                                              _dislocation_immobile_before_update,
                                              _previous_substep_dislocation_immobile,
                                              _DislocationDensity_tol);
	if(!flagDislocationDensityConverged)
	{
	  return false;
	}

	return true;
}

void
CrystalPlasticityDislocationEdgeScrew::updateSubstepConstitutiveVariableValues()
{
  _previous_substep_slip_resistance = _slip_resistance[_qp];
  _previous_substep_dislocation_mobile = _dislocation_mobile[_qp];
  _previous_substep_dislocation_immobile = _dislocation_immobile[_qp];
}

// Calculate Dislocation Density increment
void
CrystalPlasticityDislocationEdgeScrew::getDDIncrements()
{  
  Real small2 = 1.0e-5;

  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {  
	_dislocation_mobile_increment[i] = 0.00;
	_dislocation_immobile_increment[i] = 0.00;
	
	if (std::abs(_slip_increment[_qp][i]) > small2) {
	_dislocation_immobile_increment[i] = (_C_trap/_burgers_vector_mag)*std::pow(_dislocation_mobile[_qp][i],0.5)* std::abs(_slip_increment[_qp][i])
		- (_C_im_ann*_dislocation_immobile[_qp][i])* std::abs(_slip_increment[_qp][i]);	
	 
	} 
	else
	 {
		
	  _dislocation_mobile_increment[i] = 0.0;
	  _dislocation_immobile_increment[i] = 0.0;
	  
	 }
 }
}

void
CrystalPlasticityDislocationEdgeScrew::getDisloVelocity()
{ 
	//  compute velocity for each slip system
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
   {
	   tau_eff[i] = 0.00;
		t_wait[i]  = 0.00;
		t_run[i]   = 0.00;
		vel_run[i] = 0.00;
		_dislo_velocity_edge[_qp][i] =0.00;
		_dislo_velocity_screw[_qp][i] =0.00;
		_dv_dtau[i] = 0.00;
   }
   
   for (unsigned int i = 0; i < _number_slip_systems; ++i)
   {
	   slip_r[i]  = _slip_resistance[_qp][i];

	   tau_eff[i] = (_tau[_qp][i] - _tau_b[i]); 
	   tau_effAbs[i] = std::abs(tau_eff[i]);
	   tau_effSign[i] = std::copysign(1.0, tau_eff[i]);
   }
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  { 
	  if (tau_effAbs[i]>small2)
	  {
		  xi0[i] = B0*vcrit / (2*_burgers_vector_mag*tau_effAbs[i]);
		  vel_run[i] = vcrit*(std::pow((std::pow(xi0[i],2)+1),0.5) - xi0[i]);
	  }
	  else
	  vel_run[i] = 0.00;
  }
   deltaG0 = g0*mu*std::pow(_burgers_vector_mag,3)*1.0e-3;
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {	   
	  inner = 1.0 - std::pow((tau_effAbs[i] / slip_r[i] ),q1);
	  if (inner > 0.00)
	  {
		  deltaG  = deltaG0*( std::pow(inner,q2) );
		  exp_arg = deltaG / (boltz*temp);
			t_wait[i] = (exp(exp_arg) - 1.0) / omega0;
	  }
	  else
		 t_wait[i] = 0.00;
  }
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	if (vel_run[i] > small2)
	   {  
		t_run[i] = _L_bar[i] / vel_run[i];
		_dislo_velocity_edge[_qp][i] = tau_effSign[i]*_L_bar[i]/ (t_wait[i] + t_run[i]);

		inner = 1.0 - std::pow((tau_effAbs[i] / slip_r[i] ),q1);
		deltaG  = deltaG0*( std::pow(inner,q2) );
		exp_arg = deltaG / (boltz*temp);
		dtw_dtau = (-1)*q1*q2*deltaG0/(omega0*boltz*temp) * exp(exp_arg) * ( std::pow(inner,q2-1) ) * std::pow((tau_effAbs[i] / slip_r[i] ),q1-1) *(tau_effSign[i]/slip_r[i]);
		dtr_dtau = (_L_bar[i]*B0/_burgers_vector_mag)*(tau_effSign[i]/std::pow(tau_effAbs[i],2));
		_dv_dtau[i] = -1*_L_bar[i]*std::pow((t_wait[i] + t_run[i]),-2) * (dtw_dtau + dtr_dtau);
		}
	else
	  {
	  _dislo_velocity_edge[_qp][i] = 0.00;
	  _dv_dtau[i] = 0.00;
	  }
  }

}


//----Compute the back stress based on Dislocation Deformation Compatibility(DDC)
void
CrystalPlasticityDislocationEdgeScrew::DDCUpdate()
{	

	Real dislocationsPositive, dislocationsNegative;
  
		for (unsigned int i = 0; i < _number_slip_systems; ++i)
		  {	
			slip_direction_rotated = _crysrot[_qp] * _slip_direction[i];
			slip_plane_normal_rotated = _crysrot[_qp] * _slip_plane_normal[i];
			_L_bar[i] = std::pow((_dislocation_mobile[_qp][i] + _dislocation_immobile[_qp][i]),-0.5); 
			dislocationsPositive = _DD_EdgePositive[_qp][i];
			dislocationsPositive = _DD_EdgeNegative[_qp][i];

			_kappa_grad[i](0) = (_DD_EdgePositive_Grad[_qp](i) - _DD_EdgeNegative_Grad[_qp](i))*_dislo_density_factor_CDT;
			_kappa_grad[i](1) = (_DD_EdgePositive_Grad[_qp](i+_number_slip_systems) - _DD_EdgeNegative_Grad[_qp](i+_number_slip_systems))*_dislo_density_factor_CDT;
			_kappa_grad[i](2) = (_DD_EdgePositive_Grad[_qp](i+2*_number_slip_systems) - _DD_EdgeNegative_Grad[_qp](i+2*_number_slip_systems))*_dislo_density_factor_CDT;
			_tau_b[i] = (( mu * std::pow(_L_bar[i],1))/(2*3.141*(1-nu)))*_burgers_vector_mag * (_kappa_grad[i]*slip_direction_rotated);
		}	
}