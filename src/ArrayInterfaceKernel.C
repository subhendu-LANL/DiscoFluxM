//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#include "ArrayInterfaceKernel.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"

template <typename T>
InputParameters
InterfaceKernelTempl02<T>::validParams()
{
  InputParameters params = InterfaceKernelBase::validParams();
  if (std::is_same<T, RealEigenVector>::value)
    params.registerBase("ArrayInterfaceKernel");
  else
    ::mooseError("unsupported InterfaceKernelTempl02 specialization");
  return params;
}

template <typename T>
InterfaceKernelTempl02<T>::InterfaceKernelTempl02(const InputParameters & parameters)
  : InterfaceKernelBase(parameters),
    NeighborMooseVariableInterface<T>(this,
                                      false,
                                      Moose::VarKindType::VAR_NONLINEAR,
                                      std::is_same<T, Real>::value
                                          ? Moose::VarFieldType::VAR_FIELD_STANDARD
                                          : Moose::VarFieldType::VAR_FIELD_ARRAY),
    _var(*this->mooseVariable()),
    _normals(_assembly.normals()),
    _u(_is_implicit ? _var.sln() : _var.slnOld()),
    _grad_u(_is_implicit ? _var.gradSln() : _var.gradSlnOld()),
    _phi(_assembly.phiFace(_var)),
    _grad_phi(_assembly.gradPhiFace(_var)),
    _test(_var.phiFace()),
    _grad_test(_var.gradPhiFace()),
    _neighbor_var(*getVarHelper<MooseVariableFE<T>>("neighbor_var", 0)),
    _neighbor_value(_is_implicit ? _neighbor_var.slnNeighbor() : _neighbor_var.slnOldNeighbor()),
    _grad_neighbor_value(_neighbor_var.gradSlnNeighbor()),
    _phi_neighbor(_assembly.phiFaceNeighbor(_neighbor_var)),
    _grad_phi_neighbor(_assembly.gradPhiFaceNeighbor(_neighbor_var)),
    _test_neighbor(_neighbor_var.phiFaceNeighbor()),
    _grad_test_neighbor(_neighbor_var.gradPhiFaceNeighbor()),
	_count(_var.count()),
	_work_vector(_count)

{
  addMooseVariableDependency(this->mooseVariable());

  if (!parameters.isParamValid("boundary"))
    mooseError(
        "In order to use an interface kernel, you must specify a boundary where it will live.");

  if (parameters.isParamSetByUser("save_in"))
  {
    if (_save_in_strings.size() != _save_in_var_side.size())
      mooseError("save_in and save_in_var_side must be the same length");
    else
    {
      for (unsigned i = 0; i < _save_in_strings.size(); ++i)
      {
        ArrayMooseVariable * var = &_subproblem.getArrayVariable(_tid, _save_in_strings[i]);

        if (_sys.hasVariable(_save_in_strings[i]))
          mooseError("Trying to use solution variable " + _save_in_strings[i] +
                     " as a save_in variable in " + name());

        if (_save_in_var_side[i] == "m")
        {
          if (var->feType() != _var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the primary side nonlinear "
                "variable this interface kernel object is acting on.");
          _primary_save_in_residual_variables.push_back(var);
        }
        else
        {
          if (var->feType() != _neighbor_var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the secondary side nonlinear "
                "variable this interface kernel object is acting on.");
          _secondary_save_in_residual_variables.push_back(var);
        }

        var->sys().addVariableToZeroOnResidual(_save_in_strings[i]);
        addMooseVariableDependency(var);
      }
    }
  }

  _has_primary_residuals_saved_in = _primary_save_in_residual_variables.size() > 0;
  _has_secondary_residuals_saved_in = _secondary_save_in_residual_variables.size() > 0;

  if (parameters.isParamSetByUser("diag_save_in"))
  {
    if (_diag_save_in_strings.size() != _diag_save_in_var_side.size())
      mooseError("diag_save_in and diag_save_in_var_side must be the same length");
    else
    {
      for (unsigned i = 0; i < _diag_save_in_strings.size(); ++i)
      {
        ArrayMooseVariable * var = &_subproblem.getArrayVariable(_tid, _diag_save_in_strings[i]);

        if (_sys.hasVariable(_diag_save_in_strings[i]))
          mooseError("Trying to use solution variable " + _diag_save_in_strings[i] +
                     " as a save_in variable in " + name());

        if (_diag_save_in_var_side[i] == "m")
        {
          if (var->feType() != _var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the primary side nonlinear "
                "variable this interface kernel object is acting on.");
          _primary_save_in_jacobian_variables.push_back(var);
        }
        else
        {
          if (var->feType() != _neighbor_var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the secondary side nonlinear "
                "variable this interface kernel object is acting on.");
          _secondary_save_in_jacobian_variables.push_back(var);
        }

        var->sys().addVariableToZeroOnJacobian(_diag_save_in_strings[i]);
        addMooseVariableDependency(var);
      }
    }
  }

  _has_primary_jacobians_saved_in = _primary_save_in_jacobian_variables.size() > 0;
  _has_secondary_jacobians_saved_in = _secondary_save_in_jacobian_variables.size() > 0;
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeElemNeighResidual(Moose::DGResidualType type)
{
  bool is_elem;
  if (type == Moose::Element)
    is_elem = true;
  else
    is_elem = false;

  const TemplateVariableTestValue & test_space = is_elem ? _test : _test_neighbor;

  if (is_elem)
    prepareVectorTag(_assembly, _var.number());
  else
    prepareVectorTagNeighbor(_assembly, _neighbor_var.number());


  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    initQpResidual(type);
    for (_i = 0; _i < test_space.size(); _i++)
	{
		_work_vector.setZero();
		_work_vector = computeQpResidual(type);
		mooseAssert(_work_vector.size() == _count,
                  "Size of local residual is not equal to the number of array variable compoments");
		_work_vector *= _JxW[_qp] * _coord[_qp];
		_assembly.saveLocalArrayResidual(_local_re, _i, test_space.size(), _work_vector);
	}
  }
	
  accumulateTaggedLocalResidual();

  if (_has_primary_residuals_saved_in && is_elem)
  {
    Threads::spin_mutex::scoped_lock lock(_resid_vars_mutex);
    for (const auto & var : _primary_save_in_residual_variables)
    {
      auto * avar = dynamic_cast<ArrayMooseVariable *>(var);
      if (avar)
        avar->addSolution(_local_re);
      else
        mooseError("Save-in variable of Residual for an array kernel must be an array variable");
    }
  }
  else if (_has_secondary_residuals_saved_in && !is_elem)
  {
    Threads::spin_mutex::scoped_lock lock(_resid_vars_mutex);
    for (const auto & var : _secondary_save_in_residual_variables)
    {
      auto * avar = dynamic_cast<ArrayMooseVariable *>(var);
      if (avar)
        avar->addSolutionNeighbor(_local_re);
      else
        mooseError("Save-in variable of neighbor Residual for an array kernel must be an array variable");
    }
  }
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeResidual()
{
  // in the gmsh mesh format (at least in the version 2 format) the "sideset" physical entities are
  // associated only with the lower-dimensional geometric entity that is the boundary between two
  // higher-dimensional element faces. It does not have a sidedness to it like the exodus format
  // does. Consequently we may naively try to execute an interface kernel twice, one time where _var
  // has dofs on _current_elem *AND* _neighbor_var has dofs on _neighbor_elem, and the other time
  // where _var has dofs on _neighbor_elem and _neighbor_var has dofs on _current_elem. We only want
  // to execute in the former case. In the future we should remove this and add some kind of "block"
  // awareness to interface kernels to avoid all the unnecessary reinit that happens before we hit
  // this return  
  if (!_var.activeOnSubdomain(_current_elem->subdomain_id()) ||
      !_neighbor_var.activeOnSubdomain(_neighbor_elem->subdomain_id()))
    return;
    
  precalculateResidual();
  // Compute the residual for this element
  computeElemNeighResidual(Moose::Element);

  // Compute the residual for the neighbor
  computeElemNeighResidual(Moose::Neighbor);
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeElemNeighJacobian(Moose::DGJacobianType type)
{
  const TemplateVariableTestValue & test_space =
      (type == Moose::ElementElement || type == Moose::ElementNeighbor) ? _test : _test_neighbor;
  const TemplateVariableTestValue & loc_phi =
      (type == Moose::ElementElement || type == Moose::NeighborElement) ? _phi : _phi_neighbor;

  unsigned int ivar, jvar;

  switch (type)
  {
    case Moose::ElementElement:
      ivar = jvar = _var.number();
      break;
    case Moose::ElementNeighbor:
      ivar = _var.number(), jvar = _neighbor_var.number();
      break;
    case Moose::NeighborElement:
      ivar = _neighbor_var.number(), jvar = _var.number();
      break;
    case Moose::NeighborNeighbor:
      ivar = _neighbor_var.number(), jvar = _neighbor_var.number();
      break;
    default:
      mooseError("Unknown DGJacobianType ", type);
  }

  if (type == Moose::ElementElement)
    prepareMatrixTag(_assembly, ivar, jvar);
  else
    prepareMatrixTagNeighbor(_assembly, ivar, jvar, type);

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    for (_i = 0; _i < test_space.size(); _i++)
      for (_j = 0; _j < loc_phi.size(); _j++)
		{
        _work_vector = computeQpJacobian(type) * _JxW[_qp] * _coord[_qp];
        _assembly.saveDiagLocalArrayJacobian(
            _local_ke, _i, test_space.size(), _j, loc_phi.size(), _var.number(), _work_vector);
		}

  accumulateTaggedLocalMatrix();

  if (_has_primary_jacobians_saved_in && type == Moose::ElementElement)
  {
    DenseVector<Number> diag = _assembly.getJacobianDiagonal(_local_ke);
    Threads::spin_mutex::scoped_lock lock(_jacoby_vars_mutex);
    for (const auto & var : _primary_save_in_jacobian_variables)
    {
      auto * avar = dynamic_cast<ArrayMooseVariable *>(var);
      if (avar)
        avar->addSolution(diag);
      else
        mooseError("Save-in variable of Jacobian for an array kernel must be an array variable");
    }
  }
  else if (_has_secondary_jacobians_saved_in && type == Moose::NeighborNeighbor)
  {
    DenseVector<Number> diag = _assembly.getJacobianDiagonal(_local_ke);
    Threads::spin_mutex::scoped_lock lock(_jacoby_vars_mutex);
    for (const auto & var : _secondary_save_in_jacobian_variables)
    {
      auto * avar = dynamic_cast<ArrayMooseVariable *>(var);
      if (avar)
        avar->addSolution(diag);
      else
        mooseError("Save-in variable of neighbor Jacobian for an array kernel must be an array variable");
    }
  }
}


template <typename T>
void
InterfaceKernelTempl02<T>::computeJacobian()
{
  // in the gmsh mesh format (at least in the version 2 format) the "sideset" physical entities are
  // associated only with the lower-dimensional geometric entity that is the boundary between two
  // higher-dimensional element faces. It does not have a sidedness to it like the exodus format
  // does. Consequently we may naively try to execute an interface kernel twice, one time where _var
  // has dofs on _current_elem *AND* _neighbor_var has dofs on _neighbor_elem, and the other time
  // where _var has dofs on _neighbor_elem and _neighbor_var has dofs on _current_elem. We only want
  // to execute in the former case. In the future we should remove this and add some kind of "block"
  // awareness to interface kernels to avoid all the unnecessary reinit that happens before we hit
  // this return
  if (!_var.activeOnSubdomain(_current_elem->subdomain_id()) ||
      !_neighbor_var.activeOnSubdomain(_neighbor_elem->subdomain_id()))
    return;

  computeElemNeighJacobian(Moose::ElementElement);
  computeElemNeighJacobian(Moose::NeighborNeighbor);
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeOffDiagElemNeighJacobian(Moose::DGJacobianType type,
                                                         unsigned int jvar)
{
  const TemplateVariableTestValue & test_space =
      (type == Moose::ElementElement || type == Moose::ElementNeighbor) ? _test : _test_neighbor;
  const TemplateVariableTestValue & loc_phi =
      (type == Moose::ElementElement || type == Moose::NeighborElement) ? _phi : _phi_neighbor;

  unsigned int ivar;

  if (type == Moose::ElementElement || type == Moose::ElementNeighbor)
    ivar = _var.number();
  else
    ivar = _neighbor_var.number();

  if (type == Moose::ElementElement)
    prepareMatrixTag(_assembly, ivar, jvar);
  else
    prepareMatrixTagNeighbor(_assembly, ivar, jvar, type);

  // Prevent calling of Jacobian computation if jvar doesn't lie in the current block
  if ((_local_ke.m() == test_space.size()) && (_local_ke.n() == loc_phi.size()))
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      for (_i = 0; _i < test_space.size(); _i++)
        for (_j = 0; _j < loc_phi.size(); _j++)
          _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpOffDiagJacobian(type, jvar);

  accumulateTaggedLocalMatrix();
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeElementOffDiagJacobian(unsigned int jvar)
{
  // in the gmsh mesh format (at least in the version 2 format) the "sideset" physical entities are
  // associated only with the lower-dimensional geometric entity that is the boundary between two
  // higher-dimensional element faces. It does not have a sidedness to it like the exodus format
  // does. Consequently we may naively try to execute an interface kernel twice, one time where _var
  // has dofs on _current_elem *AND* _neighbor_var has dofs on _neighbor_elem, and the other time
  // where _var has dofs on _neighbor_elem and _neighbor_var has dofs on _current_elem. We only want
  // to execute in the former case. In the future we should remove this and add some kind of "block"
  // awareness to interface kernels to avoid all the unnecessary reinit that happens before we hit
  // this return
  if (!_var.activeOnSubdomain(_current_elem->subdomain_id()) ||
      !_neighbor_var.activeOnSubdomain(_neighbor_elem->subdomain_id()))
    return;

  bool is_jvar_not_interface_var = true;
  if (jvar == _var.number())
  {
    computeElemNeighJacobian(Moose::ElementElement);
    is_jvar_not_interface_var = false;
  }
  if (jvar == _neighbor_var.number())
  {
    computeElemNeighJacobian(Moose::ElementNeighbor);
    is_jvar_not_interface_var = false;
  }

  if (is_jvar_not_interface_var)
  {
    computeOffDiagElemNeighJacobian(Moose::ElementElement, jvar);
    computeOffDiagElemNeighJacobian(Moose::ElementNeighbor, jvar);
  }
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeNeighborOffDiagJacobian(unsigned int jvar)
{
  // in the gmsh mesh format (at least in the version 2 format) the "sideset" physical entities are
  // associated only with the lower-dimensional geometric entity that is the boundary between two
  // higher-dimensional element faces. It does not have a sidedness to it like the exodus format
  // does. Consequently we may naively try to execute an interface kernel twice, one time where _var
  // has dofs on _current_elem *AND* _neighbor_var has dofs on _neighbor_elem, and the other time
  // where _var has dofs on _neighbor_elem and _neighbor_var has dofs on _current_elem. We only want
  // to execute in the former case. In the future we should remove this and add some kind of "block"
  // awareness to interface kernels to avoid all the unnecessary reinit that happens before we hit
  // this return
  if (!_var.activeOnSubdomain(_current_elem->subdomain_id()) ||
      !_neighbor_var.activeOnSubdomain(_neighbor_elem->subdomain_id()))
    return;

  bool is_jvar_not_interface_var = true;
  if (jvar == _var.number())
  {
    computeElemNeighJacobian(Moose::NeighborElement);
    is_jvar_not_interface_var = false;
  }
  if (jvar == _neighbor_var.number())
  {
    computeElemNeighJacobian(Moose::NeighborNeighbor);
    is_jvar_not_interface_var = false;
  }

  if (is_jvar_not_interface_var)
  {
    computeOffDiagElemNeighJacobian(Moose::NeighborElement, jvar);
    computeOffDiagElemNeighJacobian(Moose::NeighborNeighbor, jvar);
  }
}

template class InterfaceKernelTempl02<RealEigenVector>;
