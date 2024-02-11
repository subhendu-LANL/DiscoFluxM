[Mesh]
 [./MeshFile]
  type = FileMeshGenerator
  file = 'Mesh_Bi-Crystal.e'
 [../]
[]

[Outputs]
    csv = true
    exodus = true
[]

[UserObjects]
  [./prop_read]
    type = PropertyReadFile
    prop_file_name = 'euler_ang_file.inp'
    nprop = 3
    read_type = block
    nblock= 50
  [../]
[]

[Variables]
  [./DD_EdgePositive]
    components = 12
  [../]
  [./DD_EdgeNegative]
    components = 12
  [../]
  [./DD_ScrewPositive]
    components = 12
  [../]
  [./DD_ScrewNegative]
    components = 12
  [../]
[]
  
[ICs]
  [./IC_DD_EdgePositive]
    type = ArrayFunctionIC  
    variable = DD_EdgePositive
    function = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
  [../]  
  [./IC_DD_EdgeNegative]
    type = ArrayFunctionIC  
    variable = DD_EdgeNegative
   function = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
  [../] 
  [./IC_DD_ScrewPositive]
    type = ArrayFunctionIC  
    variable = DD_ScrewPositive
    function = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
  [../]  
  [./IC_DD_ScrewNegative]
    type = ArrayFunctionIC  
    variable = DD_ScrewNegative
    function = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
  [../]
[]

[Kernels]
  [./dot_DD_EdgePositive]
    type = ArrayTimeDerivative  
    variable = DD_EdgePositive
    time_derivative_coefficient = AdvCoeff
  [../]
  [./Advection_DD_EdgePositive]
    type = ArrayAdvectionCDT  
    variable = DD_EdgePositive
    dislocation_character = edge
    dislocation_sign = positive
    advection_coefficient = AdvCoeff
  [../]
  [./Source_DD_EdgePositive]
    type = ArraySourceDislocationVolume  
    variable = DD_EdgePositive
    dislocation_character = edge
    dislocation_sign = positive
  [../]
  
  [./dot_DD_EdgeNegative]
    type = ArrayTimeDerivative  
    variable = DD_EdgeNegative
    time_derivative_coefficient = AdvCoeff
  [../]
  [./Advection_DD_EdgeNegative]
    type = ArrayAdvectionCDT
    variable = DD_EdgeNegative
    dislocation_character = edge
    dislocation_sign = negative
    advection_coefficient = AdvCoeff
  [../]
  [./Source_DD_EdgeNegative]
    type = ArraySourceDislocationVolume  
    variable = DD_EdgeNegative
    dislocation_character = edge
    dislocation_sign = negative
  [../]
  
  [./dot_DD_ScrewPositive]
    type = ArrayTimeDerivative  
    variable = DD_ScrewPositive
    time_derivative_coefficient = AdvCoeff
  [../]
  [./Advection_DD_ScrewPositive]
    type = ArrayAdvectionCDT  
    variable = DD_ScrewPositive
    dislocation_character = screw
    dislocation_sign = positive
    advection_coefficient = AdvCoeff
  [../]
  [./Source_DD_ScrewPositive]
    type = ArraySourceDislocationVolume  
    variable = DD_ScrewPositive
    dislocation_character = screw
    dislocation_sign = positive
  [../]
  
  [./dot_DD_ScrewNegative]
    type = ArrayTimeDerivative  
    variable = DD_ScrewNegative
    time_derivative_coefficient = AdvCoeff
  [../]
  [./Advection_DD_ScrewNegative]
    type = ArrayAdvectionCDT
    variable = DD_ScrewNegative
    dislocation_character = screw
    dislocation_sign = negative
    advection_coefficient = AdvCoeff
  [../]
  [./Source_DD_ScrewNegative]
    type = ArraySourceDislocationVolume  
    variable = DD_ScrewNegative
    dislocation_character = screw
    dislocation_sign = negative
  [../]
[]

[BCs]
  [./BC_Loading]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top bottom '
    function = '0.1*y*t'
  [../]

  [./BC_BottomLeft_X]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'bottom'
    function = 0.0
  [../]
   [./BC_BottomLeft_Y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'bottom'
    function = 0.0
  [../]
   [./BC_BottomBack_Z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'bottom'
    function = 0.0
  [../]

  [./BC_FreeFlow_EdgePositive]
    type = ArrayOutflowBCCDT
    variable = DD_EdgePositive
    dislocation_character = edge
    dislocation_sign = positive
    boundary = 'left right front back top bottom '
  [../]
  [./BC_FreeFlow_EdgeNegative]
    type = ArrayOutflowBCCDT
    variable = DD_EdgeNegative
    dislocation_character = edge
    dislocation_sign = negative
    boundary = 'left right front back top bottom '
  [../]
  [./BC_FreeFlow_ScrewPositive]
    type = ArrayOutflowBCCDT
    variable = DD_ScrewPositive
    dislocation_character = Screw
    dislocation_sign = positive
    boundary = 'left right front back top bottom '
  [../]
  [./BC_FreeFlow_ScrewNegative]
    type = ArrayOutflowBCCDT
    variable = DD_ScrewNegative
    dislocation_character = Screw
    dislocation_sign = negative
    boundary = 'left right front back top bottom '
  [../]
[]

[Materials]
  [./AdvCoeff]
    type = GenericConstantArray
    prop_name = 'AdvCoeff'
    prop_value = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  new_system = true
  formulation = total
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  [./compute_stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'CP_DiscoFlux'
    rtol = 1.0e-02
  [../]
  [./CP_DiscoFlux]
    type = CrystalPlasticityDislocationEdgeScrew
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.inp
    lattice_friction = 200
    Coeff_hardening = 0.55
    dislo_density_initial = 1.0e+05
    slip_increment_tolerance = 2.0e-2
    DD_EdgePositive = DD_EdgePositive
    DD_EdgeNegative = DD_EdgeNegative
    DD_ScrewPositive = DD_ScrewPositive
    DD_ScrewNegative = DD_ScrewNegative
  [../]
  [./compute_stress_wrapper]
    type = ComputeLagrangianWrappedStress
  [../]
[]

[GlobalParams]
    dislo_density_factor_CDT = 1.0e+05
    C_multi = 8.96e-06
    C_trap = 9.0e-03
    C_m_ann = 0.5
    C_im_ann = 0.5
    burgers_vector_mag = 2.52e-07
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  line_search = 'none'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  nl_abs_tol = 1e-7 
  nl_max_its = 10
  l_abs_tol = 1e-7
  l_max_its = 20

  dtmax = 0.002
  dtmin = 0.00000001
  end_time = 1
  [./TimeStepper]
    type = ConstantDT
    dt = 0.001
    growth_factor = 1.01
  [../]
[]

[AuxVariables]
  [./pk2_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./E_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./shear_stress_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
[./back_stress_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
[./DD_mobile_00]
   order = CONSTANT
   family = MONOMIAL
[../]
[]

[AuxKernels]
  [./pk2_yy] 
   type = RankTwoAux
   variable = pk2_yy
   rank_two_tensor = stress
   index_j = 1
   index_i = 1
  [../]
  [./E_yy]
   type = RankTwoAux
   variable = E_yy
   rank_two_tensor = total_strain
   index_j = 1
   index_i = 1
  [../]
[./shear_stress_00]
   type = MaterialStdVectorAux
   variable = shear_stress_00
   property = applied_shear_stress  
   index = 0
  [../]
[./back_stress_00]
   type = MaterialStdVectorAux
   variable = back_stress_00
   property = back_stress
   index = 0
  [../]
[./DD_mobile_00]
   type = MaterialStdVectorAux
   variable = DD_mobile_00
   property = dislocation_immobile
   index = 0
  [../]
[]

[Postprocessors]
  [./pk2_yy]
    type = ElementAverageValue
    variable = pk2_yy
  [../]
  [./E_yy]
    type = ElementAverageValue
    variable = E_yy
  [../]
[]