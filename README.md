Dislocation Transport based Crystal Plasticity material model(DiscoFluxM)
=====
This is a Dislocation Transport-based Crystal Plasticity Material Model(DiscoFlux), implemented within MOOSE framework(https://mooseframework.inl.gov/index.html).

# Features
1. The material model represents the dislocation within a crystalline material in the form of density. The dislocation density at any material point is defined as total_line_length/volume.
2. Total dislocation density is split into four categories: edge_positive, edge_negative, screw_positive, screw_negative.
3. Dislocation densities in each of the slip_systems are represented as a independent solution variable(coupled). As an example, this leads to 48 solution variables for fcc crystal plus three displacements. Total 51 DOFs per node.
4. Dislocation densities have both local and nonlocal evolution terms in it. The local evolution is represented in the rate form. The nonlocal evolution is relaized by the transport of dislocation within the grain and transfer of dislocation across the grain boundary(GB).
5. At the GB, dislocations can get transferred from one slip_system to a different one acording to the mis-orientation of the two grains. A geometric criterion, incorporating slip_direction, slip_plane_normal and gb_plane_normal in the rotated configuration, is used to determine the direction of transfer.



# Installation instruction
1. Download/clone the modified version of moose from LANL-Gitlab.
2. Set-up the conda and required libraries according to the instruction given in [moose website](https://mooseframework.inl.gov/getting_started/installation/conda.html)

# Instruction to run simulation

# Things to keep in mind
1. The gradient based terms in the computation of brack_stress are sensitive to physical dimension of the domain and mesh size. This may cause convergence issue, if numerical parameters are not set properly. 
2. It's worthwhile to do some experiment to find proper preconditioner when using Jacobian Free(kind of) solvers (e.g., PJFNK).


# Known Issues

# Featured under development




