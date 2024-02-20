Dislocation Transport based Crystal Plasticity material model(DiscoFluxM)
=====
This is a Dislocation Transport-based Crystal Plasticity Material Model(DiscoFlux), implemented within MOOSE framework(https://mooseframework.inl.gov/index.html).

# Features
1. The material model represents the dislocation within a crystalline materil in the form of density. The dislocation density at any material point is defined as total_line_length/volume.
2. Total dislocation density is split into four categories, edge_positive, edge_negative, screw_positive, screw_negative.
3. Dislocation densities in each of the slip_systems are epresented as a independent solution variable(coupled). As an example, this leads to 48 solution variables for fcc crystal. 
4. Dislocation densities has both local and nonlocal evolution terms in it. The local evolution is represented in the rate form. The nonlocal evolution is relaized by the transport of dislocation within the grain and across the grain boundary(GB).
5. At the GB, dislocation can transfer one slip_system to a different slip_system acording to the mis-orientation of the two grains. A geometric criterion, incorporating slip_direction, slip_plane_normal and gb_plane_normal in the rotated configuration, is used to determine the direction of transfer.



# Installation instruction


# Instruction to run simulation

# Things to keep in mind
1. This gradient based terms are sensitive to dimension of the domain and mesh size and may cause convergence issue. 
2. It's worthwhile to do some experiment to find proper preconditioner when using Jacobian Free(kind of) solvers (e.g., PJFNK).


# Known Issues

# Featured under development

This is a Dislocation Transport-based Crystal Plasticity Material Model(DiscoFlux), implemented within MOOSE framework(https://mooseframework.inl.gov/index.html).


The material model consists of several solution variables that are several order of magnitude in difference. For an FCC materal there will be  51 solution variables per node. Twelve for each of the Edge_Positive, Edge_Negative, Screw_Positive, Screw_Negative and three displacement degrees of freedom. 
The dislocation density has both the local evolution in the rate-form as well as gradient based evolution(nonlocal). This gradient based terms are sensitive to dimension of the domain and mesh size and may cause convergence issue. It's worthwhile to do some experiment to find proper preconditioner when using PJFNK.



