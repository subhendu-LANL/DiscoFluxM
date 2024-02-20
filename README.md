Dislocation Transport based Crystal Plasticity material model(DiscoFluxM)
=====

# Features


# Installation instruction


# Instruction to run simulation

# Necessary parts of the model

# Known Issues

# Featured under development

This is a Dislocation Transport-based Crystal Plasticity Material Model(DiscoFlux), implemented within MOOSE framework(https://mooseframework.inl.gov/index.html).


The material model consists of several solution variables that are several order of magnitude in difference. For an FCC materal there will be  51 solution variables per node. Twelve for each of the Edge_Positive, Edge_Negative, Screw_Positive, Screw_Negative and three displacement degrees of freedom. 
The dislocation density has both the local evolution in the rate-form as well as gradient based evolution(nonlocal). This gradient based terms are sensitive to dimension of the domain and mesh size and may cause convergence issue. It's worthwhile to do some experiment to find proper preconditioner when using PJFNK.



