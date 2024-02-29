[Mesh]
 [./MeshFile]
  type = FileMeshGenerator 
  file = './Mesh/block_n10_01_SM.inp'
 [../]
  [surface_GB]
    type = GrainBoundaryIndentifier
    input = MeshFile
    new_boundary = 'surface_GB'
    #all_blocks = '0 1' #'2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17'
  [../]

[]

[Outputs]
    csv = true
    exodus = true
[]

