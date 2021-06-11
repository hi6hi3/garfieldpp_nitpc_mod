gmsh gemcell.geo -3 -optimize -order 2
ElmerGrid 14 2 gemcell.msh -autoclean
ElmerSolver gemcell.sif
