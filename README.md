# qp1qc

An R package for solving (possibly non-convex) quadratic programming problems with 1 quadratic constraint. That is, this package solves problems of the form

<a href="https://www.codecogs.com/eqnedit.php?latex=\min_{x\in&space;R^p}\,\,&&space;x^\top&space;Ax&space;&plus;&space;x^\top&space;a&space;\\&space;\text{such&space;that}\,\,\,&space;&&space;x^\top&space;Bx&space;&plus;&space;x^\top&space;b&space;&plus;k&space;\leq&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\min_{x\in&space;R^p}\,\,&&space;x^\top&space;Ax&space;&plus;&space;x^\top&space;a&space;\\&space;\text{such&space;that}\,\,\,&space;&&space;x^\top&space;Bx&space;&plus;&space;x^\top&space;b&space;&plus;k&space;\leq&space;0" title="\min_{x\in R^p}\,\,& x^\top Ax + x^\top a \\ \text{such that}\,\,\, & x^\top Bx + x^\top b +k \leq 0" /></a>

For this implementation, one of the matrices A or B must be positive definite, but the other need not be positive definite or semi-positive definite.
