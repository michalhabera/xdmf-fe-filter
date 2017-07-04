#!/usr/bin/python
"""
Script to generate a range of XDMF files with Finite Element "Index + Data" format
"""

from dolfin import *


def output_scalar_xdmf(mesh, scalar_expression, file_suffix, elem_family, elem_degree):
    Q = FunctionSpace(mesh, elem_family, elem_degree)
    F = Function(Q)
    F.interpolate(Expression(scalar_expression, degree=1))
    xdmf = XDMFFile("scalar_" + elem_family + str(elem_degree) + "_" + file_suffix + ".xdmf")
    xdmf.write_checkpoint(F, "f", 0, XDMFFile.Encoding_ASCII)


def output_vector_xdmf(mesh, vector_expression, file_suffix, elem_family, elem_degree):
    Q = VectorFunctionSpace(mesh, elem_family, elem_degree)
    F = Function(Q)
    F.interpolate(Expression(vector_expression, degree=1))
    xdmf = XDMFFile("vector_" + elem_family + str(elem_degree) + "_" + file_suffix + ".xdmf")
    xdmf.write_checkpoint(F, "f", 0, XDMFFile.Encoding_ASCII)


# 2D
mesh = UnitSquareMesh(10, 10)
for elem_family in ["CG", "DG"]:
    for elem_degree in [1, 2]:
        output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "2D", elem_family, elem_degree)

# 2D manifold
mesh = UnitDiscMesh(mpi_comm_world(), 5, 1, 3)
mesh.init_cell_orientations(Constant((0, 0, 1)))
for elem_family in ["CG", "DG"]:
    for elem_degree in [1, 2]:
        output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "2D3", elem_family, elem_degree)
        output_vector_xdmf(mesh, ("sin(2*pi*x[0])", "cos(2*pi*x[1])", "1.0"), "2D3", elem_family, elem_degree)

# 3D
mesh = UnitCubeMesh(5, 5, 5)
for elem_family in ["CG", "DG"]:
    for elem_degree in [1]:
        output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "3D", elem_family, elem_degree)
        output_vector_xdmf(mesh, ("sin(2*pi*x[0])", "cos(2*pi*x[1])", "1.0"), "3D", elem_family, elem_degree)
for elem_family in ["RT"]:
    for elem_degree in [1]:
        output_scalar_xdmf(mesh, ("sin(2*pi*x[0])", "cos(2*pi*x[1])", "1.0"), "3D", elem_family, elem_degree)