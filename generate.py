#!/usr/bin/python
"""
Script to generate a range of XDMF files with Finite Element "Index + Data" format
"""

from dolfin import *

encoding = XDMFFile.Encoding_HDF5
savedir = "generated_files"

def output_scalar_xdmf(mesh, scalar_expression, file_suffix, elem_family, elem_degree):
    Q = FunctionSpace(mesh, elem_family, elem_degree)
    F = Function(Q)
    F.interpolate(Expression(scalar_expression, degree=1))
    xdmf = XDMFFile(savedir + "/scalar_" + elem_family + str(elem_degree) + "_" + file_suffix + ".xdmf")
    xdmf.write_checkpoint(F, "f", 0, encoding)


def output_vector_xdmf(mesh, vector_expression, file_suffix, elem_family, elem_degree):
    Q = VectorFunctionSpace(mesh, elem_family, elem_degree)
    F = Function(Q)
    F.interpolate(Expression(vector_expression, degree=1))
    xdmf = XDMFFile(savedir + "/vector_" + elem_family + str(elem_degree) + "_" + file_suffix + ".xdmf")
    xdmf.write_checkpoint(F, "f", 0, encoding)

def output_tensor_xdmf(mesh, tensor_expression, file_suffix, elem_family, elem_degree):
    Q = TensorFunctionSpace(mesh, elem_family, elem_degree)
    F = Function(Q)
    F.interpolate(Expression(tensor_expression, degree=1))
    xdmf = XDMFFile(savedir + "/tensor_" + elem_family + str(elem_degree) + "_" + file_suffix + ".xdmf")
    xdmf.write_checkpoint(F, "f", 0, encoding)

# 2D
mesh = UnitSquareMesh(10, 10)
for elem_family in ["CG", "DG"]:
    elem_degrees = [1, 2, 4]
    if elem_family == "DG":
        elem_degrees.append(0)
    for elem_degree in elem_degrees:
        output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "2D", elem_family, elem_degree)
        output_vector_xdmf(mesh, ("sin(2*pi*x[0])*cos(2*pi*x[1])", "sin(2*pi*x[0])"), "2D", elem_family, elem_degree)
        output_tensor_xdmf(mesh, (("sin(2*pi*x[0])*cos(2*pi*x[1])", "sin(2*pi*x[0])"), ("1.0", "2.0")), "2D", elem_family, elem_degree)


for elem_family in ["RT"]:
    for elem_degree in [1]:
        output_scalar_xdmf(mesh, ("sin(2*pi*x[0])", "cos(2*pi*x[1])"), "2D", elem_family, elem_degree)

# 2D manifold
mesh = UnitDiscMesh.create(mpi_comm_world(), 5, 1, 3)
mesh.init_cell_orientations(Constant((0, 0, 1)))
for elem_family in ["CG", "DG"]:
    for elem_degree in [1, 2]:
        output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "2D3", elem_family, elem_degree)
        output_vector_xdmf(mesh, ("sin(2*pi*x[0])", "cos(2*pi*x[1])", "1.0"), "2D3", elem_family, elem_degree)

# 2D quadratic triangles
mesh = UnitDiscMesh.create(mpi_comm_world(), 2, 2, 2)
for elem_family in ["CG"]:
    for elem_degree in [1, 2]:
        output_scalar_xdmf(mesh, "x[0]*x[1]", "2D_quad_triangle", elem_family, elem_degree)

# 3D
mesh = UnitCubeMesh(5, 5, 5)
for elem_family in ["CG", "DG"]:
    for elem_degree in [1]:
        output_scalar_xdmf(mesh, "sin(2*pi*x[0])*cos(2*pi*x[1])", "3D", elem_family, elem_degree)
        output_vector_xdmf(mesh, ("sin(2*pi*x[0])", "cos(2*pi*x[1])", "1.0"), "3D", elem_family, elem_degree)
        output_tensor_xdmf(mesh, (("sin(2*pi*x[0])", "cos(2*pi*x[1])", "3.0"), ("4.0", "5.0", "6.0"), ("7.0", "8.0", "9.0")), "3D", elem_family, elem_degree)
for elem_family in ["RT"]:
    for elem_degree in [1]:
        output_scalar_xdmf(mesh, ("sin(2*pi*x[0])", "cos(2*pi*x[1])", "1.0"), "3D", elem_family, elem_degree)
