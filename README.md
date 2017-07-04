
# xdmf-fe-filter

ParaView filter to read Finite Element data in XDMF format

## Format
At present, the XDMF contains an Unstructured Uniform Grids, i.e. mesh topology and geometry, 
and four *Attribute* elements. Their names are in format *{functionname}_{element}_{desc}* where
*functonname* is name of the finite element function saved, *element* is an abbreviation for
finite element family and degree, i.e. one of:

 * CG1, CG2,
 * DG0, DG1, DG2,
 * RT1,

and *desc* is a description of data stored: could be:

 * cell_dofs - global dofmap (degrees of freedom map),
 * cell_dofs_x - number of degrees of freedom in each cell,
 * vector - values of degrees of freedom,
 * cells - indices of cells, i.e. ordering of cells.

In order to read correctly, use the Xdmf3 reader in ParaView.

A script `generate.py` is provided, which uses experimental XDMF output in FEniCS to generate the XDMF files.

## Testing in ParaView

 1. Open the XDMF file with the Xdmf3 reader
 2. Add a filter "ProgrammableFilter" and in the box paste: `execfile("fe-filter.py")` - you may need to add the full path.
 3. Select the output field of the filter to display and process as usual

## Future Work
