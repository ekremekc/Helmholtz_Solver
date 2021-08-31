MERGING ALL SEPERATE GEOMETRIES

This way is more convenient to generete geometry for making later changes/adds/removals more elegantly.
Each geometry is generated seperately, then merged with merge.py.
In the merge.py, flame geometry is firstly imported because "0" physical tag for flame domain
cannot be imported.

For reference issue in gmsh forum;

https://gitlab.onelab.info/gmsh/gmsh/-/issues/1535

30 Aug 2021

- The intersection planes are changed for combustion chamber(For rotation from -4 to -3(Last loop has changed to 6), for symmetry from 1 to 4.)
- Conflicting tags problem for different geometries are solved by changing starting tag for each geometries.

31 Aug 2021
- I changed the mesh formats from msh2 to msh4, but in the merge file, I still use msh2 in order to being convertible by meshio.
