MERGING ALL SEPERATE GEOMETRIES

This way is more convenient to generete geometry for making later changes/adds/removals more elegantly.
Each geometry is generated seperately, then merged with merge.py.
In the merge.py, flame geometry is firstly imported because "0" physical tag for flame domain
cannot be imported.


30 Aug 2021

- The intersection planes are changed for combustion chamber(For rotation from -4 to -3(Last loop has changed to 6), for symmetry from 1 to 4.)
- Conflicting tags problem for different geometries are solved by changing starting tag for each geometries.
