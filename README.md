
To make the grid, enter in "make -f make_getCells"
This has one main program: getCells.f90
And four modules: f1_parameters.f90, f1_functions.f90, f1_variables.f90, and addCells.f90

The main program goes roughly like follows:
Go into the directory containing the trajectories and pick a trajectory, one by one.
Use text manipulation to make it more fortran-friendly.
Read the formatted data frames; these are six-liners that describe an xyz coordinate and gradient.
Calculate the variables desired from these (right now, r1 and r2).
Bin a particular frame based on these variables and the desired spacing of the grid.
Add the state thoughout every cell and subcell of the grid it belongs to with "addState".
Repeat this for every frame of the trajectory, and every trajectory in the larger directory

"addState" is a subroutine of addCells.f90
Taking as input the values of the variables and the xyz,gradient coordinates, bin the frame.
To 'bin' a frame, it simply appends onto a file the frame data.
The file is named to represent the bin; ex. r1,r2 -> r1_r2.dat.
To keep track of the number of frames in a bin, a second file r1_r2.p is made; it has this number.
When the bin has a specific number of states (it is 'overcrowded'), "divyUp" is called.
If it has more states than this number, then it must also populate one of this bin's subcells.
These are located in a folder of the same name r1_r2; each smaller subcell in this folder is indexed
by their variable values position within the original subcell r1_r2.
In this way, a frame is added and to its deepest subcell.

"divyUp" is a subroutine of addCells.f90
Taking as input the values of the variables, the subcell name, path, and depth, it creates a
subdirectory that has smaller subcells of the original subcell and the frames belong to them.
It does this by first collecting all of the variable values and coordinates into an array.
Then it sorts the values, then it calculates which indexes belong to which bin/cell.
For each non-empty cell, it creates a new file r1_r2.dat where r1 and r2 now represent
the index of the new, smaller subcell within the original subcell.


