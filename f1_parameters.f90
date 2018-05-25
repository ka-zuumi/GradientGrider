module f1_parameters
implicit none


!PATHS
!Path to the trajectories
character(50),parameter :: path1 = "/home/ruisun/proj/ruisun/B0/"
!Path to the f90s and bash scripts
character(50),parameter :: path2 = "/home/kazuumi/lus/B0/branch1/GradientGrider/"
!Path to the grid
character(50),parameter :: path3 = "/home/kazuumi/lus/B0/branch1/grid/"
!Path for the following files (see below)
character(50),parameter :: path4 = "/home/kazuumi/lus/B0/branch1/"



!FILES
!File that will keep the trajectory folder names
character(20),parameter :: trajectories = "f1_trajectories.txt" !Used to be file1
!File that writes the progress of the program
character(20),parameter :: progressfile = "f1_progress.txt" !Used to not exist
!File to write to for system calls
character(20),parameter :: temporaryfile1 = "tmp1.txt"
character(20),parameter :: temporaryfile2 = "tmp2.txt"
character(20),parameter :: temporaryfile3 = "tmp3.txt"




!If we need to create a grid from the trajectories in path1, then set this true
logical :: start_from_scratch = .true.



!FORMATS
character(13),parameter :: FMT1 = "(3(1x,F11.6))"           !For three variables
character(14),parameter :: FMT2 = "(36(1x,F11.6))"          !For six atoms, full state
character(14),parameter :: FMT3 = "(18(1x,F11.6))"          !For six atoms, coords
character(4),parameter :: FMT4 = "(I4)"                     !For subcell names



!VARIABLES

!The number of atoms in the system
integer,parameter :: Natoms = 6

!The number of variables in use
integer,parameter :: Nvar = 3

!The spacing is the dimension of the gridlines; currently var3 is not grided
real,parameter :: spacing1 = 1.0
real,parameter :: spacing2 = 1.0
real,parameter :: spacing3 = 0.1

!The threshold of "overcrowded" grid
integer,parameter :: overcrowd = 1000

!The scaling is the amount that is resolved for an overcrowded grid
! (10,10,10) = x1000 magnification
! This might be a little too many
integer,parameter :: scaling1 = 10
integer,parameter :: scaling2 = 10
integer,parameter :: scaling3 = 10

!There are some outliers; making a maximum throws these away
real,parameter :: max_var1 = 80.0
real,parameter :: max_var2 = 80.0

end module f1_parameters



