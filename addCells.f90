
module addCells
implicit none


contains


!The following grids the cell into smaller cells
!This ASSUMES the cell is NOT already divided up

recursive subroutine divyUp(var1,var2,var3,order,currentsubcell,Nstates)
use f1_parameters
use f1_functions
implicit none
integer, intent(in) :: order, Nstates
integer :: i,j,k,l,index1_1,index1_2,index2_1,index2_2
logical :: state1
integer,intent(in) :: var1,var2,var3
real :: gap1, gap2, gap3
real,allocatable :: coords(:,:),vals(:,:)
integer, dimension(scaling1) :: grid1
integer, dimension(scaling2) :: grid2
integer, allocatable :: indexer(:,:)
character(50),intent(in) :: currentsubcell
character(50) :: subcell1,subcell2,line_data
character(4) :: descriptor1,descriptor2

!This would be what the file is called in the grid or subcell
write(descriptor1,FMT="(I4)") var1
write(descriptor2,FMT="(I4)") var2
subcell1 = trim(currentsubcell)//trim(adjustl(descriptor1))//"_"//&
                trim(adjustl(descriptor2))

open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) "   Subdividing: ", trim(subcell1)
close(80)

!!!!!!!!!!!!!!!!!!!!!!!! NEXT TWO CHECKS MAY NOT BE NECESSARY !!!!!!!!!!!!!!!
!Checks if the file even exists; if not, exit
!inquire(file=trim(path3)//trim(subcell1)//".dat",exist=state1)
!if (.not.(state1)) return
! 
!(meaning it was already organized); if so it then stops
!call system("ls -p "//trim(path3)//trim(currentsubcell)//" | grep '"//trim(subcell1)//"/' > "//&
!!        trim(path4)//trim(temporaryfile2))
!inquire(file=trim(path4)//trim(temporaryfile2),size=i)
!if (i /= 0) return
!inquire(file=trim(path3)//trim(subcell1)//".dat",exist=state1)
!if (.not.(state1)) return
! 
!(meaning it was already organized); if so it then stops
!call system("ls -p "//trim(path3)//trim(currentsubcell)//" | grep '"//trim(subcell1)//"/' > "//&
!        trim(path4)//trim(temporaryfile2))
!inquire(file=trim(path4)//trim(temporaryfile2),size=i)
!if (i /= 0) return
!!!!!!!!!!!!!!!!!!!!!!!! but really are good practice in general !!!!!!!!!!

!Open up the file,read the variables, coordinates, gradients
!vals    ---  holds the to-be-sorted variables (three for now)
!coords  ---  holds the xyz coordinates and gradients (6*Natoms)
!indexer ---  holds the to-be-sorted indexes to access coords
open(72,file=trim(path3)//trim(subcell1)//".dat")
allocate(vals(Nstates,Nvar))
allocate(coords(Nstates,6*Natoms))
allocate(indexer(Nstates,1))
do j=1, Nstates
        read(72,FMT=FMT1,advance="no") (vals(j,i),i=1,Nvar)
        read(72,FMT=FMT2) (coords(j,i),i=1,6*Natoms)
        indexer(j,1) = j
end do
close(72)

!The size of the 'gaps' represents the length of the smaller subcells
!inside the larger (current!) subcell
gap1 = spacing1*(1.0/scaling1)**(order+1)
gap2 = spacing2*(1.0/scaling2)**(order+1)
gap3 = 1.0*(1.0/scaling3)**(order+1)

!Make a new folder with the name var1_var2/ inside the larger path3/currentsubcell/ folder
!This will hold all of the smaller subcells stemming from the current subcell
call system("mkdir -p "//trim(path3)//trim(subcell1))

!Sort the indexed states by the first variable (into columns); then grid it
!This sorts both vals and indexer; indexer can then be used to access coords
call qsort(vals,indexer,Nstates,Nvar,1,Nstates,1)
call grider(grid1,vals,gap1,var1*gap1*scaling1,scaling1,Nstates,Nvar,1,Nstates,1,scaling1,1)
 
index1_1 = grid1(1)
do i = 1, scaling1
        !We always accept these grid elements (cells) in pairs of indexed states [p1,p2)
        !where indexed state p2 is not in the cell i
        !For the last element (cell) the last indexed state must be in it (Nstate)
        !For the first variable, we can think of these are "columns"
        if (i == scaling1) then
                index1_2 = Nstates+1
        else
                index1_2 = grid1(i+1)
        end if

        !If the column is empty, don't even bother
        if (index1_2 == index1_1) then
                index1_1 = index1_2
                cycle
        end if

        !Sort the indexed states in this grid element (into cells); then grid it
        call qsort(vals,indexer,Nstates,Nvar,index1_1,index1_2-1,2)
        call grider(grid2,vals,gap2,var2*gap2*scaling2,scaling2,Nstates,Nvar,&
                        index1_1,index1_2-1,1,scaling2,2)
        index2_1 = grid2(1)
        do j = 1, scaling2
                !For the second variable, we can think of these as "cells"
                if (j == scaling2) then
                        index2_2 = index1_2
                else
                        index2_2 = grid2(j+1)
                end if                
        
                !If there's nothing in the cell, don't even bother
                if (index2_1 == index2_2) then
                        index2_1 = index2_2
                        cycle
                end if

                !Write these states into their respective cell in the grid
                !The array index needs to start at 1 so a subtraction is
                !involved
                write(descriptor1,FMT="(I4)") i-1
                write(descriptor2,FMT="(I4)") j-1
                subcell2 = trim(subcell1)//"/"//trim(adjustl(descriptor1))//&
                                &"_"//trim(adjustl(descriptor2))
                open(72,file=trim(path3)//trim(subcell2)//".dat",status="new")
                do k = index2_1, index2_2-1
                        write(72,FMT=FMT1,advance="no") (vals(k,l),l=1,Nvar)
                        write(72,FMT=FMT2)(coords(indexer(k,1),l),l=1,6*Natoms)
                end do
                close(72)

                !And we also want to keep track of how many states are in this
                !new subcell.
                open(72,file=trim(path3)//trim(subcell2)//".p",status="new")
                write(72,FMT=FMT4) index2_2-index2_1
                close(72)
               
                !The next cell pair [p2,p3) starts at the end of [p1,p2)
                index2_1 = index2_2
        end do

        !The next cell pair [p2,p3) starts at the end of [p1,p2)
        index1_1 = index1_2
end do

end subroutine divyUp





!addState adds a state/frame (its coordinates and gradients stored in coords)
!to all subcells corresponding to var1, var2, var3
!where var1, var2, var3 are stored in vals

subroutine addState(vals,coords)
use f1_parameters
implicit none
integer :: state1,i,j
logical :: flag1
real :: gap1,gap2, var1,var2,var3, var1_old, var2_old, var3_old
integer :: var1_new,var2_new,var3_new,order,line_num
real, dimension(Nvar), intent(in) :: vals
real, dimension(6*Natoms), intent(in) :: coords
character(50) :: line_data, descriptor1, descriptor2, subcell0, subcell1

!var1, var2, var3 stay static to provide a reference (not used)
var1 = vals(1)
var2 = vals(2)
var3 = vals(3)

!var_old changes, essentially only store the trailing digits
var1_old = var1
var2_old = var2

!As we go deeper in the subcells, the length of the subcell gets smaller
!And var_new represents the index of the subcell in the grid
gap1 = spacing1
gap2 = spacing2
var1_new = floor(var1/gap1)
var2_new = floor(var2/gap2)
var3_new = 1

!The order keeps track of how deep we go into subcells
!It is used in the subroutine call for divyUp
order = 0

!The file is named by the indexes of the rounded variables
write(descriptor1,FMT=FMT4) var1_new
write(descriptor2,FMT=FMT4) var2_new
subcell0 = ""
subcell1 = trim(adjustl(descriptor1))//"_"//trim(adjustl(descriptor2))

do
        !Check whether there is some subcell1 (var1,var2) in path3/
        !and reads the var1_var2.p file, increments it
        !to keep tabs on the number of states in var1_var2.dat
        inquire(file=trim(path3)//trim(subcell1)//".dat",exist=flag1)
        if (flag1) then
                open(72,file=trim(path3)//trim(subcell1)//".p",status="old")
                read(72,FMT="(I9)") i
                rewind(72)
                write(72,FMT="(I9)") i + 1
                close(72)
                descriptor1 = "old"
        !If not, we specify two new files: var1_var2.dat, var1_var2.p
        !And we keep tabs on its number of states with var1_var2.p
        else
                open(72,file=trim(path3)//trim(subcell1)//".p",status="new")
                i = 0
                write(72,FMT="(I9)") i + 1
                close(72)
                descriptor1 = "new"
        end if

        !Add this state to its designated subcell in path3/subcell1
        !subcell1 contains the parent subdirectories of subcell1 as well
        open(72,file=trim(path3)//trim(subcell1)//".dat",position="append",status=trim(descriptor1))
        write(72,FMT="(3(1x,F11.6))",advance="no") (vals(j),j=1,Nvar)
        write(72,FMT="(36(1x,F11.6))")(coords(j),j=1,6*Natoms)
        close(72)

        !If there are not many states, this is the 'deepest' subcell; exit
        !If there are many states, then we can go deeper (into a subdirectory)
        !If there is just the perfect number of states, we need to divy it up
        if (i < overcrowd) exit

        if (i == overcrowd) then 
                call divyUp(var1_new,var2_new,var3_new,order,subcell0,i)
                exit
        end if

        !In the other case where there IS a subdirectory then
        !we need to figure out what subcell the state is in

        !We lop off the heading digit used in the previous subcell
        !because we only make use of the trailing one
        !Consequently, we resolve and the subcell gap size shrinks
        var1_old = var1_old - var1_new*gap1
        var2_old = var2_old - var2_new*gap2
        gap1 = gap1/(scaling1)
        gap2 = gap2/(scaling2)
        var1_new = floor(var1_old/gap1)
        var2_new = floor(var2_old/gap2)

        !This new subcell will be named var1_var2
        !And now the path to get there is lengthened
        write(descriptor1,FMT=FMT4) var1_new
        write(descriptor2,FMT=FMT4) var2_new
        subcell0 = subcell1
        subcell1 = trim(subcell1)//"/"//trim(adjustl(descriptor1))//"_"//&
                        trim(adjustl(descriptor2))
 
        !Order represents how "deep" we are in subdirectories
        order = order+1
 
end do



end subroutine addState



end module addCells
