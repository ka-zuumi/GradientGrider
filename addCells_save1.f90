
module addCells
implicit none


contains



!divyUp takes a certain var1, var2, var3 (ex. 1, 5, 22)
recursive subroutine divyUp(var1,var2,var3,order,currentsubcell)
use f1_parameters
use f1_functions
implicit none
integer, intent(in) :: order
integer :: Nstates,Nlength,i,j,k,l,index1_1,index1_2,index2_1,index2_2,order_new
logical :: state1
integer,intent(in) :: var1,var2,var3
real :: gap1, gap2, gap3
real,allocatable :: coords(:,:),vals(:,:),temp_coords(:,:),temp_vals(:,:)
integer, dimension(scaling1) :: grid1
integer, dimension(scaling2) :: grid2
integer, allocatable :: indexer(:,:),temp_indexer(:,:)
character(50),intent(in) :: currentsubcell
character(50) :: subcell1,subcell2,line_data
character(4) :: descriptor1,descriptor2

Nlength = overcrowd*2

!This would be what the file is called in the grid or subcell
write(descriptor1,FMT="(I4)") var1
write(descriptor2,FMT="(I4)") var2
subcell1 = trim(currentsubcell)//trim(adjustl(descriptor1))//"_"//trim(adjustl(descriptor2))

!Checks if the file even exists
inquire(file=trim(path3)//trim(subcell1)//".dat",exist=state1)
if (.not.(state1)) return

!Checks if the file already has a folder named after it
!(meaning it was already organized); if so it then stops
call system("ls -p "//trim(path3)//trim(currentsubcell)//" | grep '"//subcell1//"/' > "//&
        trim(path2)//"/"//trim(temporaryfile))
inquire(file=trim(path2)//"/"//trim(temporaryfile),size=i)
if (i /= 0) return

!Open up the file,read the variables, coordinates, gradients
open(72,file=trim(path3)//trim(subcell1)//".dat")
allocate(vals(Nlength,Nvar),coords(Nlength,6*Natoms),indexer(Nlength,1))
order_new = order + 1
Nstates = 0
do
        Nstates = Nstates + 1
                !If there is an overflow of states, increase the array length
                !warning: this is pretty spaghetti
                if (Nstates > Nlength) then
                Nlength = Nlength*2
                allocate(temp_vals(Nlength,Nvar),temp_coords(Nlength,6*Natoms),temp_indexer(Nlength,1))
                temp_vals(1:Nstates,:) = vals
                temp_coords(1:Nstates,:) = coords
                temp_indexer(1:Nstates,:) = indexer
                deallocate(vals,coords,indexer)
                allocate(vals(Nlength,Nvar),coords(Nlength,6*Natoms),indexer(Nlength,1))
                vals = temp_vals
                coords = temp_coords
                indexer = temp_indexer
                deallocate(temp_vals,temp_coords,temp_indexer)
                end if
        read(72,FMT=FMT1,advance="no",iostat=j) (vals(Nstates,i),i=1,Nvar)
        if (j /= 0) then
                Nstates = Nstates - 1
                exit
        end if
        read(72,FMT=FMT2) (coords(Nstates,i),i=1,6*Natoms)
        indexer(Nstates,1) = Nstates
end do

if (Nstates < overcrowd) return

!The size of the 'gaps' is one over the number of gaps, raised to the order
!where the order is how far we have magnified
gap1 = spacing1*(1.0/scaling1)**order_new
gap2 = spacing2*(1.0/scaling2)**order_new
gap3 = 1.0*(1.0/scaling3)**order_new

!The following grids the cell into smaller cells
!This ASSUMES the cell is NOT already divided up

!Make a new folder with the name var1_var2/ inside the larger path3/currentsubcell/ folder
call system("mkdir -p "//trim(path3)//trim(subcell1))

!Sort the indexed states by the first variable (into columns); then grid it
call qsort(vals,indexer,Nlength,Nvar,1,Nstates,1)
call grider(grid1,vals,gap1,var1*gap1*scaling1,scaling1,Nlength,Nvar,1,Nstates,1,scaling1,1)
!print *, "The first-level grid looks like:", grid1
!print *, "With cell spacings of:", gap1
!print *, "With values starting at:", var1*gap1*scaling1
index1_1 = grid1(1)
print *, ""
do i = 1, scaling1
        !We always accept these grid elements (cells) in pairs of indexed states [p1,p2)
        !where indexed state p2 is not in the cell i
        !For the last element (cell) the last indexed state must  be in it (Nstate)
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
        call qsort(vals,indexer,Nlength,Nvar,index1_1,index1_2-1,2)
        call grider(grid2,vals,gap2,var2*gap2*scaling2,scaling2,Nlength,Nvar,&
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
               
                !If this cell is overcrowded, just call divyUp again 
                if (index2_2-index2_1 > overcrowd) then
                        call divyUp(i-1,j-1,var3,order_new,trim(subcell1))
                end if

                index2_1 = index2_2
        end do

        index1_1 = index1_2
end do

end subroutine divyUp





!addState adds a state (its coordinates and gradients stored in coords)
!to all subcells corresponding to var1, var2, var3
!Thus, var1, var2, var3 need to be calculated beforehand

!WARNING, not yet updated to reflect new filenames and paths

subroutine addState(var1,var2,var3,vals,coords)
use f1_parameters
implicit none
integer :: state1,i
logical :: flag1
real, intent(in) :: var1, var2, var3
real :: gap1,gap2
integer :: var1_new,var2_new,var3_new,var1_old,var2_old,var3_old,order,line_num
real, dimension(Nvar), intent(in) :: vals
real, dimension(6*Natoms), intent(in) :: coords
character(50) :: line_data, descriptor1, descriptor2, subcell1

var1_old = var1
var2_old = var2
var3_old = var3

gap1 = spacing1
gap2 = spacing2
var1_new = floor(var1/gap1)
var2_new = floor(var2/gap2)
var3_new = 1
order = 0

write(descriptor1,FMT="(A50)") var1_new
write(descriptor2,FMT="(A50)") var2_new
subcell1 = trim(descriptor1)//"_"//trim(descriptor2)

do
        !Check whether var1_var2 is in path3/subcell1
        !If not, we specify this is a new file
        inquire(file=trim(path3)//subcell1//".dat",exist=flag1)
        if (flag1) then
                descriptor1 = "old"
        else
                descriptor1 = "new"
        end if

        !Add this state to var1_var2 in path3/subcell1
        open(72,file=trim(path3)//trim(subcell1)//".dat",position="append",status=descriptor1)
        write(72,FMT="(3(1x,F11.6))",advance="no") (vals(i),i=1,Nvar)
        write(72,FMT="(36(1x,F11.6))")(coords(i),i=1,6*Natoms)
        close(72)

        !We check for a subcell in a subdirectory var1_var2
        call system("find "//trim(path3)//trim(subcell1)//"/ > "//trim(path4)//&
                    trim(temporaryfile))
        inquire(file=trim(path4)//trim(temporaryfile),size=i)
        
        !START OF IF STATEMENT
        !If there isn't a subdirectory for var1_var2, we need to check if
        !there is overcrowding and make one
        if (i == 0) then
                line_num = 0
                open(72,file=trim(path3)//trim(subcell1)//".dat")
                do
                        line_num = line_num + 1
                        read(72,FMT="(A50)",iostat=state1) line_data
                        if (state1 /= 0) exit
                end do
                close(72)

        !Just call divyUp to sort that one out
        if (line_num > overcrowd) then
                order = order + 1
                call divyUp(var1_new,var2_new,var3_new,order,trim(subcell1))
        end if

        !This would be the end of the recursive subroutine
        !As we go deeper and deeper, subcell1 grows longer and longer
        !Because of ordinality, eventually this ends!
        exit
        end if
        !END OF IF STATEMENT

        !In the other case where there IS a subdirectory then
        !we need to figure out what subcell the state is in
        var1_old = var1 - var1_new*gap1
        var2_old = var2 - var2_new*gap2
        gap1 = gap1/(scaling1)
        gap2 = gap2/(scaling2)
        var1_new = floor(var1_old/gap1)
        var2_new = floor(var2_old/gap2)

        !It will be named var1_var2
        write(descriptor1,FMT="(A50)") var1_new
        write(descriptor2,FMT="(A50)") var2_new
        subcell1 = subcell1//"/"//trim(descriptor1)//"_"//trim(descriptor2)
 
        !Order represents how "deep" we are in subdirectories
        order = order+1
 
end do



end subroutine addState



end module addCells
