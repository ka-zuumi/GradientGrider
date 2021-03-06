program getHeatMapData
use f1_parameters
implicit none
logical :: flag1
integer :: i, j, cells_var1, cells_var2,states,state1
integer,allocatable :: cells(:,:)
character(20) :: descriptor1, descriptor2, file2, file3, line_data

!Name of the file
file3 = "population.dat"

!How many cells are in the first-level grid
cells_var1 = floor(max_var1/spacing1)
cells_var2 = floor(max_var2/spacing2)
allocate(cells(cells_var1,cells_var2))

!Iterate through these cells
do i = 1, cells_var1
do j = 1, cells_var2
states = 0

!Get the file corresponding to this cell
write(descriptor1,"(I4)") i
write(descriptor2,"(I4)") j
file2 = trim(adjustl(descriptor1))//"_"//trim(adjustl(descriptor2))

!Figure out whether this file exists
inquire(file=trim(path3)//trim(file2)//".dat",exist=flag1)

!If it does, figure out how many states are in it
if (flag1) then
open(72,file=trim(path3)//trim(file2)//".dat")
do
read(72,FMT="(A20)",iostat=state1) line_data
if (state1 /= 0) exit
states = states + 1
end do
close(72)

!Otherwise, this cell is empty
else
end if

cells(i,j) = states
print *, states

end do
end do



!Then print the data onto file3
open(72,file=trim(path2)//trim(file3))
write(descriptor1,"(I9)") cells_var1
do j = 1, cells_var2
write(72,FMT="("//trim(adjustl(descriptor1))//"I9)") cells(:,j)
end do


end program getHeatMapData
