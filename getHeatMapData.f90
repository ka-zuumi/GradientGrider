program getHeatMapData
use f1_parameters
implicit none
logical :: flag1
integer :: i, j, cells_var1, cells_var2,states,state1
integer,allocatable :: cells(:,:)
character(20) :: descriptor1, descriptor2, subcell, line_data


!How many cells are in the first-level grid
cells_var1 = floor(max_var1/spacing1)
cells_var2 = floor(max_var2/spacing2)
allocate(cells(cells_var1,cells_var2))

!Iterate through these cells
do i = 1, cells_var1
do j = 1, cells_var2

!Get the file corresponding to this cell
write(descriptor1,FMT=FMT4) i
write(descriptor2,FMT=FMT4) j
subcell = trim(adjustl(descriptor1))//"_"//trim(adjustl(descriptor2))

!Figure out whether this file exists
inquire(file=trim(path3)//trim(subcell)//".dat",exist=flag1)

!If it does, figure out how many states are in it
if (flag1) then
open(72,file=trim(path3)//trim(subcell)//".p")
read(72,FMT="(I9)") states
close(72)

!Otherwise, this cell is empty
else
states = 0
end if

cells(i,j) = states

end do
end do



!Then print the data onto temporaryfile3
open(72,file=trim(path4)//trim(temporaryfile3))
write(descriptor1,FMT=FMT4) cells_var1
do j = 1, cells_var2
write(72,FMT="("//trim(adjustl(descriptor1))//"I9)") cells(:,j)
end do




!Then call gnuplot to make the image
call system("gnuplot < makeHeatMap")

end program getHeatMapData
