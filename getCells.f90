program getCells
use f1_parameters
use f1_variables
use f1_functions
use addCells
implicit none
integer :: Ntraj,Nstates,state1,state2,state3,line_num,skips,i,j,k
logical :: flag1, flag2
real :: var1,var2,var3
integer :: var1_int, var2_int, var3_int
real :: t1,t2
real,allocatable :: coords(:),vals(:)
character(50) :: current_path,line_data
character(20) :: subcell,some_folder,descriptor1,descriptor2


!Instead of printing to the terminal, we print to the progressfile
!First we check if it already exists and just wipe it
!We want to start from a clean slate
inquire(file=trim(path4)//trim(progressfile),exist=flag1)
if (flag1) call system("rm "//trim(path4)//trim(progressfile))
open(70,file=trim(path4)//trim(progressfile),status="new")
write(70,FMT="(A50)") "Let's start"
write(70,FMT="(A50)") ""
close(70)

!All trajectory folders are formatted as a number
!So search for these numbered folders and read them
call system("ls -p "//trim(path1)//" | grep '[0123456789]/' > "//trim(path4)//trajectories)
open(70,file=trim(path4)//trim(trajectories),action="read")

!coords  ---  this array keeps the coordinates and gradients of a state
!vals    ---  this array keeps var1, var2, var3 of a state
allocate(coords(6*Natoms))
allocate(vals(Nvar))
skips = Natoms+1

!We will keep track of how many trajectories and states we encounter
call CPU_time(t1)
Ntraj = 0
Nstates = 0
do

        if (.not. (start_from_scratch)) exit

        !Fetch the name of one folder (a trajectory)
        !Format its contents with sed into tmp.txt
        read(70,FMT="(A20)",iostat=state1) some_folder
        if (state1 /= 0) exit
        call system(trim(path2)//"trajectory_sed.sh "//trim(path1)//&
                trim(some_folder)//"kkk.out "//trim(path4)//trim(temporaryfile))
        Ntraj = Ntraj+1



!To track the progress, open up 80
open(80,file=trim(path4)//trim(progressfile),position="append")
write(80,*) trim(path2)//"trajectory_sed.sh "//trim(path1)//&
            trim(some_folder)//"kkk.out "//trim(path4)//trim(temporaryfile)
write(80,*) "Now accessing folder:", some_folder
close(80)


        
        !Open the now-formatted trajectory, discard the first line
        open(71,file=trim(path4)//trim(temporaryfile))
        read(71,FMT="(A50)", iostat=state2) line_data
        line_num = 1

        do
                if (line_num == skips) then

                        read(71,FMT="(A50)", iostat=state2) line_data
                        if (state2 /= 0) exit
                        line_num = 1

                        !With the fully described state, calculate the
                        !variables wanted
                        call getVar1(coords(1:3*Natoms),Natoms,var1)
                        call getVar2(coords(1:3*Natoms),Natoms,var2)
                        call getVar3(coords(1:3*Natoms),Natoms,var3)

                        !If they are outliers, just skip this cycle
                        if ((var1 > max_var1).or.(var2 > max_var2)) then
                                Nstates = Nstates - 1
                                cycle
                        end if

                        vals(1) = var1
                        vals(2) = var2
                        vals(3) = var3

                        !Find which appropriate cell the state is in
                        var1_int = floor(var1/spacing1)
                        var2_int = floor(var2/spacing2)
                        var3_int = floor(var3/spacing3)

                        !Make the filename for the cell
                        write(descriptor1,FMT=FMT4) var1_int
                        write(descriptor2,FMT=FMT4) var2_int
                        subcell = trim(adjustl(descriptor1))//"_"//trim(adjustl(descriptor2))//".dat"
                        inquire(file=trim(path3)//trim(subcell),exist=flag1)

                        !Write to the file the variables, coordinates, and gradients
                        if (flag1) then
                                open(72,file=trim(path3)//trim(subcell),position="append")
                                write(72,FMT=FMT1,advance="no") (vals(i),i=1,Nvar)
                                write(72,FMT=FMT2)(coords(i),i=1,6*Natoms)
                                close(72)
                        else
                                open(72,file=trim(path3)//trim(subcell),position="append",status="new")
                                write(72,FMT=FMT1,advance="no") (vals(i),i=1,Nvar)
                                write(72,FMT=FMT2)(coords(i),i=1,6*Natoms)
                                close(72)
                        end if

                else

                        !If the state is not yet fully described, continue
                        !adding coordinates; (line_num keeps track of atom #)
                        call read_coords(coords,Natoms,line_num)
                        line_num = line_num + 1
                end if
        end do              
        close(71)
!       if (Ntraj == 20) exit        (if we wanted to end data collection prematurely)
end do
call CPU_time(t2)
close(70)

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*), ""
write(70,*) "The reading and collecting took:", t2-t1, "seconds"
write(70,*) "The number of states is:", Nstates
close(70)
t1 = t2

deallocate(vals,coords)


















!Now we just need to organize every folder
!Note: probably more efficient do this inside the loop
!maybe by using the addState subroutine instead
line_data = ""
do i = 1, floor(max_var1/spacing1)-1
do j = 1, floor(max_var2/spacing1)-1

!open(70,file=trim(path4)//trim(progressfile),position="append")
!write(70,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
!write(70,*) ""
!write(70,*) "divying up   var1:", i, "var2:", j
!write(70,*) ""
!close(70)

!Access divyUp from the addCells module
call divyUp(i,j,1,0,line_data)

end do
end do




call CPU_time(t2)

open(70,file=trim(path4)//trim(progressfile),position="append")
write(70,*), ""
write(70,*) "The divyUp took:", t2-t1, "seconds"
close(70)






end program getCells






!maybe keep all distance subroutine into one mod (e.g. f1_variables.f90)
subroutine read_coords(coords,Natoms,line_num)
implicit none
integer, intent(in) :: Natoms, line_num
!integer, intent(out) :: stat
integer :: i,j
character(11) :: cvar
real, dimension(6*Natoms), intent(out) :: coords

! again, call the dis funct. here

!come back to later !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!this only has the coordinates of one atom (and its gradients)
!when line_num ==6 I can call the distance function, admittedly


i = 3*line_num
j = 3*Natoms
read(71,FMT="(10x,3(1x,F10.6),1x,3(1x,F10.6))") coords(i-2), coords(i-1), &
coords(i), coords(j+i-2), coords(j+i-1), coords(j+i)


end subroutine
