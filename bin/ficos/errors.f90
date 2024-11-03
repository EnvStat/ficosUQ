MODULE wq_error
!! Error codes

integer :: error_num_input_args=8       ! Incorrect number of input arguments
integer :: error_num_input_args_ncin=18 ! Incorrect number of input arguments in ncin
integer :: error_num_no_stratigraphy=19 ! Stratigraphy flag not defined
integer :: error_num_ini_advection=31   ! Advection set wrongly in ini-file
integer :: error_num_ini_integration=32 ! Integration set wrongly in ini-file

!! Subroutines
CONTAINS

SUBROUTINE wqkiirikki_error(error_number,error_number2)
!Error handler for fatal errors
!Exits the program with the error number
IMPLICIT NONE

integer :: error_number
integer,optional,intent(in) :: error_number2

select case(error_number)
case (8)
   write(0,*)"FATAL ERROR: incorrect number of input arguments in wqficos"
   write(0,*)"wqficos usage: wqficos <parameterfile> <resultdir> <hdfile> <loadfile> &
              &<ini-file> <modelname> <start timestep> <end timestep>"
   write(0,*)"Result directory will be created if does not exist."
   call exit(error_number)
case (9)
   write(0,*)"FATAL ERROR: error in ODE-solver"
   write(0,*)"ODE solver error number:",error_number2
   call exit(error_number)
case(18)
   write(0,*)"FATAL ERROR: incorrect number of input arguments in ncin"
   write(0,*)"ncin usage: ncin <start date> <end date> <ncddir> &
              &<boundaryconditions> <solarradiation> <modelpointsfile> [<ini-file>]"
   write(0,*)"Start date and end date format YYYYMM"
   write(0,*)"ncdir is the directory where the Netcdf-files are."
   call exit(error_number)
case(19)
   write(0,*)"FATAL ERROR: no stratigraphy flag defined"
   write(0,*)"It should be one of: 1: thermocline, 2: halocline or 3: pycnocline"
   call exit(error_number)
case(31)
   write(0,*)"INI-FILE ERROR: advection error:",error_number2," Advection set to 1"
   return
case(32)
   write(0,*)"INI-FILE ERROR: integration error:",error_number2," Integration set to 1"
   return
case default
   write(0,*)"FATAL ERROR: unexpected error"
   call exit(error_number)
end select

END SUBROUTINE

END MODULE wq_error
