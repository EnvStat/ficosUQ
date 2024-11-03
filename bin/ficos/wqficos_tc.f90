!==================================================================
!
! Water quality module for the Finnish coastal modelling system
!
!------------------------------------------------------------------
!
! License: (TBD) (CC BY 4.0?, MIT?)
! 
! Finnish Environment Institute
!
! Code authors:
! * Kai Rasmus, Luode Consulting Oy, 2015-2016
! * Janne Ropponen, Finnish Environment Institute, 2017-2019
!==================================================================

! Notes regarding input data

! BC parameter order 
! cc-surf,din-surf,dip-surf,din-bot,dip-bot,totn-surf,totp-surf,totn-bot,totp-bot,ca-surf

! Loading-file order 
! Last value is suspect? 
! din,dip,totn,totp 
! Note that the layer number is in the attributes (0 surf, 1 bot) 
! The values can be cyclic so that 1 year is used in all

! Initial value order. Note C-style array index 
! [1][0] din surface 
! [2][0] dip surface 
! [5][0] totp surface 
! [6][0] totn surface 
! Column [1] is bottom

! Note that in the original model the algae are not allowed to go into the lower layer. 
! Should a flag be used that either allows or denies this opportunity for the  algae.

PROGRAM wqficos
!------------------------------------------------------------------------- 
! Water quality module based on Kiirikki et al 2001 with some modifications 
! 
! Translated into Fortran90 from Python-code developed by some Swedish dudes. 
! 
! Version history: 0.3 This is a hardcoded version 
! 25.1.2016 0.2 
! 26.1.2016 0.3
! 28.1.2016 0.4 Implicit vertical advection
! 29.1.2016 0.5
! 1.2.2016 0.6 DLSODE
! 5.2.2016 0.7
! 11.2.2016 0.8 Thicknesses written to file
! 24.2.2016 0.9 Thermocline version
! 12.5.2016 0.91 Bug fixes
! 22.11.2017         TIMESTEPMULT set as 0.7 for stability at high resolutions
! 4.12.2017          Cleaning up console writes
! 2019               Major cleanup and refactoring of code
! 2024               (Karel Kaurila, University of Helsinki) Save cumulative N2 fixation to results.
!
! All data is read into arrays at the start of the model run. This might be 
! somewhat memory intensive. 
! 
! TODO What license? 
! 
! Input arguments: 
! (Unused arguments need to be replaced by dummies)
! runid  ! Implemented --> paramsfile
! rundir  ! Write results to rundir as results.hdf5
! hdfile ! Implemented 
! loadfile  ! Implemented
! inifile   ! Location of ini-file
! file_restart  ! Needed as initial conditions (not implemented at this time)
! model_name  ! Implemented
! Optional arguments:
! t0   ! Timestep of start timestep. 
! tfinal  ! Timestep to end simulation
!The model can start from some other timestep than the first one. 
! 
! 
! Kai Rasmus, Luode Consulting Oy, 2015 
! Janne Ropponen, Finnish Environment Institute, 2017-2019
!------------------------------------------------------------------------- 

use hdf5 
use hdfio
use prec
use wqficos_helper ! wqficos_helper.f90
use wq_error ! Error codes and error handling routines 
use wqficos_names ! Model state variables and parameter names

implicit none

external fex_ficos

!!-------------------------------------------- 
!! Part 1 Definitions 
!!--------------------------------------------

!---Version
 character(len=16) :: version="20191118 FICOS"

!---Debug and feature flags (1=on, 0=off for integer flags)
 logical, parameter :: enableoutput = .TRUE.
 integer :: debugadvection = 1 ! 0/! Disable/Enable advection
 integer :: debugintegration = 1 ! 0/1 Disable/Enable integration
 logical, parameter :: enableboundaryconditions = .TRUE.
 logical, parameter :: explicitverticaladvection = .FALSE. ! Explicit vertical advection
 integer, parameter :: debugloading = 0 ! Explicit loading
 logical, parameter :: asciioutput = .FALSE.
 logical, parameter :: enablevvmix = .TRUE. ! Virtual vertical mixing
 logical, parameter :: debugmode = .FALSE. ! Enable debug output in e.g. equations

 logical :: debugrk4 = .FALSE. ! Enable debug in integration (flag value can be determined within run)
 logical :: saveifacefluxes = .FALSE. ! Whether to calculate total nutrient fluxes over user selectable boundaries, detected by checking if fluxfile exists

 integer :: cline=0 ! Depth of constant cline (0 means that the cline is calculated)


 !integer, parameter :: nvars = 8 ! Number of water quality model variables

 integer, parameter :: nvars = 9 ! Number of water quality model variables

!! Define variables
 integer, parameter :: MAX_IFACES=1000 ! maximum number of interfaces in interface list for a block
 real(kind=dp), parameter :: TIMESTEPMULT = 0.7 ! Multiply theoretical maximum timestep by this value for stability
 real(kind=dp), parameter :: TIMESTEPEPSILON = 1.0E-9 ! Epsilon in time stepping calculations
 integer :: inarg ! Number of command line arguments

 integer :: ii, jj, kk, nn, mm, vv ! Loop counters (x,y,z)

! Time variables
 integer :: tt ! Time loop counter
 integer :: ttb,tte ! Begin and end timestep
 integer :: ntims ! Number of timesteps
 real(kind=dp), allocatable, dimension(:) :: deltat ! Timesteps
 real(kind=dp) :: deltatt
 real(kind=dp) :: tottime=0 ! Modeltime elapsed since start of run
 real(kind=dp) :: tottime2
 integer :: starttt=0,stoptt=0 ! Start and stop timesteps
 character(len=6) :: ctt !
 
!! Allocatable result arrays (t,x,z) x is number of blocks 
!! TODO: note that z is 2  at this point
! real(kind=dp),allocatable,dimension(:,:,:) :: s_cc ! Concentration of cyanobacteria [g m-2]
! real(kind=dp),allocatable,dimension(:,:,:) :: s_ca ! Concentration of other algae [g m-2]
! real(kind=dp),allocatable,dimension(:,:,:) :: s_cdin ! Concentration of dissolved nitrogen [mg m-3]
! real(kind=dp),allocatable,dimension(:,:,:) :: s_cdip ! Concentration of dissolved phosphorous [mg m-3]
! real(kind=dp),allocatable,dimension(:,:,:) :: s_cndet ! Concentration of detritus nitrogen [mg m-3]
! real(kind=dp),allocatable,dimension(:,:,:) :: s_cpdet ! Concentration of detritus phosphorous [mg m-3]
! real(kind=dp),allocatable,dimension(:,:,:) :: s_ctotn ! Concentration of total nitrogen
! real(kind=dp),allocatable,dimension(:,:,:) :: s_ctotp ! Concentration of total phosphorous
 real(kind=dp),allocatable,dimension(:,:,:,:) :: s_c ! Concentrations
 
! real(kind=dp),allocatable,dimension(:,:) :: s_tcc ! Concentration of cyanobacteria [g m-2]
! real(kind=dp),allocatable,dimension(:,:) :: s_tca ! Concentration of other algae [g m-2]
! real(kind=dp),allocatable,dimension(:,:) :: s_tcdin ! Concentration of dissolved nitrogen [mg m-3]
! real(kind=dp),allocatable,dimension(:,:) :: s_tcdip ! Concentration of dissolved phosphorous [mg m-3]
! real(kind=dp),allocatable,dimension(:,:) :: s_tcndet ! Concentration of detritus nitrogen [mg m-3]
! real(kind=dp),allocatable,dimension(:,:) :: s_tcpdet ! Concentration of detritus phosphorous [mg m-3]
! real(kind=dp),allocatable,dimension(:,:) :: s_tctotn ! Concentration of total nitrogen
! real(kind=dp),allocatable,dimension(:,:) :: s_tctotp ! Concentration of total phosphorous

 real(kind=dp),allocatable,dimension(:,:,:) :: s_dc ! Concentration delta within a timestemp

 integer(kind=dp),allocatable,dimension(:) :: fluxblocks ! Blocks into/out of which fluxes are saved, Used if saveifacefluxes == .TRUE.
 real(kind=dp),allocatable,dimension(:,:,:,:,:) :: ifacefluxes ! Used if saveifacefluxes == .TRUE.
 real(kind=dp),allocatable,dimension(:,:) :: fluxin, fluxout ! Used if saveifacefluxes == .TRUE.

 real(kind=dp),allocatable,dimension(:,:) :: voldif ! Block layer volume differences after horizontal transport

 integer :: nlayers ! Number of layers (2 in this example) ! Forcing arrays
 real(kind=dp),allocatable,dimension(:) :: f_itot ! Total solar radiation
 real(kind=dp),allocatable,dimension(:,:,:) :: f_avr ! Right advection at interface
 real(kind=dp),allocatable,dimension(:,:,:) :: f_avl ! Left advection at interface
 real(kind=dp),allocatable,dimension(:,:,:) :: f_diff ! Diffusion

 integer,allocatable,dimension(:,:) :: ifaces_right ! Indii of right interfaces
 integer,allocatable,dimension(:,:) :: ifaces_left ! Indii of left interfaces
 integer,allocatable,dimension(:) :: n_ifaces_right ! Number of right interfaces
 integer,allocatable,dimension(:) :: n_ifaces_left ! Number of left interfaces
!! real(kind=dp),allocatable,dimension(:,:,:) :: iface_area ! Area of interface. Not used anywhere/commented out 2019-08-28/JR

! Boundary condition arrays
 real(kind=dp),allocatable,dimension(:,:,:) :: bnddata ! Boundary condition data (nbndblocks,ntims,nparameters=11)
 integer,allocatable,dimension(:) :: bndblocks  ! List of boundary blocks (nbndblocks) 
 integer :: nbndblocks    ! Number of boundary blocks 
 integer,allocatable,dimension(:) ::bndmap ! List of boundary block indii in whole block list
 integer :: index ! Output index
 logical :: bndblockleft, bndblockright
 logical :: isBoundblock
 
!  logical, intent(out) :: error
  logical :: error
  CHARACTER value*(40)
  
  character(128) :: sim_name,paramsfile ! More length to parameters 2017-09-04 JR

!
! HDF-related variables
  character(len=128) :: rundir ! Run directory
  CHARACTER(LEN=128) :: outfilename ! = "results.hdf5"
  character(len=128) :: infilename  !="../test/data/hd_files.hdf5"
  character(len=128) :: loadingfilename !="../test/data/loading.hdf5"
  CHARACTER(LEN=3) , PARAMETER :: dataset = "DS1"
  INTEGER          , PARAMETER :: dim0     = 4
  INTEGER          , PARAMETER :: dim1     = 7

  INTEGER(SIZE_T)  , PARAMETER :: sdim      = 10

  INTEGER :: hdferr

! File I/O
  character(len=128) :: inifilename
  character(len=11)  :: logfile="logfile.txt"
  character(len=128) :: fluxfile = "fluxboundaryblocks.dat" ! Used if saveborderfluxes == .TRUE.

! Handles
  INTEGER(HID_T) :: infile,outfile,loadingfile

! Additional variables
  integer :: blidint  ! block id integer

! Variables used in calculation of advection and diffusion
  integer :: rb,lb,rb_i,lb_i ! Indii of blocks
  real(kind=dp),allocatable,dimension(:,:) :: delta_adv ! Advection changes

! Buffers for variables
  character(len=10),allocatable,dimension(:), target :: block_list
  integer :: nblocks ! Number of blocks
  integer :: nifaces ! Number of interfaces
  integer, allocatable, dimension(:,:), target :: inter_list ! List of interfaces between blocks
  character(len=21), allocatable, dimension(:), target :: s_inter_list ! List of interfaces in string format
  character(len=10) :: init_datum ! Simulation initilisation date
  integer, dimension(6) :: start_date
  integer,dimension(6) :: model_date ! Continuous model date
  integer,allocatable,dimension(:,:) :: iface_block_list

! Block properties
  real(kind=dp),allocatable,dimension(:,:,:) :: b_volume ! Volume of block
  real(kind=dp),allocatable,dimension(:,:,:) :: b_area  ! Area of block (Volume and area can change with time)
  real(kind=dp),allocatable,dimension(:,:) :: b_ad  ! Diffusion within block (TODO: vertical diffusion??)
  real(kind=dp),allocatable,dimension(:,:) :: b_awdw  ! Advection down
  real(kind=dp),allocatable,dimension(:,:) :: b_awup  ! Advection up
  real(kind=dp),allocatable,dimension(:,:,:) :: b_s  ! Salinity
  real(kind=dp),allocatable,dimension(:,:,:) :: b_t  ! Temperature
  real(kind=dp),allocatable,dimension(:,:) :: b_tbot  ! Bottom temperature
  real(kind=dp),allocatable,dimension(:,:) :: b_ubot  ! Bottom velocity
  real(kind=dp),allocatable,dimension(:,:,:) :: b_thick  ! Thickness of block (volume/area)

! Loadings
  real(kind=dp),allocatable,dimension(:,:,:,:) :: loadings  ! Loading (ntims,nblocks,nlayers,nparameters (=4))

! ODE solver equation specific
  integer :: neq=38   ! Number of equations
  real(kind=dp),dimension(38) :: y, y2 ! Equations

!! INITIALIZATION
IF (debugmode) THEN
   ! Remove old debugfile if it exists
   OPEN(unit=199,iostat=ii,file="./ficos_debug_rk4.txt",status='OLD')
   IF (ii .EQ. 0) CLOSE(199, status='delete')
END IF

!!-----------------------------------------------------
!! Part 1.2 Process command line arguments and ini-file
!!-----------------------------------------------------
!!
inarg = command_argument_count()
IF (inarg .LT. 8 .OR. inarg .GT. 8) THEN
   call wqkiirikki_error(error_num_input_args)
END IF

write(*,*) "Number of input arguments: ", inarg

call get_command_argument(1,paramsfile)  ! parameter file name
call get_command_argument(2,rundir)  ! Run directory
call get_command_argument(3,infilename)  ! Inputfilename
call get_command_argument(4,loadingfilename) ! Loading filename
call get_command_argument(5,inifilename) ! Ini file
call get_command_argument(6,sim_name)  ! Simulation name
call get_command_argument(7,ctt)            ! Start timestep
read(ctt,*)starttt
call get_command_argument(8,ctt)            ! Stop timestep
read(ctt,*)stoptt

! Outputfile name
outfilename = trim(rundir)//"/results.hdf5"
write (*,*) "Run directory: ", rundir
call system("mkdir -p "//TRIM(rundir)) ! Make output directory if necessary
write (*,*) "Writing results to file: ", outfilename

!! Initialize parameter values
write (*,*) "Wqficos ini-file: ", inifilename
call initvalues(inifilename,paramsfile)
call listparams

nlayers = 2 ! Number of layers
write(*,*) "Number of layers: ", nlayers

!---HDF5-related stuff
! Initialize Fortran-interface
CALL h5open_f(hdferr)

! Open input hydrodynamics file
CALL h5fopen_f(infilename, H5F_ACC_RDWR_F,infile, hdferr)

! Open results file (output)
CALL h5fcreate_f(outfilename, H5F_ACC_TRUNC_F, outfile, hdferr)

! Open input loading file
CALL h5fopen_f(loadingfilename, H5F_ACC_RDWR_F, loadingfile, hdferr)

!---Read parameters from ini-file wqficos.ini or what is set in inifilename
call setIniFilename(inifilename)
call loadOptions

call getValue('model','layers', value, error)
if(error.eqv..FALSE.)read(value,'(I3)')nlayers

call getValue('model','neq', value, error)
if(error.eqv..FALSE.)read(value,'(I3)')neq

call getValue('model','advection', value, error)
if(error.eqv..FALSE.)then
!  write(*,*)value
  read(value,'(I1)')debugadvection
  if(debugadvection.lt.0.and.debugadvection.gt.1)then
    call wqkiirikki_error(error_num_ini_advection,debugadvection)
    debugadvection=1
  end if
end if

call getValue('model','integration', value, error)
if(error.eqv..FALSE.)then
!  write(*,*)value
  read(value,'(I1)')debugintegration
  if(debugintegration.lt.0.and.debugintegration.gt.1)then
    call wqkiirikki_error(error_num_ini_integration,debugintegration)
    debugintegration=1
  end if
end if

! JR note 2019-11-15: Thermocline depth is not currently used anywhere!
call getValue('model','thermoclinedepth', value, error)
IF (error .EQV. .FALSE.) READ(value,'(I3)') cline

!!
!!--------------------------
!! Part 2 Initial file I/O
!!--------------------------

!! Note to self and to anyone else who may be listening. All of the input and forcing data is read into memory
!! at model initialisation. If this is not possible due to memory limitation issues, then I suggest
!! that we write them to fortran binary files that I know how to use. The HDF-interface is a bit cumbersome.
!! The other option is to learn how to use the HDF-interface. This might be useful to do in any case.
!! The HDF-interface isn't exactly rocket science.

!! Initial date
call wqhdf_getinitialdate(init_datum,start_date,infile)
call incrementdate(start_date,0) ! Add day of year value to start_date(6)
write(*,*)"Simulation start date: ",start_date

! Block list
call wqhdf_getnblocks(nblocks,infile)
write(*,*)"Number of blocks: ",nblocks

  ! Allocate memory for block list
allocate(block_list(1:nblocks))

call wqhdf_getblocklist(nblocks,block_list,infile)

! Interfaces
call wqhdf_getnifaces(nifaces,infile)
write(*,*)"Number of interfaces: ",nifaces

  ! Allocate memory for interface lists
ALLOCATE(inter_list(1:2,1:nifaces))
allocate(s_inter_list(1:nifaces))

allocate(ifaces_right(nblocks,MAX_IFACES))
allocate(ifaces_left(nblocks,MAX_IFACES)) 
allocate(n_ifaces_left(nblocks))
allocate(n_ifaces_right(nblocks))

call wqhdf_getifacelist(nifaces,inter_list,s_inter_list,infile)

call getifaceindii(nblocks,nifaces,ifaces_right,ifaces_left,&
        n_ifaces_right,n_ifaces_left,inter_list,block_list,MAX_IFACES)

! Generate interface block-list table
write(*,*)"Generating interface block-list table"
allocate(iface_block_list(nifaces,2))
do nn=1,nifaces
      ! rb right block, lb left block
      lb=inter_list(1,nn) ! block numbers
      rb=inter_list(2,nn)
      call getblockindex(nblocks,block_list,lb,lb_i)
      iface_block_list(nn,1)=lb_i
      call getblockindex(nblocks,block_list,rb,rb_i)
      iface_block_list(nn,2)=rb_i
end do

! Timesteps
call wqhdf_getntimesteps(ntims,nifaces,s_inter_list,infile)

write(*,*)"Number of timesteps: ",ntims

!! Interface areas (not used anywhere! Commented out 2019-08-28/JR)
!!allocate(iface_area(ntims,nifaces,nlayers))
!!
!!! Iface areas returned as zero values for some reason??
!!call wqhdf_getifacearea(nifaces,nlayers,ntims,iface_area,s_inter_list,infile)

! Allocate memory for timesteps
allocate(deltat(1:ntims))

call wqhdf_gettimesteps(ntims,deltat,nifaces,s_inter_list,infile)

! Allocate and get block properties
allocate(b_volume(ntims,nblocks,nlayers+1))! Volume of block ! JR: Check why nlayers+1 is needed!
allocate(b_area(ntims,nblocks,nlayers)) ! Area of block
allocate(b_ad(ntims,nblocks))  ! Diffusion within block (TODO: vertical diffusion??)
allocate(b_awdw(ntims,nblocks)) ! Advection down
allocate(b_awup(ntims,nblocks)) ! Advection up
allocate(b_s(ntims,nblocks,nlayers)) ! Salinity
allocate(b_t(ntims,nblocks,nlayers)) ! Temperature
allocate(b_tbot(ntims,nblocks)) ! Bottom temperature
allocate(b_ubot(ntims,nblocks)) ! Bottom velocity
allocate(b_thick(ntims,nblocks,nlayers)) ! Thickness of block

call wqhdf_getblockproperties(ntims,nblocks,nlayers,b_volume,b_area,&
                              b_ad,b_awdw,b_awup,b_s,b_t,b_tbot,b_ubot,&
                              block_list,infile)

! Calculate block thickness
WHERE (b_area(:,:,:).GT.0.0)
   b_thick(:,:,:) = b_volume(:,:,1:nlayers)/b_area(:,:,:)
ELSEWHERE
   b_thick(:,:,:) = 0.0
END WHERE
 
! Allocate result arrays (t,x,z) for state variables
!allocate(s_cc(ntims,nblocks,nlayers))
!allocate(s_ca(ntims,nblocks,nlayers))
!allocate(s_cdin(ntims,nblocks,nlayers))
!allocate(s_cdip(ntims,nblocks,nlayers))
!allocate(s_cndet(ntims,nblocks,nlayers))
!allocate(s_cpdet(ntims,nblocks,nlayers))
!allocate(s_ctotn(ntims,nblocks,nlayers))
!allocate(s_ctotp(ntims,nblocks,nlayers))
allocate(s_c(ntims,nblocks,nlayers,nvars))

allocate(voldif(nblocks,nlayers))

!allocate(s_tcc(nblocks,nlayers))
!allocate(s_tca(nblocks,nlayers))
!allocate(s_tcdin(nblocks,nlayers))
!allocate(s_tcdip(nblocks,nlayers))
!allocate(s_tcndet(nblocks,nlayers))
!allocate(s_tcpdet(nblocks,nlayers))
!allocate(s_tctotn(nblocks,nlayers))
!allocate(s_tctotp(nblocks,nlayers))

allocate(s_dc(nblocks,nlayers,nvars))

! If fluxfile exists, automatically enable saving of interface fluxes
INQUIRE(FILE=TRIM(fluxfile), EXIST=saveifacefluxes)
IF (saveifacefluxes) THEN
  fluxblocks = getFluxBorders(TRIM(fluxfile)) ! fluxblocks(:) integer array, automatically allocated in Fortran 2003

  IF (SIZE(fluxblocks).LE.0) THEN
     WRITE(*,*) "No blocks found in ", fluxfile, ", disabling flux logging."
     saveifacefluxes = .FALSE.
  ELSE
     ALLOCATE(ifacefluxes(nifaces,ntims,nlayers,2,nvars)) ! interfaces, timeteps, layers, direction, variables
     ifacefluxes = 0.0
  END IF

END IF

!---Initialise arrays
!s_cc     = 0.0; s_ca     = 0.0; s_cdin   = 0.0; s_cdip   = 0.0
!s_cndet  = 0.0; s_cpdet  = 0.0; s_ctotn  = 0.0; s_ctotp  = 0.0
s_c = 0.0
voldif   = 0.0
!s_tcc    = 0.0; s_tca    = 0.0; s_tcdin  = 0.0; s_tcdip  = 0.0
!s_tcndet = 0.0; s_tcpdet = 0.0; s_tctotn = 0.0; s_tctotp = 0.0
s_dc = 0.0

!---Initial conditions
!s_cc(:,:,1) = FCrefuge
!s_ca(:,:,1) = Arefuge
!s_cndet = 3.0 ! obsolete?
!s_cpdet = 0.2 ! obsolete?
!!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_cc,1,block_list,infile)
!!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_ca,2,block_list,infile)
!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_cdin,2,block_list,infile)
!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_cdip,3,block_list,infile)
!!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_cndet,4,block_list,infile)
!!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_cpdet,5,block_list,infile)
!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_ctotn,7,block_list,infile)
!call wqhdf_getinitcond(ntims,nblocks,nlayers,s_ctotp,6,block_list,infile)

s_c(:,:,1,1) = FCrefuge
s_c(:,:,1,2) = Arefuge
s_c(:,:,:,5) = 3.0 ! obsolete?
s_c(:,:,:,6) = 0.2 ! obsolete?
call wqhdf_getinitcond(ntims,nblocks,nlayers,s_c(:,:,:,3),2,block_list,infile)
call wqhdf_getinitcond(ntims,nblocks,nlayers,s_c(:,:,:,4),3,block_list,infile)
call wqhdf_getinitcond(ntims,nblocks,nlayers,s_c(:,:,:,7),7,block_list,infile)
call wqhdf_getinitcond(ntims,nblocks,nlayers,s_c(:,:,:,8),6,block_list,infile)

! Allocate and get forcing arrays
! Radiation
allocate(f_itot(ntims))
call wqhdf_getradiation(ntims,f_itot,infile)

! Advection and diffusion
allocate(f_avl(nifaces,ntims,nlayers))
allocate(f_avr(nifaces,ntims,nlayers))
allocate(f_diff(nifaces,ntims,nlayers))
call wqhdf_getadvection(ntims,nifaces,nlayers,f_avl,f_avr,f_diff,s_inter_list,infile)

! Advection and diffusion changes
allocate(delta_adv(nblocks,nlayers)) 

! Allocate and get boundary conditions
call wqhdf_getnbndblocks(nbndblocks,infile)
write(*,*) "Number of boundary blocks: ", nbndblocks
allocate(bndblocks(nbndblocks))
allocate(bnddata(nbndblocks,ntims,11))
allocate(bndmap(nbndblocks))
call wqhdf_getbnddata(ntims,nbndblocks,11,bndblocks,bnddata,infile)

! find boundary blocks in block list and save their indii
do ii=1,nbndblocks
   call getblockindex(nblocks,block_list,bndblocks(ii),index)
   bndmap(ii)=index
end do

! Create output datasets
if (enableoutput) then
  if(starttt.eq.0.and.stoptt.eq.0)then
    CALL wqhdf_createoutput(nblocks,ntims+1,block_list,outfile)
  else
    CALL wqhdf_createoutput(nblocks,stoptt-starttt+3,block_list,outfile)
  end if
end if

! Allocate and get loadings
allocate(loadings(ntims,nblocks,nlayers,4)) ! 4 parameters
loadings=0.0
call wqhdf_loadings(ntims,nblocks,nlayers,4,block_list,loadings,loadingfile) ! 4 parameters


!!
!!--------------------------
!! Part 3 Time loop
!!--------------------------

write(*,*) "Water quality model code version: ",version

! Initialize model_date
model_date = start_date

! Start from offset start index
ttb=1
if (starttt.ne.0) then
 do ttb=1,starttt
  call incrementdate(model_date,1)
 end do
end if

! Set stop index
tte=ntims
if (stoptt.ne.0) then
 tte=stoptt
end if

! The interface advection and diffusion are resolved with the same timestep
! as the differential equations with regard to time are solved.

! Write initial values
IF (enableoutput) THEN
!  call hdfwq_writeoutput("cDIN_0",s_cdin(1,:,1),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("cDIN_1",s_cdin(1,:,2),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("cDIP_0",s_cdip(1,:,1),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("cDIP_1",s_cdip(1,:,2),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("Ndet_0",s_cndet(1,:,1),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("Ndet_1",s_cndet(1,:,2),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("Pdet_0",s_cpdet(1,:,1),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("Pdet_1",s_cpdet(1,:,2),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("totN_0",s_ctotn(1,:,1),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("totN_1",s_ctotn(1,:,2),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("totP_0",s_ctotp(1,:,1),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("totP_1",s_ctotp(1,:,2),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("cC",s_cc(1,:,1),nblocks,block_list,1,outfile)
!  call hdfwq_writeoutput("cA",s_ca(1,:,1),nblocks,block_list,1,outfile)

  call hdfwq_writeoutput("cC",s_c(1,:,1,1),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("cA",s_c(1,:,1,2),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("cDIN_0",s_c(1,:,1,3),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("cDIN_1",s_c(1,:,2,3),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("cDIP_0",s_c(1,:,1,4),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("cDIP_1",s_c(1,:,2,4),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("Ndet_0",s_c(1,:,1,5),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("Ndet_1",s_c(1,:,2,5),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("Pdet_0",s_c(1,:,1,6),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("Pdet_1",s_c(1,:,2,6),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("totN_0",s_c(1,:,1,7),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("totN_1",s_c(1,:,2,7),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("totP_0",s_c(1,:,1,8),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("totP_1",s_c(1,:,2,8),nblocks,block_list,1,outfile)
  ! ---------------------- KK 16.04.24 save cumulative N2fix
  call hdfwq_writeoutput("rN2fixFC",s_c(1,:,1,9),nblocks,block_list,1,outfile)
  call hdfwq_writeoutput("cN2fixFC",s_c(1,:,2,9),nblocks,block_list,1,outfile)
  ! -------------------------------------------------------------------
  
  ! Update date and write it to file
  call hdfwq_writedate(model_date,1,outfile)

  ! Write initial thicknesses
  call hdfwq_writethickness(nblocks,b_thick(1,:,1),b_thick(1,:,2),&
                          block_list,1,outfile)
END IF

! Note that all data is provided with a 24h timestep

!=================
! TIME LOOP: DAILY
!=================
Dayloop: DO tt=ttb,tte ! Usually 1,ntims

!if (tt.ge.3) STOP !JR DEBUG

!--- Copy old values to current values
  IF (tt.GT.1) THEN
!     s_cc(tt,:,:)    = s_cc(tt-1,:,:)
!     s_ca(tt,:,:)    = s_ca(tt-1,:,:)
!     s_cdin(tt,:,:)  = s_cdin(tt-1,:,:)
!     s_cdip(tt,:,:)  = s_cdip(tt-1,:,:)
!     s_cndet(tt,:,:) = s_cndet(tt-1,:,:)
!     s_cpdet(tt,:,:) = s_cpdet(tt-1,:,:)
!     s_ctotn(tt,:,:) = s_ctotn(tt-1,:,:)
!     s_ctotp(tt,:,:) = s_ctotp(tt-1,:,:)

     s_c(tt,:,:,:) = s_c(tt-1,:,:,:)
! ----- (KK 08.08.2024) ------
! Initialize N2fix values 
! to ensure stored values are daily cumulants 
     s_c(tt,:,:,9) = 0.0
! -----------------------
  END IF
!---Multiply maximum inner loop timestep deltatt by TIMESTEPMULT (0..1) for stability
  deltatt = deltat(tt)*TIMESTEPMULT

!---Console logging messages
  write(*,*) ""
  write(*,'(A4,I6,A4,I6,A30,F7.1)') "Day ",tt," of ",ntims," with intra-day timestep [s]: ", deltatt
  write(*,'(A40,4F10.3)') "TOTAL LOADING [TONS] DIN,DIP,TOTN,TOTP: ", &
       sum(loadings(tt,1:nblocks,1,1)+loadings(tt,1:nblocks,2,1))/1000.0, &
       sum(loadings(tt,1:nblocks,1,2)+loadings(tt,1:nblocks,2,2))/1000.0, &
       sum(loadings(tt,1:nblocks,1,3)+loadings(tt,1:nblocks,2,3))/1000.0, &
       sum(loadings(tt,1:nblocks,1,4)+loadings(tt,1:nblocks,2,4))/1000.0

!---Log file
   open(unit=10112,file=logfile)
   write(10112,'(AI5AI5AI5)')"Start timestep: ",ttb," Stop timestep: ",tte," Current timestep: ",tt
   close(10112)

!=====================
! TIME LOOP: Intra-day
!=====================

  tottime2=0.0 ! Elapsed time within current day
  Intradayloop: DO WHILE (tottime2.LT.86400.0)
     deltatt = MIN(deltatt,86400.0-tottime2+TIMESTEPEPSILON)
     tottime2 = tottime2 + deltatt

!Update global time
     tottime=tottime+deltatt

!!
!! Update boundary conditions
!!
  IF (enableboundaryconditions) THEN
    DO mm=1,nbndblocks
        ! algal_N [ug/L] = algal_ww [ug/L] * 0.11 * 16*14/(106*12) = algal_ww [ug/L) * 0.019371 = algal_ww [g/m3] * 19.371
!        s_cc(tt,bndmap(mm),1)    = bnddata(mm,tt,1)*19.371 ! ww g/m3 -> mg N/m3 ! JR 2017-10-25
!        s_cdin(tt,bndmap(mm),1)  = bnddata(mm,tt,2)
!        s_cdip(tt,bndmap(mm),1)  = bnddata(mm,tt,3)
!        s_cdin(tt,bndmap(mm),2)  = bnddata(mm,tt,5)
!        s_cdip(tt,bndmap(mm),2)  = bnddata(mm,tt,6)
!
!        s_ctotn(tt,bndmap(mm),1) = bnddata(mm,tt,7)
!        s_ctotp(tt,bndmap(mm),1) = bnddata(mm,tt,8)
!        s_ctotn(tt,bndmap(mm),2) = bnddata(mm,tt,9)
!        s_ctotp(tt,bndmap(mm),2) = bnddata(mm,tt,10)
!        s_ca(tt,bndmap(mm),1)    = bnddata(mm,tt,11)*19.371 ! ww g/m3 -> mg N/m3 ! JR 2017-10-25

        IF (bndmap(mm).EQ.0) THEN
           WRITE(*,'(A,I10,A)',advance="no") "WARNING: Skipping block ", bndblocks(mm), " which has been defined"
           WRITE(*,*) "as boundary block but does not exist."
           CYCLE
        END IF

        s_c(tt,bndmap(mm),1,1)    = bnddata(mm,tt,1)*19.371 ! ww g/m3 -> mg N/m3 ! JR 2017-10-25
        s_c(tt,bndmap(mm),1,3)  = bnddata(mm,tt,2)
        s_c(tt,bndmap(mm),1,4)  = bnddata(mm,tt,3)
        s_c(tt,bndmap(mm),2,3)  = bnddata(mm,tt,5)
        s_c(tt,bndmap(mm),2,4)  = bnddata(mm,tt,6)

        s_c(tt,bndmap(mm),1,7) = bnddata(mm,tt,7)
        s_c(tt,bndmap(mm),1,8) = bnddata(mm,tt,8)
        s_c(tt,bndmap(mm),2,7) = bnddata(mm,tt,9)
        s_c(tt,bndmap(mm),2,8) = bnddata(mm,tt,10)
        s_c(tt,bndmap(mm),1,2)    = bnddata(mm,tt,11)*19.371 ! ww g/m3 -> mg N/m3 ! JR 2017-10-25
    END DO
  END IF

  do ii=1,nblocks
    ! Check if this is a boundary block
    isBoundblock = .FALSE.
    READ(block_list(ii),*) blidint ! Read block id character string into an integer variable
    IF (ANY(bndblocks.EQ.blidint)) isBoundblock = .TRUE.

    !Point source [kg/d]
    if (debugloading.eq.1 .AND. .NOT. isBoundBlock) then
      do jj=1,nlayers
        IF (b_volume(tt,ii,jj).gt.0.0) THEN
          ! loadings are in kg/d, internal values in mg/m3, deltatt in s
!          s_cdin(tt,ii,jj)=s_cdin(tt,ii,jj)+&
!                ((1.0E6/86400.0)*loadings(tt,ii,jj,1)/b_volume(tt,ii,jj))*deltatt
!          s_cdip(tt,ii,jj)=s_cdip(tt,ii,jj)+&
!                ((1.0E6/86400.0)*loadings(tt,ii,jj,2)/b_volume(tt,ii,jj))*deltatt
!          s_ctotn(tt,ii,jj)=s_ctotn(tt,ii,jj)+&
!                ((1.0E6/86400.0)*loadings(tt,ii,jj,3)/b_volume(tt,ii,jj))*deltatt
!          s_ctotp(tt,ii,jj)=s_ctotp(tt,ii,jj)+&
!                ((1.0E6/86400.0)*loadings(tt,ii,jj,4)/b_volume(tt,ii,jj))*deltatt
          s_c(tt,ii,jj,3)=s_c(tt,ii,jj,3)+&
                ((1.0E6/86400.0)*loadings(tt,ii,jj,1)/b_volume(tt,ii,jj))*deltatt
          s_c(tt,ii,jj,4)=s_c(tt,ii,jj,4)+&
                ((1.0E6/86400.0)*loadings(tt,ii,jj,2)/b_volume(tt,ii,jj))*deltatt
          s_c(tt,ii,jj,7)=s_c(tt,ii,jj,7)+&
                ((1.0E6/86400.0)*loadings(tt,ii,jj,3)/b_volume(tt,ii,jj))*deltatt
          s_c(tt,ii,jj,8)=s_c(tt,ii,jj,8)+&
                ((1.0E6/86400.0)*loadings(tt,ii,jj,4)/b_volume(tt,ii,jj))*deltatt
        END IF
      end do
    end if
  end do


!------------------------
! Advection and diffusion
!------------------------

  !---Horizontal advection
  hadv: IF (debugadvection.eq.1) THEN
    ! Temporary concentration delta arrays for current interface
!    s_tcc    = 0.0 !s_cc
!    s_tca    = 0.0 !s_ca
!    s_tcdin  = 0.0 !s_cdin
!    s_tcdip  = 0.0 !s_cdip
!    s_tcndet = 0.0 !s_cndet
!    s_tcpdet = 0.0 !s_cpdet
!    s_tctotn = 0.0 !s_ctotn
!    s_tctotp = 0.0 !s_ctotp
    s_dc = 0.0

    voldif = 0.0

    ifloop: DO nn=1,nifaces ! Loop over interfaces
      ! rb right block, lb left block
      lb = inter_list(1,nn) ! Left block id in interface
      rb = inter_list(2,nn) ! Right block id in interface

      lb_i = iface_block_list(nn,1) ! Left block index
      rb_i = iface_block_list(nn,2) ! Right block index

      ! Check if either block is a boundary block.
      ! Boundary block values are not updated in the advection routines.
 !     bndblockright=0
!      READ(block_list(rb_i),*) blidint ! Read block id character string into an integer variable
!      IF (ANY(bndblocks.EQ.blidint)) bndblockright=1
!
!      bndblockleft=0
!      READ(block_list(lb_i),*) blidint ! Read block id character string into an integer variable
!      IF (ANY(bndblocks.EQ.blidint)) bndblockleft=1
      bndblockleft = .FALSE.
      bndblockright = .FALSE.
      IF (ANY(bndblocks.EQ.lb)) bndblockleft = .TRUE.
      IF (ANY(bndblocks.EQ.rb)) bndblockright = .TRUE.

      ! Do not update algae over interfaces with boundary blocks
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_cc,s_tcc,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)
!
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_ca,s_tca,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)
!
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_cdin,s_tcdin,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)
!
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_cdip,s_tcdip,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)
!
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_cndet,s_tcndet,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)
!
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_cpdet,s_tcpdet,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)
!
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_ctotn,s_tctotn,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)
!
!      call calcadv(ntims,nlayers,nblocks,nifaces,&
!                s_ctotp,s_tctotp,f_avl,f_avr,f_diff,deltatt,&
!                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
!                bndblockleft,bndblockright)

      call calcadv(ntims,nlayers,nblocks,nifaces,nvars,&
                s_c,s_dc,f_avl,f_avr,f_diff,deltatt,&
                tt,rb_i,lb_i,nn,b_volume(tt,lb_i,1:nlayers),b_volume(tt,rb_i,1:nlayers),&
                bndblockleft,bndblockright)

      IF (saveifacefluxes) THEN
        !---Save mass fluxes to both directions over an interface to ifacefluxes.
        !---We combine fluxes from each intra-day timestep.
        varloop: DO vv = 1,nvars
          ! Left block --> right block [kg]
          ifacefluxes(nn,tt,:,1,vv) = ifacefluxes(nn,tt,:,1,vv) &
                                   + s_c(tt,lb_i,:,vv) * (-f_avr(nn,tt,:)) * deltatt
          ! Left block <-- right block [kg]
          ifacefluxes(nn,tt,:,2,vv) = ifacefluxes(nn,tt,:,2,vv) &
                                   + s_c(tt,rb_i,:,vv) * f_avl(nn,tt,:) * deltatt
        END DO varloop

      END IF

      ! Calculates net volume difference (voldif) in block layers introduced by horizontal advection
      IF (enablevvmix) THEN
         call calcvoldif(ntims,nlayers,nblocks,nifaces,&
                   f_avl,f_avr,voldif,deltatt,tt,&
                   rb_i,lb_i,nn,bndblockleft,bndblockright)
      END IF

    END DO ifloop ! end loop over interfaces

  !---Update concentrations with changes from this timestep
!  s_cdin(tt,:,:)  = s_cdin(tt,:,:)  + s_tcdin
!  s_cdip(tt,:,:)  = s_cdip(tt,:,:)  + s_tcdip
!  s_cc(tt,:,:)    = s_cc(tt,:,:)    + s_tcc
!  s_ca(tt,:,:)    = s_ca(tt,:,:)    + s_tca
!  s_cndet(tt,:,:) = s_cndet(tt,:,:) + s_tcndet
!  s_cpdet(tt,:,:) = s_cpdet(tt,:,:) + s_tcpdet
!  s_ctotn(tt,:,:) = s_ctotn(tt,:,:) + s_tctotn ! Comment out for vvmix2, vvmix3, vvmix4
!  s_ctotp(tt,:,:) = s_ctotp(tt,:,:) + s_tctotp

  s_c(tt,:,:,:) = s_c(tt,:,:,:) + s_dc

  !---Virtual vertical mixing [vvmix 1b: after adding mass transports]
  IF (enablevvmix) THEN
     !CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_cc, b_volume, voldif)
     !CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_ca, b_volume, voldif)
     !Do not mix algae!
!     CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_cdin, b_volume, voldif)
!     CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_cdip, b_volume, voldif)
!     CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_cndet, b_volume, voldif)
!     CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_cpdet, b_volume, voldif)
!     CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_ctotn, b_volume, voldif)
!     CALL virtualverticalmixing(ntims, nblocks, nlayers, tt, s_ctotp, b_volume, voldif)

     CALL virtualverticalmixing(ntims, nblocks, nlayers, nvars, tt, 3, 8, s_c, b_volume, voldif)
  END IF

  END IF hadv

!---Explicit vertical advection and diffusion
!---Note that the volume has not been corrected here for time-varying volume
!---Normally vertical advection is implicit during integration
explicit_vadv: IF (explicitverticaladvection) THEN
!   ! Algae are not allowed to go to the bottom layer in this implementation
!   !call calcvertadv(ntims,nlayers,nblocks,&
!   !                s_cc,deltatt,tt,&
!   !                b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
!   !call calcvertadv(ntims,nlayers,nblocks,&
!   !                s_ca,deltatt,tt,&
!   !                b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
!   call calcvertadv(ntims,nlayers,nblocks,&
!                   s_cdin,deltatt,tt,&
!                   b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
!   call calcvertadv(ntims,nlayers,nblocks,&
!                   s_cdip,deltatt,tt,&
!                   b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
!   call calcvertadv(ntims,nlayers,nblocks,&
!                   s_cndet,deltatt,tt,&
!                   b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
!   call calcvertadv(ntims,nlayers,nblocks,&
!                   s_cpdet,deltatt,tt,&
!                   b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
!   call calcvertadv(ntims,nlayers,nblocks,&
!                   s_ctotn,deltatt,tt,&
!                   b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
!   call calcvertadv(ntims,nlayers,nblocks,&
!                   s_ctotp,deltatt,tt,&
!                   b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)

   call calcvertadv(ntims,nlayers,nblocks,nvars,3,8,s_c,deltatt,tt,&
                   b_volume(:,:,1:nlayers),b_ad,b_awdw,b_awup)
END IF explicit_vadv


!---------------------
!Integration over time
!---------------------

BlockLoop: DO mm=1,nblocks

 ! Check if this is a boundary block
 isBoundBlock = .FALSE.
 READ(block_list(mm),*) blidint ! Read block id character string into an integer variable
 IF (ANY(bndblocks.EQ.blidint)) isBoundblock = .TRUE.

 integrate: IF (debugintegration.eq.1 .AND. .NOT. isBoundblock) THEN
   y=0.0
   y(1) = s_c(tt,mm,1,1) ! cC 1
   y(2) = s_c(tt,mm,1,2) ! cA 1

   y(3) = s_c(tt,mm,1,3) ! cdin 1
   y(4) = s_c(tt,mm,1,4) ! cdip 1
   y(5) = s_c(tt,mm,2,3) ! cdin 2
   y(6) = s_c(tt,mm,2,4) ! cdip 2

   y(7) = s_c(tt,mm,1,5) ! detN 1
   y(8) = s_c(tt,mm,1,6) ! detP 1
   y(9) = s_c(tt,mm,2,5) ! detN 2
   y(10) = s_c(tt,mm,2,6) ! detP 2

   y(11) = s_c(tt,mm,1,7) ! totN 1
   y(12) = s_c(tt,mm,1,8) ! totP 1
   y(13) = s_c(tt,mm,2,7) ! totN 2
   y(14) = s_c(tt,mm,2,8) ! totP 2
   
   !---------------------------------------------
   ! N2 fixation stored in spare vars (KK 2024-08-07)
   y(34) = s_c(tt,mm,1,9) 
   y(35) = s_c(tt,mm,2,9)
   !----------------------------------------------


!   y(1)=s_cc(tt,mm,1)
!   y(2)=s_ca(tt,mm,1)
!   y(3)=s_cdin(tt,mm,1)
!   y(4)=s_cdip(tt,mm,1)
!   y(5)=s_cdin(tt,mm,2)
!   y(6)=s_cdip(tt,mm,2)
!   y(7)=s_cndet(tt,mm,1)
!   y(8)=s_cpdet(tt,mm,1)
!   y(9)=s_cndet(tt,mm,2)
!   y(10)=s_cpdet(tt,mm,2)
!   y(11)=s_ctotn(tt,mm,1)
!   y(12)=s_ctotp(tt,mm,1)
!   y(13)=s_ctotn(tt,mm,2)
!   y(14)=s_ctotp(tt,mm,2)

   y(15)=b_t(tt,mm,1)
   y(16)=b_t(tt,mm,2)
   y(17)=f_itot(tt)

   IF (b_t(tt,mm,1).le.0.0) THEN
        y(18)=1.0 ! ice if temperature equal to zero or below
   ELSE
        y(18)=0.0 ! No ice otherwise
   END IF

   y(19)=b_thick(tt,mm,1)
   y(20)=b_thick(tt,mm,2)
   y(21)=0.0D0 ! Spare var, was: real(jj,dp) ! Cast layer number to real(kind=dp)
   y(22)=b_awdw(tt,mm)*86400.0 ! Vertical advection [m/d]
   y(23)=b_awup(tt,mm)*86400.0
   y(24)=b_volume(tt,mm,1)
   y(25)=b_volume(tt,mm,2)

   IF (y(24).eq.0.0) THEN
      WRITE(*,*) "FATAL ERROR: Layer 1 volume equal to 0 in block: ", block_list(mm)
      STOP
   END IF

   y(38)=b_ubot(tt,mm)

   ! Loadings [mg/d]
   y(26)=1.0D6*loadings(tt,mm,1,1) ! din layer 1
   y(27)=1.0D6*loadings(tt,mm,2,1) ! din layer 2
   y(28)=1.0D6*loadings(tt,mm,1,2) ! dip layer 1
   y(29)=1.0D6*loadings(tt,mm,2,2) ! dip layer 2
   y(30)=1.0D6*loadings(tt,mm,1,3) ! totn layer 1
   y(31)=1.0D6*loadings(tt,mm,2,3) ! totn layer 2
   y(32)=1.0D6*loadings(tt,mm,1,4) ! totp layer 1
   y(33)=1.0D6*loadings(tt,mm,2,4) ! totp layer 2

   ! Spare vars

   y(36) = 0.0D0
   y(37) = REAL(model_date(6),kind=dp) ! Day of year

   !==================================
   ! INTEGRATE WATER QUALITY EQUATIONS
   !==================================

   y2 = 0.0

   !---Debug/check for errors before integration

   ! Layer 1 thickness & volume check (debug)
    IF ( y(19) .LE. 0.0 .OR. y(24) .LE. 0.0) THEN
       WRITE(*,*) "FATAL Epic fail: Surface layer thickness cannot be zero."
       CALL EXIT()
    END IF

   ! Layer 2 thickness & volume check (debug)
   IF ( (y(20) .GT. 0.0 .AND. y(25) .LE. 0.0) .OR. (y(20) .LE. 0.0 .AND. y(25) .GT. 0.0) ) THEN
       WRITE(*,*) "FATAL Epic fail: Bottom layer thickness does not match bottom volume"
       CALL EXIT()
    END IF

   ! Make sure that concentrations are positive before integration
   IF (ANY(y(1:14).LT.0.0)) THEN
      WRITE(*,*) "Warning: removing negative y concentrations before integration in block ", &
         block_list(mm)
      WRITE(*,*) "Offending values: ", y(1:14)
      WHERE (y(1:14).LT.0.0) y(1:14) = 0.0
   END IF

   !---INTEGRATE: Fourth order Runge-Kutta (RK4) call (explicit method)
   debugrk4 = debugmode .AND. block_list(mm).EQ."1900000099"
   CALL rk4vec(0.0,neq,y,deltatt/86400.0,fex_ficos,y2,debugrk4)

   ! Make sure that concentrations are positive after integration
   WHERE (y2(1:14).LT.0.0) y2(1:14) = 0.0

   !---Update variables
!   s_cc(tt,mm,1) = max(FCrefuge,y2(1))
!   s_ca(tt,mm,1) = max(Arefuge,y2(2))
!
!   s_cdin(tt,mm,1) = y2(3)
!   s_cdip(tt,mm,1) = y2(4)
!   s_cdin(tt,mm,2) = y2(5)
!   s_cdip(tt,mm,2) = y2(6)
!
!   s_cndet(tt,mm,1) = y2(7)
!   s_cpdet(tt,mm,1) = y2(8)
!   s_cndet(tt,mm,2) = y2(9)
!   s_cpdet(tt,mm,2) = y2(10)
! 
!   s_ctotn(tt,mm,1) = y2(11)
!   s_ctotp(tt,mm,1) = y2(12)
!   s_ctotn(tt,mm,2) = y2(13)
!   s_ctotp(tt,mm,2) = y2(14)

   s_c(tt,mm,1,1) = max(FCrefuge,y2(1)) ! cC
   s_c(tt,mm,1,2) = max(Arefuge,y2(2)) ! cA

   s_c(tt,mm,1,3) = y2(3) ! din 1
   s_c(tt,mm,1,4) = y2(4) ! dip 1
   s_c(tt,mm,2,3) = y2(5) ! din 2
   s_c(tt,mm,2,4) = y2(6) ! dip 2

   s_c(tt,mm,1,5) = y2(7) ! detN 1
   s_c(tt,mm,1,6) = y2(8) ! detP 1
   s_c(tt,mm,2,5) = y2(9) ! detN 2
   s_c(tt,mm,2,6) = y2(10) ! detP 2

   s_c(tt,mm,1,7) = y2(11) ! totN 1
   s_c(tt,mm,1,8) = y2(12) ! totP 1
   s_c(tt,mm,2,7) = y2(13) ! totN 2
   s_c(tt,mm,2,8) = y2(14) ! totP 2
   
   ! ---- KK 16.04.24 save cumulative N2fix
   s_c(tt,mm,1,9) = y2(34) ! N2fix [mg]
   s_c(tt,mm,2,9) = y2(35) ! N2fix [mg N/m3]=[mu g N L-1]
   ! ------------------------------------------------
 END IF integrate

END DO BlockLoop ! Block loop

END DO IntradayLoop ! Inner timeloop

!!
!! Write data to file every 24h
!!
IF (enableoutput) THEN
!   call hdfwq_writeoutput("cC",s_cc(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("cA",s_ca(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("cDIN_0",s_cdin(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("cDIN_1",s_cdin(tt,:,2),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("cDIP_0",s_cdip(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("cDIP_1",s_cdip(tt,:,2),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("Ndet_0",s_cndet(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("Ndet_1",s_cndet(tt,:,2),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("Pdet_0",s_cpdet(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("Pdet_1",s_cpdet(tt,:,2),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("totN_0",s_ctotn(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("totN_1",s_ctotn(tt,:,2),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("totP_0",s_ctotp(tt,:,1),nblocks,block_list,tt+1,outfile)
!   call hdfwq_writeoutput("totP_1",s_ctotp(tt,:,2),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("cC",    s_c(tt,:,1,1),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("cA",    s_c(tt,:,1,2),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("cDIN_0",s_c(tt,:,1,3),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("cDIN_1",s_c(tt,:,2,3),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("cDIP_0",s_c(tt,:,1,4),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("cDIP_1",s_c(tt,:,2,4),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("Ndet_0",s_c(tt,:,1,5),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("Ndet_1",s_c(tt,:,2,5),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("Pdet_0",s_c(tt,:,1,6),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("Pdet_1",s_c(tt,:,2,6),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("totN_0",s_c(tt,:,1,7),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("totN_1",s_c(tt,:,2,7),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("totP_0",s_c(tt,:,1,8),nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("totP_1",s_c(tt,:,2,8),nblocks,block_list,tt+1,outfile)
   ! ---- KK 8.8.24 save cumulative N2 fixation in [tn]
   call hdfwq_writeoutput("rN2fixFC",s_c(tt,:,1,9)/1.0E9,nblocks,block_list,tt+1,outfile)
   call hdfwq_writeoutput("cN2fixFC",s_c(tt,:,2,9),nblocks,block_list,tt+1,outfile)
   ! --------------------------------------------------
END IF

! Update date and write it to file
call incrementdate(model_date,1)
call hdfwq_writedate(model_date,tt+1,outfile)

! Write thickness to file
call hdfwq_writethickness(nblocks,b_thick(tt,:,1),b_thick(tt,:,2),&
                          block_list,tt+1,outfile)

write(*,'(A,3I5)')"Date after timestep: ",model_date(1),model_date(2),model_date(3)
!write(*,'(A27A10A1F20.3)')"Surface cC max at block    ",block_list(maxloc(s_cc(tt,:,1))),":",maxval(s_cc(tt,:,1))
!write(*,'(A27A10A1F20.3)')"Surface cA max at block    ",block_list(maxloc(s_ca(tt,:,1))),":",maxval(s_ca(tt,:,1))
!write(*,'(A27A10A1F20.3)')"Surface cDIN max at block  ",block_list(maxloc(s_cdin(tt,:,1))),":",maxval(s_cdin(tt,:,1))
!write(*,'(A27A10A1F20.3)')"Surface cDIP max at block  ",block_list(maxloc(s_cdip(tt,:,1))),":",maxval(s_cdip(tt,:,1))
!write(*,'(A27A10A1F20.3)')"Surface cTOTN max at block ",block_list(maxloc(s_ctotn(tt,:,1))),":",maxval(s_ctotn(tt,:,1))
!write(*,'(A27A10A1F20.3)')"Surface cTOTP max at block ",block_list(maxloc(s_ctotp(tt,:,1))),":",maxval(s_ctotp(tt,:,1))
write(*,'(A27A10A1F20.3)')"Surface cC max at block    ",block_list(maxloc(s_c(tt,:,1,1))),":",maxval(s_c(tt,:,1,1))
write(*,'(A27A10A1F20.3)')"Surface cA max at block    ",block_list(maxloc(s_c(tt,:,1,2))),":",maxval(s_c(tt,:,1,2))
write(*,'(A27A10A1F20.3)')"Surface cDIN max at block  ",block_list(maxloc(s_c(tt,:,1,3))),":",maxval(s_c(tt,:,1,3))
write(*,'(A27A10A1F20.3)')"Surface cDIP max at block  ",block_list(maxloc(s_c(tt,:,1,4))),":",maxval(s_c(tt,:,1,4))
write(*,'(A27A10A1F20.3)')"Surface cTOTN max at block ",block_list(maxloc(s_c(tt,:,1,7))),":",maxval(s_c(tt,:,1,7))
write(*,'(A27A10A1F20.3)')"Surface cTOTP max at block ",block_list(maxloc(s_c(tt,:,1,8))),":",maxval(s_c(tt,:,1,8))

!IF (ANY(s_cc(tt,:,2).NE.0.0)) WRITE(*,*) "ALERT: cC not 0 in layer 2."
!IF (ANY(s_ca(tt,:,2).NE.0.0)) WRITE(*,*) "ALERT: cA not 0 in layer 2."
IF (ANY(s_c(tt,:,2,1).NE.0.0)) WRITE(*,*) "ALERT: cC not 0 in layer 2."
IF (ANY(s_c(tt,:,2,2).NE.0.0)) WRITE(*,*) "ALERT: cA not 0 in layer 2."

!!! End time loop
END DO DayLoop ! End time-loop

!------------------------------
! Part 4 END MODEL AND FINALIZE
!------------------------------
WRITE(*,*) "WQFICOS DONE, finalizing and closing files."

!---Close HDF files and interface
  CALL h5fclose_f(infile, hdferr)
  CALL h5fclose_f(outfile, hdferr)
  call h5fclose_f(loadingfile, hdferr)
  CALL h5close_f(hdferr)

!---Write mass fluxes over selected interfaces to file if
!---fluxboundaryblocks.dat is present. 
IF (saveifacefluxes) THEN
  WRITE(*,*) "Saving fluxes over interfaces"
  ALLOCATE(fluxin(nlayers,nvars))
  ALLOCATE(fluxout(nlayers,nvars))

  OPEN(UNIT=201,FILE=TRIM(rundir)//"/masstransfer_to_fluxblocks.txt")
  OPEN(UNIT=202,FILE=TRIM(rundir)//"/masstransfer_from_fluxblocks.txt")

  WRITE(201,*) "columns: time <vars_top vars_deep>; nvars:", nvars, "; fluxblocks:", fluxblocks
  WRITE(202,*) "columns: time <vars_top vars_deep>; nvars:", nvars, "; fluxblocks:", fluxblocks
  fluxtimeloop: DO tt=1,ntims
    fluxin = 0.0
    fluxout = 0.0
    ! Sum fluxes in/out of any user defined blocks for timestep (day)
    fluxblockloop: DO ii=1,SIZE(fluxblocks)
      fluxifaceloop: DO nn=1,nifaces

        lb = inter_list(1,nn) ! Left block id in interface
        rb = inter_list(2,nn) ! Right block id in interface
        lb_i = iface_block_list(nn,1) ! Left block index
        rb_i = iface_block_list(nn,2) ! Right block index

        ! Only sum interfaces that have a flux block on one side of the interface
        IF (lb .EQ. fluxblocks(ii) .AND. (.NOT. rb .EQ. fluxblocks(ii))) THEN
           fluxin = fluxin + ifacefluxes(nn,tt,:,1,:)
           fluxout = fluxout + ifacefluxes(nn,tt,:,2,:)
!           WRITE(*,*) "C [kg]", fluxblocks(ii), "-->", SUM(ifacefluxes(nn,:,:,1,1))/19371.0 ! in ww kg/m3
!           WRITE(*,*) "C [kg]", fluxblocks(ii), "<--", SUM(ifacefluxes(nn,:,:,2,1))/19371.0 ! in ww kg/m3
        ELSEIF (rb .EQ. fluxblocks(ii) .AND. (.NOT. lb .EQ. fluxblocks(ii))) THEN
           fluxin = fluxin + ifacefluxes(nn,tt,:,2,:)
           fluxout = fluxout + ifacefluxes(nn,tt,:,1,:)
        END IF

      END DO fluxifaceloop
    END DO fluxblockloop

    WRITE(201,'(I0*(XG0))') tt, fluxin(1,:), fluxin(2,:)
    WRITE(202,'(I0*(XG0))') tt, fluxout(1,:), fluxout(2,:)
!    WRITE(*,*) "IN", tt, fluxin(1,:), fluxin(2,:)
!    WRITE(*,*) "OUT", tt, fluxout(1,:), fluxout(2,:)
  END DO fluxtimeloop
  CLOSE(201)
  CLOSE(202)
ENDIF

!---Write ASCII output
writeasciifiles: IF (asciioutput) THEN
  write(*,*)"Writing ascii output"

  write(*,*)"Writing CC"
  open(unit=101,file=TRIM(rundir)//"/cc_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,1),kk=1,nblocks)
  end do
  close(101)
!  open(unit=101,file="cc_2_out.txt")
!  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
!  do tt=1,ntims
!    write(101,*)(s_cc(tt,kk,2),kk=1,nblocks)
!  end do
!  close(101)

  write(*,*)"Writing CA"
  open(unit=101,file=TRIM(rundir)//"/ca_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,2),kk=1,nblocks)
  end do
  close(101)
!  open(unit=101,file="ca_2_out.txt")
!  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
!  do tt=1,ntims
!    write(101,*)(s_ca(tt,kk,2),kk=1,nblocks)
!  end do
!  close(101)

  write(*,*)"Writing CDIN"
  open(unit=101,file=TRIM(rundir)//"/cdin_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,3),kk=1,nblocks)
  end do
  close(101)
  open(unit=101,file=TRIM(rundir)//"/cdin_2_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,2,3),kk=1,nblocks)
  end do
  close(101)

  write(*,*)"Writing CDIP"
  open(unit=101,file=TRIM(rundir)//"/cdip_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,4),kk=1,nblocks)
  end do
  close(101)
  open(unit=101,file=TRIM(rundir)//"/cdip_2_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,2,4),kk=1,nblocks)
  end do
  close(101)

  write(*,*)"Writing CNDET"
  open(unit=101,file=TRIM(rundir)//"/cndet_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,5),kk=1,nblocks)
  end do
  close(101)
  open(unit=101,file=TRIM(rundir)//"/cndet_2_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,2,5),kk=1,nblocks)
  end do
  close(101)

  write(*,*)"Writing CPDET"
  open(unit=101,file=TRIM(rundir)//"/cpdet_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,6),kk=1,nblocks)
  end do
  close(101)
  open(unit=101,file=TRIM(rundir)//"/cpdet_2_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,2,6),kk=1,nblocks)
  end do
  close(101)

  write(*,*)"Writing TOTP"
  open(unit=101,file=TRIM(rundir)//"/totp_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,7),kk=1,nblocks)
  end do
  close(101)
  open(unit=101,file=TRIM(rundir)//"/totp_2_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,2,7),kk=1,nblocks)
  end do
  close(101)

  write(*,*)"Writing TOTN"
  open(unit=101,file=TRIM(rundir)//"/totn_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,8),kk=1,nblocks)
  end do
  close(101)
  open(unit=101,file=TRIM(rundir)//"/totn_2_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,2,8),kk=1,nblocks)
  end do
  close(101)
  
  ! ----------- KK 16.04.24 --------------------
  ! -- writing cumulative N2fixFC to ASCII-file 
  write(*,*)"Writing rN2fixFC"
  open(unit=101,file=TRIM(rundir)//"/rN2fixFC_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,1,9),kk=1,nblocks)
  end do
  close(101)
  
  write(*,*)"Writing cN2fixFC"
  open(unit=101,file=TRIM(rundir)//"/cN2fixFCout.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
  do tt=1,ntims
    write(101,*)(s_c(tt,kk,2,9),kk=1,nblocks)
  end do
  close(101)
  
  ! --------------------------------------------

  write(*,*)"Writing chla"
  open(unit=101,file=TRIM(rundir)//"/chla_1_out.txt")
  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)

  model_date=start_date ! Re-initialise model_date to start_date
  do tt=1,ntims
    ! chla-a conversion (JR FIX 2017-10-27)
    !    months 1-5:  Chl-a = 5.67857*AN/30.0
    !    months 6-12: Chl-a = 5.67857*(AN+FCN)/20.0
    if (model_date(2).ge.1.and.model_date(2).le.5) then ! Months 1-5
      write(101,*) (((s_c(tt,kk,1,2)/30.0)*5.67857), kk=1, nblocks)
    else ! Months 6-12
      write(101,*) (((s_c(tt,kk,1,1)/20.0+s_c(tt,kk,1,2)/20.0)*5.67857), kk=1, nblocks)
    endif
    call incrementdate(model_date,1) ! Increment model date by 1 day
  end do
  close(101)

!  write(*,*)"Writing height_0"
!  open(unit=101,file="h0_out.txt")
!  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
!  do tt=1,ntims
!    write(101,*)(b_thick(tt,kk,1),kk=1,nblocks)
!  end do
!  close(101)
!
!  write(*,*)"Writing height_1"
!  open(unit=101,file="h1_out.txt")
!  write(101,FMT='(*(A11))')(block_list(kk),kk=1,nblocks)
!  do tt=1,ntims
!    write(101,*)(b_thick(tt,kk,2),kk=1,nblocks)
!  end do
!  close(101)

END IF writeasciifiles

!! End model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                     !!
!! END OF MAIN PROGRAM !!
!!                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM wqficos
