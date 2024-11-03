!================================================
! Helper subroutines for the water quality models
!================================================
MODULE wqficos_helper

CONTAINS

SUBROUTINE INITVALUES(inifilename,paramsfile)
! Subroutine to initialize parameter and variable values
! First give all parameters some default value and then check
! the ini-file to see whether there are some user-defined values
!
! Kai Rasmus 2015
! Janne Ropponen 2017 (FICOS)

USE wqficos_names
IMPLICIT NONE

!---Input arguments
 character(len=128) :: inifilename
 character(len=128) :: paramsfile

!---Local variables
 logical            :: error
 CHARACTER(len=40)  :: value
 integer            :: kk ! Loop counter
 character(len=12)  :: userdefparam
 character(len=3)   :: skk

!==============================================================
! Default FICOS Lignell model constants
! Names are introduced in module wqficos_names.f90
!==============================================================
amaxAN      = 4.2857143 ! = 0.0025*1000.0*24.0/14.0  ! [L ugN-1 d-1] DIN uptake affinity
amaxFCN     = 0.5142857 ! = 0.0003*1000.0*24.0/14.0  ! [L ugN-1 d-1] DIN uptake affinity
Arefuge     = 0.00966   ! = 14.0*13.8/(20.0*1000.0)  ! [ugN/L] minimum concentration of algae A
FCrefuge    = 0.4844    ! = 14.0*692.0/(20.0*1000.0) ! [ugN/L] minimum concentration of algae FC
Klight      = 10.0      ! Ripa algae fix 2019-12-13
KlightA     = 15.0      ! [MJ m-2 d-1] half saturation coefficent of radiation for algae A
KlightFC    = 20.0      ! [MJ m-2 d-1] half saturation coefficent of radiation for algae FC
LightAttIce = 0.0       ! [-] radiation attenuation by ice
LightN2fix  = 15.0      ! [MJ m-2 d-1] Threshold light availability for N2 fixation (FC)
LightThres  = 10.0      ! Ripa algae fix 2019-12-13
maxN2fix    = 0.192     ! = 0.008*24.0  ! [d-1] maximum N2 fixation rate by FC
nutLimitStart = 0.5     ! [-] minimum Michaelis Menten nutrient limitation term
PNratioA    = 0.1383929 ! = 31.0/(14.0*16.0)    ! [-] P:N mass ratio (algae A)
PNratioFC   = 0.1383829 ! = 31.0/(14.0*16.0)    ! [-] P:N mass ratio (algae FC)
RAmax       = 0.144     ! = 0.006*24.0  ! [d-1] maximal loss rate of algae A
reminNmax   = 0.0192    ! = 0.0008*24.0 ! [d-1] maximal detritus N mineralisation rate (beta_0)
reminPmax   = 0.0432    ! = 0.0018*24.0 ! [d-1] maximal detritus P minerailsation rate (gamma_0)
resuspStartRate = 0.10  ! [m/s] threshold bottom flow velocity for resuspension
RFCmax      = 0.096     ! = 0.004*24.0  ! [d-1] maximal loss rate of algae FC
sedrate     = 0.96      ! = 0.04*24.0   ! [m d-1] settling rate
springEndSedCoef = 5.0  ! [-] coefficient of sedimentation rate increase during end of spring bloom (after DIN exhaustion)
SumAFCmax   = 294.0     ! = 14.0*21.0       ! [ugN L-1] maximum total concentration of algae by weight
TcoefA      = 1.001     ! [-] Coefficient for temperature limiting factor for growth of algae A
TcoefFC     = 1.14      ! [-] Coefficient for temperature limiting factor for growth of algae FC
TcoefR      = 1.05      ! [-] Coefficient for temperature limiting factor for losses
TcoefNremin = 1.31      ! [-] Coefficient for temperature limiting factor for detritus N mineralisation
TcoefPremin = 1.6       ! [-] Coefficient for temperature limiting factor for detritus P mineralisation
ToptA       = 15        ! [degC] optimal temperature for the growth of algae A
ToptFC      = 25        ! [degC] optimal temperature for the growth of algae FC
ToptR       = 25        ! [degC] optimal temperature for losses
ToptNremin  = 18        ! [degC] optimal temperature for detritus N mineralisation
ToptPremin  = 18        ! [degC] optimal temperature for detritus P mineralisation
TSpringEnd  = 6.0       ! [degC] maximum temperature for A spring bloom
TN2fixStart = 15.0      ! [degC] threshold temperature for N2 fixation (FC)
umaxA       = 0.696     ! = 0.029*24.0  ! [d-1] maximal growth rate of algae A
umaxFC      = 0.408     ! = 0.017*24.0  ! [d-1] maximal growth rate of algae FC

! HK & JR blue algae fix, 2019-09-03
!N2fixStartDay = 123     ! N2 fixing possible within 2 months of midsummer due to available light
!N2fixEndDay   = 243     ! N2 fixing possible within 2 months of midsummer due to available light

!---Read (changed) model constants from paramsfile
call readsettingsfile(paramsfile)
  
!---Check ini-file for parameter values
call setIniFilename(inifilename)
call loadOptions

!---Parse user defined parameters
do kk=1,100
  if(kk.lt.10)then
    write(skk,'(I1)')kk
    userdefparam='userdef('//skk(1:1)//')'
  end if
  if(kk.lt.100.and.kk.ge.10)then
    write(skk,'(I2)')kk
    userdefparam='userdef('//skk(1:2)//')'
  end if
  if(kk.eq.100)then
    write(skk,'(I3)')kk
    userdefparam='userdef('//skk//')'
  end if
  call getValue('parameters',userdefparam, value, error)
  if(error.eqv..FALSE.)then
    read(value,'(f10.5)')p_userdef(kk)
  end if
end do

END SUBROUTINE




FUNCTION getFluxBorders(fluxfile) RESULT(fluxborders)
  USE prec
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: fluxfile
  INTEGER(KIND=dp), DIMENSION(:), ALLOCATABLE :: fluxborders ! Variable size Fortran 2003 return array

  !---Local variables
  INTEGER(KIND=dp), DIMENSION(:), ALLOCATABLE :: fluxtemp
  INTEGER(KIND=dp) :: temp, nborders, io

  OPEN(UNIT=150, FILE=fluxfile, STATUS='OLD', ACTION='READ')
  nborders = 0
  
  DO
    READ(150,*,IOSTAT=io) temp
    IF (io .NE. 0) EXIT
    IF (temp .GT. 0) THEN
       nborders = nborders + 1
       ALLOCATE(fluxtemp(1:nborders))
       IF ( ALLOCATED(fluxborders) ) fluxtemp(1:nborders-1) = fluxborders
       CALL MOVE_ALLOC(fluxtemp, fluxborders)
       fluxborders(nborders) = temp
    END IF
  END DO

  CLOSE(150)

  RETURN

END FUNCTION getFluxBorders




SUBROUTINE readsettingsfile(paramsfile)

USE prec
USE wqficos_names

IMPLICIT NONE

CHARACTER(len=100) :: buffer, paramname
character(len=128) :: paramsfile
real(kind=dp)      :: value
integer            :: ios = 0
logical            :: file_exists ! Added 2017-09-04 JR

!---Check if settings file exists (Added 2018-09-04 JR)
INQUIRE(FILE=TRIM(paramsfile),EXIST=file_exists)
IF (.NOT. file_exists) THEN
   WRITE(*,*) "Can't find settings file "//TRIM(paramsfile)
   STOP
ENDIF

!open(19, file='ficossettings.ini') ! Commented out 2017-09-04 JR
open(19, file=TRIM(paramsfile),iostat=ios,action='READ') ! Open file 'paramsfile'

do while (ios == 0)
  read(19,'(A)',iostat=ios) buffer
  if ((ios == 0) .and. (LEN(TRIM(buffer)).GT.3))then
    write(*,*) buffer
    read(buffer,*) paramname, value

    SELECT CASE (paramname)
      CASE("amaxAN")
        amaxAN  = value
      CASE("amaxFCN")
        amaxFCN = value
      CASE("Arefuge")
        Arefuge = value
      CASE("FCrefuge")
         FCrefuge = value
      CASE("Klight")
         Klight = value
      CASE("KlightA")
         KlightA = value
      CASE("KlightFC")
         KlightFC = value
      CASE("LightAttIce")
         LightAttIce = value
      CASE("LightN2fix")
         LightN2fix = value
      CASE("LightThres")
         LightThres = value
      CASE("maxN2fix")
         maxN2fix = value
      CASE("nutLimitStart")
         nutLimitStart = value
      CASE("PNratioA")
         PNratioA = value
      CASE("PNratioFC")
         PNratioFC = value
      CASE("RAmax")
         RAmax = value
      CASE("reminNmax")
         reminNmax = value
      CASE("reminPmax")
         reminPmax = value
      CASE("resuspStartRate")
         resuspStartRate = value
      CASE("RFCmax")
         RFCmax = value
      CASE("sedrate")
         sedrate= value 
      CASE("springEndSedCoef")
         springEndSedCoef = value
      CASE("SumAFCmax")
         SumAFCmax = value
      CASE("TcoefA")
         TcoefA= value
      CASE("TcoefFC")
         TcoefFC= value
      CASE("TcoefR")
         TcoefR= value
      CASE("TcoefNremin")
         TcoefNremin= value
      CASE("TcoefPremin")
         TcoefPremin= value
      CASE("ToptA")
         ToptA= value
      CASE("ToptFC")
         ToptFC= value
      CASE("ToptR")
         ToptR= value
      CASE("ToptNremin")
         ToptNremin= value
      CASE("ToptPremin")
         ToptPremin= value
      CASE("TSpringEnd")
         TSpringEnd= value
      CASE("TN2fixStart")
         TN2fixStart= value
      CASE("umaxA")
         umaxA= value
      CASE("umaxFC")
         umaxFC= value
!      CASE("N2fixStartDay")
!         N2fixStartDay = NINT(value)
!      CASE("N2fixEndDay")
!         N2fixEndDay = NINT(value)
      CASE DEFAULT
         WRITE(*,*) "Unknown parameter name ", paramname
         STOP
    END SELECT
  end if
end do

close(19)

END SUBROUTINE readsettingsfile


SUBROUTINE listparams()
! Subroutine to list parameter values
! Kai Rasmus 2015
USE wqficos_names
IMPLICIT NONE

write(*,'(A,F10.5)')"amaxAN: ", amaxAN
write(*,'(A,F10.5)')"amaxFCN: ", amaxFCN 
write(*,'(A,F10.5)')"Arefuge: ", Arefuge
write(*,'(A,F10.5)')"FCrefuge: ", FCrefuge 
write(*,'(A,F10.5)')"Klight: ", Klight
write(*,'(A,F10.5)')"KlightA: ", KlightA 
write(*,'(A,F10.5)')"KlightFC: ", KlightFC 
write(*,'(A,F10.5)')"LightAttIce: ", LightAttIce 
write(*,'(A,F10.5)')"LightThres: ", LightThres
write(*,'(A,F10.5)')"maxN2fix: ", maxN2fix 
write(*,'(A,F10.5)')"nutLimitStart: ", nutLimitStart 
write(*,'(A,F10.5)')"PNratioA: ", PNratioA 
write(*,'(A,F10.5)')"PNratioFC: ", PNratioFC 
write(*,'(A,F10.5)')"RAmax: ", RAmax 
write(*,'(A,F10.5)')"reminNmax: ", reminNmax 
write(*,'(A,F10.5)')"reminPmax: ", reminPmax 
write(*,'(A,F10.5)')"resuspStartRate: ", resuspStartRate 
write(*,'(A,F10.5)')"RFCmax: ", RFCmax 
write(*,'(A,F10.5)')"sedrate: ", sedrate 
write(*,'(A,F10.5)')"springEndSedCoef: ", springEndSedCoef 
write(*,'(A,F10.5)')"SumAFCmax: ", SumAFCmax 
write(*,'(A,F10.5)')"TcoefA: ", TcoefA
write(*,'(A,F10.5)')"TcoefFC: ", TcoefFC
write(*,'(A,F10.5)')"TcoefR: ", TcoefR
write(*,'(A,F10.5)')"TcoefNremin: ", TcoefNremin
write(*,'(A,F10.5)')"TcoefPremin: ", TcoefPremin
write(*,'(A,F10.5)')"ToptA: ", ToptA
write(*,'(A,F10.5)')"ToptFC: ", ToptFC
write(*,'(A,F10.5)')"ToptR: ", ToptR
write(*,'(A,F10.5)')"ToptNremin: ", ToptNremin
write(*,'(A,F10.5)')"ToptPremin: ", ToptPremin
write(*,'(A,F10.5)')"TSpringEnd: ", TSpringEnd
write(*,'(A,F10.5)')"TN2fixStart: ", TN2fixStart
write(*,'(A,F10.5)')"umaxA: ", umaxA
write(*,'(A,F10.5)')"umaxFC: ", umaxFC
!write(*,'(A,I3)')"N2fixStartDay: ", N2fixStartDay
!write(*,'(A,I3)')"N2fixEndDay: ", N2fixEndDay

END SUBROUTINE listparams


!---Function to get index of block
!---Kai Rasmus 2015
SUBROUTINE getblockindex(nblocks,block_list,block,idx)
  IMPLICIT NONE

  integer, intent(in)  :: nblocks
  character(len=10), dimension(nblocks), intent(in) :: block_list
  integer, intent(in)  :: block
  integer, intent(out) :: idx ! Output index

  integer              :: ii
  character(len=10)    :: sblock

  idx = 0
  write(sblock,'(i10)') block

  do ii=1,nblocks
    if (block_list(ii).eq.sblock) then
       idx = ii
    end if
  end do

  if( idx.eq.0) write(*,*) "BC ERROR: ", block

END SUBROUTINE getblockindex


!-----------------------------------------------------------------
! Subroutine to get indii of interfaces for specific blocks
! Kai Rasmus 2015
!-----------------------------------------------------------------
subroutine getifaceindii(nblocks,nifaces,ifaces_right,ifaces_left,&
             n_ifaces_right,n_ifaces_left,inter_list,block_list,max_ifaces)

  implicit none

! Input list
  integer :: nblocks ! Number of blocks
  integer :: nifaces ! Number of interfaces
  integer :: max_ifaces ! Maxmimum number of interfaces for a block
  integer,dimension(nblocks,max_ifaces) :: ifaces_right 
  integer,dimension(nblocks,max_ifaces) :: ifaces_left
  integer,dimension(nblocks) :: n_ifaces_right
  integer,dimension(nblocks) :: n_ifaces_left
  integer,dimension(2,nifaces) :: inter_list ! Interface list (left,right)
!  integer,dimension(nblocks) :: block_list
  character(len=10),dimension(nblocks) :: block_list

! Local variables
  integer :: ii,jj,kk,ll ! Loop counters
  integer :: iblock  

! Initialize
  n_ifaces_right=0
  n_ifaces_left=0
  ifaces_right=0
  ifaces_right=0

  do ii=1,nblocks
    read(block_list(ii),*)iblock
    kk=1
    ll=1
    do jj=1,nifaces
        ! Left interface
        if(iblock.eq.inter_list(1,jj)) then
          n_ifaces_left(ii)=n_ifaces_left(ii)+1
          ifaces_left(ii,kk)=jj                  ! Save interface index
          kk=kk+1
        end if 
        ! Right interface
        if(iblock.eq.inter_list(2,jj)) then
          n_ifaces_right(ii)=n_ifaces_right(ii)+1
          ifaces_right(ii,ll)=jj                 ! Save interface index
          ll=ll+1
        end if

    end do
  end do

end subroutine


!---------------------------------------------------------
! Subroutine to calculate advection of a substance
!
! Called from main program
!
! Values updated only for nonboundary blocks.
!
! Kai Rasmus 2015
! Janne Ropponen 2019
!---------------------------------------------------------
SUBROUTINE calcadv(ntims,nlayers,nblocks,nifaces,nvars,&
                s_c,s_dc,f_avl,f_avr,f_diff,deltat,&
                tim,rb_i,lb_i,nn,b_volume_l,b_volume_r,&
                bbl,bbr,debugifin)

USE prec
USE wqficos_names

IMPLICIT NONE

!---Input arguments
integer, intent(in) :: ntims, nlayers, nblocks, nifaces, nvars    ! Dimensions
real(kind=dp), intent(in)  :: deltat                              ! Timestep
real(kind=dp),dimension(ntims,nblocks,nlayers,nvars), intent(in) :: s_c
real(kind=dp),dimension(nblocks,nlayers,nvars), intent(inout) :: s_dc
real(kind=dp),dimension(nifaces,ntims,nlayers), intent(in) :: f_avl, f_avr, f_diff
integer, intent(in) :: tim                                 ! Timestep index
integer, intent(in) :: rb_i,lb_i                           ! Block indii
real(kind=dp),dimension(nlayers), intent(in) :: b_volume_l, b_volume_r ! Block volume, left/right
integer, intent(in) :: nn                                  ! Interface index
logical, intent(in) :: bbl, bbr ! Flag to indicate whether block is a boundary block
logical, optional, intent(in)   :: debugifin

!---Local variables
integer :: ll,vv               ! Loop counter

real(kind=dp),dimension(2) :: delta_adv ! Concentration change
real(kind=dp)  :: c_rb, c_lb       ! right/left block concentration
real(kind=dp)  :: qq_l, qq_r, qq_d ! Masses transferred in left and right advection, and diffusion
logical :: debugif

!---Mass transfer calculation by var and by layer
varloop: DO vv=1, nvars
layerloop: DO ll=1, nlayers
  delta_adv(1) = 0.0
  delta_adv(2) = 0.0

  c_rb   = s_c(tim,rb_i,ll,vv)
  c_lb   = s_c(tim,lb_i,ll,vv)

  !---Advection
  IF (f_advection.eq.1.and.b_volume_l(ll).gt.0.0.and.b_volume_r(ll).gt.0.0) THEN
    qq_l         =  f_avl(nn,tim,ll) * deltat * c_rb ! Right to left [m3*kg/m3]
    qq_r         = -f_avr(nn,tim,ll) * deltat * c_lb ! Left to right [m3*kg/m3]
    delta_adv(1) = (-qq_r + qq_l) / b_volume_l(ll) ! 1 left block, 2 right block
    delta_adv(2) = (-qq_l + qq_r) / b_volume_r(ll)
  END IF ! End advection

  !---Diffusion
  IF (f_diffusion.EQ.1.AND.b_volume_l(ll).GT.0.0.AND.b_volume_r(ll).GT.0.0) THEN
    qq_d         = f_diff(nn,tim,ll) * deltat * (c_rb - c_lb)
    delta_adv(1) = delta_adv(1) - qq_d / b_volume_l(ll)
    delta_adv(2) = delta_adv(2) + qq_d / b_volume_r(ll)
  END IF ! End diffusion

  !---Update concentration values only for non-boundary blocks.
  !---Note that concentration can get negative here, however total mass is conserved(?)
  IF (.NOT. bbr) THEN ! Right block
    s_dc(rb_i,ll,vv) = s_dc(rb_i,ll,vv) + delta_adv(2)
  END IF

  IF (.NOT. bbl) THEN ! Left block
    s_dc(lb_i,ll,vv) = s_dc(lb_i,ll,vv) + delta_adv(1)
  END IF

END DO layerloop
END DO varloop

!---START debug block
debugif = .FALSE.
IF (PRESENT(debugifin)) debugif = debugifin

IF (debugif) THEN
   WRITE(*,*) "DEBUG: calcadv", tim, lb_i, rb_i, nn, deltat, &
              s_c(tim,lb_i,1,1), s_c(tim,lb_i,2,1), &
              s_c(tim,rb_i,1,1), s_c(tim,rb_i,2,1), &
              s_dc(lb_i,1,1), s_dc(lb_i,2,1), &
              s_dc(rb_i,1,1), s_dc(rb_i,2,1), &
              f_avl(nn,tim,1), f_avl(nn,tim,2), &
              f_avr(nn,tim,1), f_avr(nn,tim,2), &
              b_volume_l(1)/1.0E9, b_volume_l(2)/1.0E9, &
              b_volume_r(1)/1.0E9, b_volume_r(2)/1.0E9
END IF
!---END debug block
  
END SUBROUTINE calcadv


!---------------------------------------------------------
! SUBROUTINE calcvoldif
!
! Calculates the net transferred volume to left and
! right block. All layers are handled independently.
! Used to determine possible "virtual" volume imbalance
! between top and bottom layers during intra-day timestep 
! after calculating mass transfers over all interfaces.
! 
! Janne Ropponen, 2019
! Finnish Environment Institute SYKE
!---------------------------------------------------------
SUBROUTINE calcvoldif(ntims,nlayers,nblocks,nifaces,&
                f_avl,f_avr,voldif,deltat,&
                tim,rb_i,lb_i,nn,bbl,bbr)

USE prec
USE wqficos_names

IMPLICIT NONE

!---Input arguments
integer, intent(in) :: ntims, nlayers, nblocks, nifaces    ! Dimensions
real(kind=dp), intent(in)  :: deltat                              ! Timestep
real(kind=dp),dimension(nifaces,ntims,nlayers), intent(in) :: f_avl, f_avr
real(kind=dp),dimension(nblocks,nlayers), intent(inout) :: voldif
integer, intent(in) :: tim                                 ! Timestep index
integer, intent(in) :: rb_i,lb_i                           ! Block indii
integer, intent(in) :: nn                                  ! Interface index
logical, intent(in) :: bbl, bbr ! Flag to indicate whether block is a boundary block

!---Mass transfer calculation by layer
IF (.NOT. bbl) voldif(lb_i,:) = voldif(lb_i,:) + (f_avl(nn,tim,:) + f_avr(nn,tim,:)) * deltat ! [m3]
IF (.NOT. bbr) voldif(rb_i,:) = voldif(rb_i,:) - (f_avl(nn,tim,:) + f_avr(nn,tim,:)) * deltat ! [m3]

END SUBROUTINE calcvoldif


!---------------------------------------------------------
! SUBROUTINE virtualverticalmixing
!
! Mixes parts of top and bottom layers in case mass transfer 
! with neighbouring blocks have caused an imbalance between 
! top and bottom layers. E.g. bottom layer has received more
! net mass and top layer has lost net mass, thus making it
! look like mass has transferred from top to bottom layer
! in horizontal transfer.
!
! Janne Ropponen, 2019
! Finnish Environment Institute SYKE
!---------------------------------------------------------
SUBROUTINE virtualverticalmixing(nt,nblocks,nlayers,nvars,t,vari0,vari1,c,vol,dv)
USE prec

IMPLICIT NONE

!---Input arguments
integer, intent(in) :: nt, nblocks, nlayers, nvars ! Dimensions
integer, intent(in) :: t, vari0, vari1 ! timestep, variable index range to mix
real(kind=dp),dimension(nt,nblocks,nlayers,nvars), intent(inout) :: c ! concentration to mix
real(kind=dp),dimension(nt,nblocks,nlayers), intent(in) :: vol ! block volumes
real(kind=dp),dimension(nblocks,nlayers), intent(in) :: dv ! volume differences

!---Local variables
real(kind=dp),dimension(nblocks,nlayers,nvars) :: c2 ! temporary var
integer :: vv

!---Virtual vertical mixing [vvmix 1b: after adding mass transports]
  c2 = c(t,:,:,:) ! Save values
  DO vv=vari0,vari1
    WHERE (vol(t,:,2).GT.0.0 .AND. (ABS(dv(:,1)).GT.0.0 .OR. ABS(dv(:,2)).GT.0.0) )
     
       c(t,:,1,vv) = c2(:,1,vv) - c2(:,1,vv)*ABS(dv(:,1))/vol(t,:,1) &
                         + ( ABS(dv(:,1))*c2(:,1,vv) + ABS(dv(:,2))*c2(:,2,vv) ) &
                         * ( ABS(dv(:,1))/(ABS(dv(:,1))+ABS(dv(:,2))) ) &
                         / vol(t,:,1)
     
       c(t,:,2,vv) = c2(:,2,vv) - c2(:,2,vv)*ABS(dv(:,2))/vol(t,:,2) &
                         + ( ABS(dv(:,1))*c2(:,1,vv) + ABS(dv(:,2))*c2(:,2,vv) ) &
                         * ( ABS(dv(:,2))/(ABS(dv(:,1))+ABS(dv(:,2))) ) &
                         / vol(t,:,2)
    END WHERE
  END DO

END SUBROUTINE virtualverticalmixing



!!
subroutine calcvertadv(ntims,nlayers,nblocks,nvars,vari0,vari1,&
                s_cc,deltat,tim,&
                b_volume,b_ad,b_awdw,b_awup)
! Subroutine to calculate vertical advection and diffusion of a substance
! Positive up
!
!  Kai Rasmus 2015
use prec
use wqficos_names

implicit none

! Input arguments
 integer :: ntims,nlayers,nblocks,nvars,vari0,vari1       ! Dimensions
 real(kind=dp) :: deltat                               ! Timestep
 real(kind=dp),dimension(ntims,nblocks,nlayers,nvars) :: s_cc
 real(kind=dp),dimension(ntims,nblocks) :: b_ad,b_awdw,b_awup
 integer :: tim                                 ! Timestep index
 real(kind=dp),dimension(nblocks,nlayers) :: b_volume  ! Block volumes
 
 integer :: ii,vv                        ! Loop counters
 real(kind=dp) :: vol_top,vol_bot
 real(kind=dp) :: qq_top,qq_bot,qq_d
 real(kind=dp) :: c_top,c_bot

 do vv=vari0,vari1
 do ii=1,nblocks
  vol_top=b_volume(ii,1)
  vol_bot=b_volume(ii,2)

  c_top=s_cc(tim,ii,1,vv)
  c_bot=s_cc(tim,ii,2,vv)


  if(f_vadvection.eq.1.and.vol_top.gt.0.0.and.vol_bot.gt.0.0)then
   qq_bot=-b_awdw(tim,ii)*deltat*c_top
   qq_top=b_awup(tim,ii)*deltat*c_bot
   s_cc(tim,ii,1,vv)=s_cc(tim,ii,1,vv)-qq_bot/vol_top+qq_top/vol_top
   s_cc(tim,ii,2,vv)=s_cc(tim,ii,2,vv)-qq_top/vol_bot+qq_bot/vol_bot
  end if ! End advection

  if(f_vdiffusion.eq.1.and.vol_top.gt.0.0.and.vol_bot.gt.0.0)then
   qq_d=b_ad(tim,ii)*deltat*(c_top-c_bot)
   s_cc(tim,ii,1,vv)=s_cc(tim,ii,1,vv)-qq_d/vol_top
   s_cc(tim,ii,2,vv)=s_cc(tim,ii,2,vv)+qq_d/vol_bot
  end if ! End diffusion

 end do ! End blocks
 end do ! end variables

end subroutine


!!
subroutine incrementdate(model_date, daystoinc)
! Subroutine to increment a dim(6) date variable by 1 day
! Kai Rasmus, 2016
! Janne Ropponen, 2019 (added day of year)
implicit none
! Input arguments
integer,dimension(6),intent(inout) :: model_date
INTEGER, intent(in) :: daystoinc

integer :: year, month, day
integer :: i, dti
integer :: leap=0 ! 0= common year, 1 = leap year
integer :: daysinmon=0 ! 0=28 (29), 1=30 and 2=31
integer, dimension(12) :: nbdays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

year  = model_date(1)
month = model_date(2)
day   = model_date(3)

! Detect leap year
if(mod(year,4).ne.0)then
        leap=0 ! Common year
else ! Year evenly divisible by 4 so let's check the leap-year cases
 if(mod(year,100).eq.0.and.mod(year,400).ne.0)then
        leap=0 ! Common year
 else
        leap=1 ! leap year
 end if
end if

IF (leap .EQ. 1) nbdays(2) = 29

if(month.eq.1.or.month.eq.3.or.month.eq.5&
        .or.month.eq.7.or.month.eq.8.or.month.eq.10.or.month.eq.12)then
        daysinmon=2
else if(month.eq.4.or.month.eq.6.or.month.eq.9&
        .or.month.eq.11)then
        daysinmon=1
else if(month.eq.2)then
        daysinmon=0
end if

dti = daystoinc
DO WHILE (dti .GT. 0)
   ! Increment one day
   model_date(3) = model_date(3) + 1
   dti = dti - 1

   ! Increment month
   if(model_date(3).eq.32.and.daysinmon.eq.2)then
           model_date(3) = 1
           model_date(2) = model_date(2) + 1
   else if(model_date(3).eq.31.and.daysinmon.eq.1)then
           model_date(3) = 1
           model_date(2) = model_date(2) + 1
   else if(model_date(3).eq.29.and.daysinmon.eq.0.and.leap.eq.0)then
           model_date(3) = 1
           model_date(2) = model_date(2) + 1
   else if(model_date(3).eq.30.and.daysinmon.eq.0.and.leap.eq.1)then
           model_date(3) = 1
           model_date(2) = model_date(2) + 1
   end if

   ! Increment year
   if (model_date(2).eq.13) then
           model_date(2) = 1
           model_date(1) = model_date(1) + 1
   end if
END DO

! Day of year
model_date(6) = 0
DO i = 1,month
   IF (i .EQ. month) THEN
      model_date(6) = model_date(6) + day
   ELSE
      model_date(6) = model_date(6) + nbdays(i)
   END IF
END DO

end subroutine

! https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
! Convert integer to string. Requires Fortran 2003.
function itoa(i) result(res)
  character(:),allocatable :: res
  integer,intent(in) :: i
  character(range(i)+2) :: tmp
  write(tmp,'(i0)') i
  res = trim(tmp)
end function


END MODULE wqficos_helper
