!-------------------------------------
! Subroutine sent to ode solver
!
! Called by ode solver
!
! Kai Rasmus 2015
! Janne Ropponen, Risto Lignell 2016-
!------------------------------------
subroutine fex_ficos(t,neq,y,ydot,debugmode)
 USE prec
 USE wqficos_names

 IMPLICIT NONE

 real(kind=dp) templimit

! Input arguments
 real(kind=dp),intent(in)                   :: t    ! Time [days] (not used at this time)
 integer, intent(in)                        :: neq  ! Number of equations
 real(kind=dp),dimension(neq),intent(inout) :: y    ! Input and output values for functions
 real(kind=dp),dimension(neq),intent(out)   :: ydot ! Changes in time of the functions
 logical, intent(in)                        :: debugmode ! If enabled, write debug output.

! Local variables
 real(kind=dp), allocatable, dimension(:) :: A, L   ! Advection (vertical), Loadings

! Lignell model factors / rates
 real(kind=dp) :: AFCshading,AnutLimit,AradLimit,ArefCorr,FCnutLimit,FCradLimit,AFCradLimit
 real(kind=dp) :: FCrefCorr,MMenFCN,MMenFCP,TcorrA,TcorrFC,TcorrR,TcorrNremin
 real(kind=dp) :: TcorrPremin,ANspringSed,N2fix,Nremin,Premin,RA,RFC,uA,uFC
 real(kind=dp) :: TcorrNremin2,TcorrPremin2,Nremin2,Premin2 ! Factors for layer 2 / JR 2019-12
 
! Sedimentation factor
 real(kind=dp) :: sed

! Allocate memory
 allocate(A(neq)) ! Advection (vertical)
 allocate(L(neq)) ! Loadings

! Initialize coefficients
 A = 0.0
 L = 0.0

! The values in y are:
! 1 cc 
! 2 ca 
! 3 cdin 
! 4 cdip 
! 5 s_cdin(tt,mm,2)
! 6 s_cdip(tt,mm,2)

! 7 cndet 
! 8 cpdet 
! 9 s_cndet(tt,mm,2)
! 10 s_cpdet(tt,mm,2)

! 11 ctotn
! 12 ctotp
! 13 s_ctotn(tt,mm,2)
! 14 s_ctotp(tt,mm,2)

! 15 T temperature 
! 16 b_t(tt,mm,2)

! 17 Irad radiation 
! 18 Ice 
! 19 h ! Thickness of block 
! 20 b_thick(mm,2)

! 21 layer number 
! 22 Downwards velocity [m/d] 
! 23 Upwards velocity [m/d] 
! 
! 24 b_volume(mm,1) 
! 25 b_volume(mm,2)

! Loadings [mg/d]
! 26 din layer 1
! 27 din layer 2
! 28 dip layer 1
! 29 dip layer 2
! 30 totn layer 1
! 31 totn layer 2
! 32 totp layer 1
! 33 totp layer 2

! N2 fixation
! 34 N2 fixation [mg]
! 35 N2 fixation [mg/m3]

! Spare variables (36-37)
! 36 
! 37 Day of year

! 38 Bottom velocity

!==============================================================
! DEBUG / checks
!==============================================================
!
!---Algae limitation checks (debug)
 IF ( y(1) .LE. 0.0 ) THEN
    WRITE(*,*) "WARNING: FC concentration too low, ", y(1)
    !y(1) = FCrefuge
 END IF
 IF ( y(2) .LE. 0.0 ) THEN
    WRITE(*,*) "WARNING: A concentration too low, ", y(2)
    !y(2) = Arefuge
 END IF

 !---Positivity check for input y (JR)
 WHERE (y(1:14).LT.0.0) y(1:14)=0.0 ! Concentrations of state variables must be positive

!==============================================================
! Vertical advection (array E)
! y(22) downwards flux [m3/d], always negative
! y(23) upwards flux [m3/d], always positive
! y(24) volume, layer 1
! y(25) volume, layer 2
!==============================================================

 A = 0.0 ! Initialise vertical advection array A(:)
 
 IF ( y(24) .GT. 0.0 .AND. y(25) .GT. 0.0 ) THEN ! Two layers
    !A(1)  = ( y(23)*y(??) + y(22)*y(1)  ) / y(24) ! FCN layer 1 ! Only 1 layer at this time
    !A(2)  = ( y(23)*y(??) + y(22)*y(2)  ) / y(24) ! AN layer 1 ! Only 1 layer at this time
    A(3)  = ( y(23)*y(5)  + y(22)*y(3)  ) / y(24) ! DIN layer 1
    A(4)  = ( y(23)*y(6)  + y(22)*y(4)  ) / y(24) ! DIP layer 1
    A(7)  = ( y(23)*y(9)  + y(22)*y(7)  ) / y(24) ! DetN layer 1
    A(8)  = ( y(23)*y(10) + y(22)*y(8)  ) / y(24) ! DetP layer 1
    A(11) = ( y(23)*y(13) + y(22)*y(11) ) / y(24) ! totN layer 1
    A(12) = ( y(23)*y(14) + y(22)*y(12) ) / y(24) ! totP layer 1

    A(5)  = ( -y(23)*y(5)  - y(22)*y(3)  ) / y(25) ! DIN layer 2
    A(6)  = ( -y(23)*y(6)  - y(22)*y(4)  ) / y(25) ! DIP layer 2
    A(9)  = ( -y(23)*y(9)  - y(22)*y(7)  ) / y(25) ! DetN layer 2
    A(10) = ( -y(23)*y(10) - y(22)*y(8)  ) / y(25) ! DetP layer 2
    A(13) = ( -y(23)*y(13) - y(22)*y(11) ) / y(25) ! totN layer 2
    A(14) = ( -y(23)*y(14) - y(22)*y(12) ) / y(25) ! totP layer 2
    !A(??) = ( -y(23)*y(??) - y(22)*y(1)  ) / y(25) ! FCN layer 2 ! Only 1 layer at this time
    !A(??) = ( -y(23)*y(??) - y(22)*y(2)  ) / y(25) ! AN layer 2 ! Only 1 layer at this time
 END IF

!==============================================================
! Loadings (array F)
! input should be in same units as the equations
!==============================================================
!
 L = 0.0 ! Initialise loadings array F(:)
 
 !---Point loads:  [mg] --> [umolN/m3] == [nmolN/L]
 IF ( y(24) .GT. 0.0 ) THEN ! Top layer
    L(3)  = y(26) / y(24) ! DIN layer 1
    L(4)  = y(28) / y(24) ! DIP layer 1
    L(11) = y(30) / y(24) ! totN layer 1
    L(12) = y(32) / y(24) ! totP layer 1
 END IF 
 
 IF ( y(25) .GT. 0.0 ) THEN ! Bottom layer
    L(5)  = y(27) / y(25) ! DIN layer 2
    L(6)  = y(29) / y(25) ! DIP layer 2
    L(13) = y(31) / y(25) ! totN layer 2
    L(14) = y(33) / y(25) ! totP layer 2
 ELSE ! If there is no bottom layer then add bottom layer loadings to top layer
    L(3)  = L(3)  + y(27) / y(24) ! DIN layer1 + layer2
    L(4)  = L(4)  + y(29) / y(24) ! DIP layer1 + layer2
    L(11) = L(11) + y(31) / y(24) ! totN layer1 + layer2
    L(12) = L(12) + y(33) / y(24) ! totP layer1 + layer2
 END IF

!==============================================================
! Define equations
!     dc/dt = ydot = A + E + F
!
!     A = change in time of the concentration
!     E = vertical advection (vertical transport of nutrients)
!     F = Loading (new external nutrients)
!==============================================================
!
 ydot = 0.0 ! Initialise time derivatives

 !---Sedimentation flag. 1 = sedimentation, 0 = no sedimentation.
 !---No sedimentation if bottom velocity is larger than 0.1 m/s.
 sed = 1.0
 IF ( y(38) .GT. resuspStartRate ) sed = 0.0

 !---LAYER 1
 !
 !---Limiting factors (layer 1)
 ! Ripa algae fix 2019-12-13:
 AFCRadLimit = ( y(17) - LightThres - LightAttIce*y(17) ) / &
               ( y(17) - LightThres - LightAttIce*y(17) + Klight )
 AFCRadLimit = MERGE(AFCRadLimit, 0.0_dp, AFCRadLimit .GE. 0.0_dp)

 AFCshading  = 1.0 - ( y(1) + y(2) ) / SumAFCmax
 AnutLimit   = MIN( (y(3) / ( y(3) + umaxA/amaxAN)), (y(4) / (y(4) + PNratioA*umaxA/amaxAN)) )
 AradLimit   = y(17)*(1 - LightAttIce) / ((y(17)*(1 - LightAttIce)) + KlightA)
 IF (y(2) .GT. Arefuge) THEN
    ArefCorr = (y(2)-Arefuge)/y(2)
 ELSE
    ArefCorr = 0.0
 ENDIF
 FCnutLimit  = MIN( (y(3) / (y(3)+umaxFC/amaxFCN)), (y(4) / (y(4) + PNratioA*umaxFC/amaxFCN)) )
 FCradLimit  = y(17)*(1-LightAttIce)/((y(17)*(1-LightAttIce))+KlightFC)
 IF ( y(1) .GT. FCrefuge) THEN
    FCrefCorr = (y(1)-FCrefuge)/y(1)
 ELSE
    FCrefCorr = 0.0
 ENDIF
 MMenFCN     = y(3) / ( y(3) + ( umaxFC/amaxFCN ) )
 MMenFCP     = y(4) / ( y(4) + ( PNratioFC*umaxFC/amaxFCN ) )

 !---Temperature correction (layer 1)
 TcorrA      = templimit(ToptA,TcoefA,y(15))
 TcorrFC     = templimit(ToptFC,TcoefFC,y(15))
 TcorrR      = templimit(ToptR,TcoefR,y(15))
 TcorrNremin = templimit(ToptNremin,TcoefNremin,y(15))
 TcorrPremin = templimit(ToptPremin,TcoefPremin,y(15))
 
 !---Other terms (layer 1)
 !---ANspringSed = Straight loss term for algae A. Should drop to layer below(?) in multi-layer model
 ANspringSed = MERGE( springEndSedCoef*(sedrate/y(19)), y(2)*0.0, &
                      AnutLimit .LT. nutLimitStart .AND. y(15) .LT. TSpringEnd )
! N2fix       = MERGE( MMenFCP*TcorrFC*maxN2fix, y(1)*0.0, &
!                      MMenFCN .LT. nutLimitStart .AND. y(15) .GT. TN2fixStart )
! Extra N2 fixing for blue-green algae (HK & JR fix 2019-09-04)
! N2fix       = MERGE( MMenFCP*TcorrFC*maxN2fix, y(1)*0.0, &
!                            MMenFCN .LT. nutLimitStart &
!                      .AND. y(15)   .GT. TN2fixStart &
!                      .AND. y(37)   .GT. N2fixStartDay &
!                      .AND. y(37)   .LT. N2fixEndDay &
!                    )

 !---Cyanobacteria N2 fixing: Ripa N2fix fix 5.12.2019/JR
 N2fix       = MERGE( MMenFCP*TcorrFC*maxN2fix, y(1)*0.0, &
                          MMenFCN .LT. nutLimitStart &
                      .AND. y(15) .GT. TN2fixStart &
                      .AND. y(17) .GT. LightN2fix )

 Nremin      = reminNmax*TcorrNremin
 Premin      = reminPmax*TcorrPremin
 
 !---Rates
 RA          = RAmax*TcorrR*ArefCorr
 RFC         = RFCmax*TcorrR*FCrefCorr
! uA          = umaxA*AnutLimit*AradLimit*TcorrA*AFCshading
! uFC         = umaxFC*FCnutLimit*FCradLimit*TcorrFC*AFCshading
! Ripa algae fix 2019-12-13
 uA          = umaxA*AnutLimit*AFCradLimit*TcorrA*AFCshading
 uFC         = umaxFC*FCnutLimit*AFCradLimit*TcorrFC*AFCshading

 !---Equations for rate of concentration change (layer 1)
 ydot(1) = (N2fix + uFC - RFC)*y(1) ! FCN layer 1
 ydot(2) = (uA - RA - ANspringSed)*y(2) ! AN layer 1

 ydot(3) = Nremin*y(7) - uFC*y(1) - uA*y(2) ! DIN layer 1
 ydot(4) = Premin*y(8) - (N2fix + uFC)*PNratioFC*y(1) - uA*PNratioA*y(2) ! DIP layer 1

 ydot(7) = RFC*y(1) + RA*y(2) - Nremin*y(7) ! - (sedrate/zmix)*y(7) ! detN layer 1 
 ydot(8) = RFC*PNratioFC*y(1) + RA*PNratioFC*y(2) - Premin*y(8) ! - (sedrate/zmix)*y(8) ! detP layer 1



 !---Settling and sedimentation of layer 1
 IF ( y(20) .EQ. 0.0 ) THEN ! If lower block depth is zero then burial from surface layer
    ydot(7)  = ydot(7) - (sed*sedrate/y(19))*y(7) ! detN layer 1 (burial)
    ydot(8)  = ydot(8) - (sed*sedrate/y(19))*y(8) ! detP layer 1 (burial)

    ydot(11) = -(sed*sedrate/y(19))*y(7) ! totN layer 1 (burial)
    ydot(12) = -(sed*sedrate/y(19))*y(8) ! totP layer 1 (burial)
 ELSE                                    ! Otherwise settling
    ydot(7)  = ydot(7) - (sedrate/y(19))*y(7) ! detN layer 1 (settling)
    ydot(8)  = ydot(8) - (sedrate/y(19))*y(8) ! detP layer 1 (settling)

    ydot(11) = -(sedrate/y(19))*y(7) ! totN layer 1 (settling)
    ydot(12) = -(sedrate/y(19))*y(8) ! totP layer 1 (settling)
 END IF

 !---LAYER 2
 !
 IF ( y(20) .GT. 0.0 ) THEN ! Bottom layer exists
    !---Limiting factors (layer 2)
    TcorrNremin2 = templimit(ToptNremin,TcoefNremin,y(16))
    TcorrPremin2 = templimit(ToptPremin,TcoefPremin,y(16))
    Nremin2      = reminNmax*TcorrNremin*y(7)
    Premin2      = reminPmax*TcorrPremin*y(8)
 
    !---Equations for rate of concentration change (layer 2)
    ydot(5) = Nremin2*y(9)
    ydot(6) = Premin2*y(10) 

    ydot(9) = -Nremin2*y(9)
    ydot(10) = -Premin2*y(10)

    !---Sedimentation (layer 2)
    ydot(9)  = ydot(9)  - (sed*sedrate/y(20))*y(9)
    ydot(10) = ydot(10) - (sed*sedrate/y(20))*y(10)

    !---Settling from above (to layer 2)
    ydot(9)  = ydot(9)  + (sedrate/y(19))*(y(24)/y(25))*y(7) ! detN layer 2
    ydot(10) = ydot(10) + (sedrate/y(19))*(y(24)/y(25))*y(8) ! detP layer 2
    ydot(13) = (sedrate/y(19))*(y(24)/y(25))*y(7) - &
               (sed*sedrate/y(20))*y(9) ! totN layer 2
    ydot(14) = (sedrate/y(19))*(y(24)/y(25))*y(8) - &
               (sed*sedrate/y(20))*y(10) ! toptP layer 2

 END IF ! Bottom layer
 
 ! ------------------------
 ! Karel Kaurila 2024
 ! N2 fixation
 ydot(34) = N2fix*y(24)*y(1) ! [mg]
 ydot(35) = N2fix*y(1) ! [mg/m3]
 !-------------------------
 
 
 !---Add vertical transport and loads to equations
 ydot(1)  = ydot(1)  + A(1)  + L(1)
 ydot(2)  = ydot(2)  + A(2)  + L(2)
 ydot(3)  = ydot(3)  + A(3)  + L(3)
 ydot(4)  = ydot(4)  + A(4)  + L(4)
 ydot(5)  = ydot(5)  + A(5)  + L(5)
 ydot(6)  = ydot(6)  + A(6)  + L(6)
 ydot(7)  = ydot(7)  + A(7)  + L(7)
 ydot(8)  = ydot(8)  + A(8)  + L(8)
 ydot(9)  = ydot(9)  + A(9)  + L(9)
 ydot(10) = ydot(10) + A(10) + L(10)
 ydot(11) = ydot(11) + A(11) + L(11)
 ydot(12) = ydot(12) + A(12) + L(12)
 ydot(13) = ydot(13) + A(13) + L(13)
 ydot(14) = ydot(14) + A(14) + L(14)

!===================================================
!DEBUG OUTPUT/JR 2019
!===================================================
IF (debugmode) THEN
   OPEN(unit=199,file="./ficos_debug_rk4.txt",action='WRITE',access='APPEND')

   IF (t.EQ.0.0) write(199,*) "---RK4-START---"
   write(199,"(140G12.4)") y(37), y(1:38), A(1:38), L(1:38), uA, RA, uFC, RFC, &
            AnutLimit, AradLimit, ArefCorr, &
            FCnutLimit, FCradLimit, FCrefCorr, &
            ANspringSed, AFCshading, N2fix, Nremin, Premin, MMenFCN, MMenFCP, &
            TcorrA, TcorrFC, TcorrNremin, TcorrPremin, &
            TcorrNremin2, TcorrPremin2, Nremin2, Premin2

   CLOSE(199)
ENDIF

!==============================================================
! Clean up
!==============================================================
!
! Deallocate memory
 deallocate(A)
 deallocate(L)

end subroutine fex_ficos


!!
function templimit(topt,at,t) result(tlimit)
! Function to calculate temperature limitation
! Kai Rasmus 2015
   USE prec

   implicit none

   real(kind=dp), intent(in) :: topt      ! Optimal temperature
   real(kind=dp), intent(in) :: at        ! Limitation parameter
   real(kind=dp), intent(in) :: t         ! Temperature
   real(kind=dp) :: tlimit                ! Temperature limitation output

   tlimit=exp(topt-t+((at/(1.0-at))*topt+t)*log(at+(1.0-at)*t/topt))

end function templimit

