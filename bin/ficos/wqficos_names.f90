!--------------------------------------------------------
! Paramater names and default flags for the water quality model
!--------------------------------------------------------
MODULE wqficos_names

USE prec

IMPLICIT NONE

!---Lignell model constants
real(kind=dp) :: amaxAN
real(kind=dp) :: amaxFCN
real(kind=dp) :: Arefuge
real(kind=dp) :: FCrefuge
real(kind=dp) :: Klight
real(kind=dp) :: KlightA
real(kind=dp) :: KlightFC
real(kind=dp) :: LightAttIce
real(kind=dp) :: LightN2fix
real(kind=dp) :: LightThres
real(kind=dp) :: maxN2fix
real(kind=dp) :: PNratioA
real(kind=dp) :: PNratioFC
real(kind=dp) :: RAmax
real(kind=dp) :: reminNmax
real(kind=dp) :: reminPmax
real(kind=dp) :: RFCmax
real(kind=dp) :: sedrate
real(kind=dp) :: SumAFCmax
real(kind=dp) :: TcoefA
real(kind=dp) :: TcoefFC
real(kind=dp) :: TcoefR
real(kind=dp) :: TcoefNremin
real(kind=dp) :: TcoefPremin
real(kind=dp) :: ToptA
real(kind=dp) :: ToptFC
real(kind=dp) :: ToptR
real(kind=dp) :: ToptNremin
real(kind=dp) :: ToptPremin
real(kind=dp) :: TN2fixStart
real(kind=dp) :: umaxA
real(kind=dp) :: umaxFC
real(kind=dp) :: nutLimitStart
real(kind=dp) :: springEndSedCoef
real(kind=dp) :: TSpringEnd
real(kind=dp) :: resuspStartRate
real(kind=dp) :: mmN, mmP

! HK & JR blue algae fix, 2019-09-03
!integer :: N2fixStartDay
!integer :: N2fixEndDay

!---User-defined parameters
real(kind=dp),dimension(100) :: p_userdef

!---Model flags
integer :: f_advection  = 1 ! Advection (0 off, 1 on, 1 default)
integer :: f_diffusion  = 0 ! Horizontal diffusion (0 off, 1 on, 0 default)
integer :: f_vdiffusion = 0 ! Vertical diffusion (0 off, 1 on, 0 default)
integer :: f_vadvection = 1 ! Vertical advection (0 off, 1 on, 1 default)

END MODULE wqficos_names
