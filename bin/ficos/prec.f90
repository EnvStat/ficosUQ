!--- Fortran real precision definitions
!--- See: http://fortranwiki.org/fortran/show/Real+precision
!
! Janne Ropponen, Finnish Environment Institute, 2019

MODULE prec

IMPLICIT NONE

integer, parameter :: sp = selected_real_kind(6, 37) ! single precision (32-bit)
integer, parameter :: dp = selected_real_kind(15, 307) ! double precision (64-bit)
integer, parameter :: qp = selected_real_kind(33, 4931) ! quad precision (128-bit)

END MODULE prec
