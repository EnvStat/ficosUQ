!************************************************************ 
! 
!  This example shows how to iterate over group members using 
!  H5Literate. 
! 
!  This file is intended for use with HDF5 Library version 1.8 
!  with --enable-fortran2003 
! 
! 
!************************************************************ 
MODULE g_iterate
  USE ISO_C_BINDING
  IMPLICIT NONE

  character(len=256),dimension(40) :: loading_stations
  integer :: numstations
  integer,dimension(40) :: stationlens

CONTAINS

!********************************************************************* 
! 
!  Callback operator function.  Prints the name and type of the object 
!  being examined. 
! 
!  This function is used by the loading subroutine wqhdf_loadings
!
! ********************************************************************

  INTEGER(KIND=C_INT) FUNCTION op_func(loc_id, name, info, operator_data) bind(C)
    USE HDF5
    USE ISO_C_BINDING
    IMPLICIT NONE

    INTEGER(HID_T), VALUE :: loc_id
    CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(1:128) :: name ! must have LEN=1 for bind(C) strings
    TYPE(C_PTR), VALUE :: info
    TYPE(C_PTR), VALUE :: operator_data
    
    INTEGER :: status, i, len

    TYPE(H5O_info_t), TARGET :: infobuf
    CHARACTER(LEN=128) :: name_string

    !
    ! The name of the object is passed to this FUNCTION by
    ! the Library.
    !

    numstations=numstations+1

    ! Get name string from list of characters
    DO i = 1, 127
       loading_stations(numstations)(i:i)=name(i)(1:1)
       name_string(i:i) = name(i)(1:1)
    ENDDO

    CALL H5Oget_info_by_name_f(loc_id, name_string, infobuf, status)

    ! Include the string up to the C NULL CHARACTER
    len = 0
    DO
       IF(name_string(len+1:len+1).EQ.C_NULL_CHAR.OR.len.GE.128) EXIT
       len = len + 1
    ENDDO

    stationlens(numstations)=len

    op_func = 0 ! return successful

  END FUNCTION op_func

END MODULE g_iterate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE hdfio

CONTAINS

subroutine wqhdf_gettimesteps(ntims,deltat,nifaces,s_inter_list,file)
! Subroutine to get timesteps
! Retrieve all timesteps from all interfaces and then find the smallest
! one per timestep
!
!  Kai Rasmus 2015

  use hdf5
  use prec
  implicit none

! Input arguments
  integer :: ntims ! Output: number of timesteps
  integer(HID_T) :: file
  integer :: nifaces
  character(len=21), dimension(nifaces) :: s_inter_list ! List of interfaces in string format
  real(kind=dp), dimension(ntims) :: deltat
  real(kind=dp), dimension(nifaces,ntims) :: tims
  real(kind=dp), dimension(ntims) :: tim

! Handles

  INTEGER(HID_T) :: space1,dset
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group,subgroup

  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of interfaces
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
  integer(hsize_t),dimension(2) :: tmin_dims

  INTEGER :: II,jj ! LOOP COUNTERs

  CALL h5gopen_f(file,"Interfaces", group, hdferr)
  do ii=1,nifaces
        CALL h5gopen_f(group,s_inter_list(ii),subgroup,hdferr)

        CALL h5dopen_f (subgroup,"Tmin", dset, hdferr)
        CALL h5dget_space_f(dset, space1, hdferr)
        CALL h5sget_simple_extent_dims_f(space1,data_dims, blndims, hdferr)
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,tim,tmin_dims,hdferr)

        do jj=1,ntims
           tims(ii,jj)=tim(jj)
        end do

        ! Close subgroup and dataset
        CALL h5gclose_f(subgroup, hdferr)
        CALL h5sclose_f(space1, hdferr)
        CALL h5dclose_f(dset, hdferr)

  end do

  deltat=minval(tims,1)

  ! Close group
  CALL h5gclose_f(group, hdferr)

end subroutine

subroutine wqhdf_getntimesteps(ntims,nifaces,s_inter_list,file)
! Subroutine to get number of timesteps
!  Kai Rasmus 2015

  use hdf5
  implicit none

! Input arguments
  integer :: ntims ! Output: number of timesteps
  integer(HID_T) :: file
  integer :: nifaces
  character(len=21), dimension(nifaces) :: s_inter_list ! List of interfaces in string format
 
! Handles
  INTEGER(HID_T) :: space1,dset
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group,subgroup

  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of interfaces
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  CALL h5gopen_f(file,"Interfaces", group, hdferr)
  CALL h5gopen_f(group,s_inter_list(nifaces),subgroup,hdferr)

  CALL h5dopen_f (subgroup,"Tmin", dset, hdferr)
  CALL h5dget_space_f(dset, space1, hdferr)
  CALL h5sget_simple_extent_dims_f(space1,data_dims, blndims, hdferr)


  ! Close dataset and dataspace
  CALL h5sclose_f(space1, hdferr)
  CALL h5dclose_f(dset, hdferr)
  CALL h5gclose_f(group, hdferr)
  CALL h5gclose_f(subgroup, hdferr)

  ! Define the number of timesteps as the larger of the 2 dimensions.
  ! In some cases the largest value is in the first dimension and in
  ! some cases it is the second dimension. Only Sir H.D.F. Fortran
  ! knows why this is.
  ntims=maxval(blndims)

end subroutine


subroutine wqhdf_getifacelist(nifaces,inter_list,s_inter_list,file)
! Subroutine to extract interface list in integer and string form
! Kai Rasmus 2015

  use hdf5
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: nifaces
  integer, dimension(2,nifaces) :: inter_list ! List of interfaces between blocks
  character(len=21), dimension(nifaces) :: s_inter_list ! List of interfaces in string format

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: hdferr

  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of interfaces
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  integer :: ii ! Loop counter

  CALL h5dopen_f (file,"Inter_list", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

! Read interface list using the default properties.
  CALL h5dread_f(dset,H5T_NATIVE_INTEGER,inter_list,data_dims, hdferr)

  ! Define string type interface list
  DO ii=1,nifaces
        write(s_inter_list(ii),'(I0,A1,I0)')inter_list(1,ii),"_",inter_list(2,ii)
!write(*,*)s_inter_list(ii)
  end do

  ! Close dataset and dataspace
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)


end subroutine



subroutine wqhdf_getnifaces(nifaces,file)
! Subroutine to extract number of interfaces
! Kai Rasmus 2015
  use hdf5
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: nifaces

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: hdferr

  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of interfaces
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  CALL h5dopen_f (file,"Inter_list", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

  ! Close dataset and dataspace
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)

  nifaces=blndims(2) 
end subroutine


!!
!!----------------------------------------------------------------------------
!!
subroutine wqhdf_getnblocks(nblocks,file) 
! Subroutine to extract number of blocks 
! Kai Rasmus 2015
  use hdf5
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: nblocks

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: hdferr

  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of blocks in block list
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  CALL h5dopen_f(file,"Block_list", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space, data_dims, blndims, hdferr)

  ! Close dataset and dataspace
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)

  nblocks=blndims(1)
end subroutine

!!
!!----------------------------------------------------------------------------
!!
subroutine wqhdf_getblocklist(nblocks,block_list,file)
! Subroutine to extract block list
! Kai Rasmus 2015
  use hdf5
  implicit none

! Input arguments
  integer(HID_T) :: file

!  integer, dimension(nblocks) :: block_list
  character(len=10),dimension(nblocks) :: block_list
  integer :: nblocks

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: size,hdferr,filetype2
  INTEGER(HID_T)  :: filetype,memtype
  INTEGER(HSIZE_T), DIMENSION(2) :: maxdims
  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of blocks in block list
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  integer,dimension(nblocks) :: iblock ! Handle to put integer block name in

  integer :: ii ! Loop counter

  CALL h5dopen_f (file,"Block_list", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

  ! Get the datatype and its size.
  CALL H5Dget_type_f(dset, filetype, hdferr)
  CALL H5Tget_size_f(filetype, size, hdferr)

 ! Get dataspace
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims,maxdims, hdferr)

  ! Get filetype type (H5T_INTEGER_F or H5T_STRING_F)
  call h5tget_class_f(filetype,filetype2,hdferr)

  ! Create the memory datatype.
! From https://www.hdfgroup.org/hdf5-quest.html#str1
! HDF5 has a string datatype (H5T_C_S1) that can be modified to create a fixed length
! or variable length array of characters (a string).
!
! Using this datatype without modifying it will create a one-character string. To create
! a string longer than one character, you must get a copy or instance of the datatype,
! and modify it with the H5Tset_size API.

  CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
  CALL H5Tset_size_f(memtype, size, hdferr)

! Read block using the default properties.

  if(filetype2.eq.H5T_INTEGER_F) then
    write(*,*)"Found integer-type blocklist"
    CALL h5dread_f(dset,H5T_NATIVE_INTEGER,iblock,data_dims,hdferr)
    do ii=1,nblocks
       write(block_list(ii),'(I10)')iblock(ii)
    end do
  else
    write(*,*)"Found string-type blocklist"
    CALL h5dread_f(dset,memtype,block_list,data_dims,hdferr)
  end if

  ! Close dataset and dataspace
  call h5tclose_f(memtype,hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)
 
end subroutine

!!
!!-------------------------------------------------------------------
!!
subroutine wqhdf_getinitialdate(init_datum,start_date,file)
! Subroutine to extract inital date from HDF-file
! Kai Rasmus 2015
  use hdf5
  !use iso_c_binding

  implicit none

!Input arguments
  INTEGER(HID_T) :: file
  character(len=10) :: init_datum ! Simulation intialisation date
  integer, dimension(6) :: start_date ! Simulation start date

  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: size,hdferr
  INTEGER(HID_T)  :: filetype,memtype ! Handle
  INTEGER(HSIZE_T), DIMENSION(2) :: maxdims
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

!
! Open dataset "Init_datum"
!
if(1.gt.2) then ! N0 initial date in HD-file generated by ncin
  CALL h5dopen_f (file,"Init_datum", dset, hdferr)

  ! Get the datatype and its size.
  CALL H5Dget_type_f(dset, filetype, hdferr)
  CALL H5Tget_size_f(filetype, size, hdferr)

 ! Get dataspace
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims,maxdims, hdferr)

  ! Create the memory datatype.
! From https://www.hdfgroup.org/hdf5-quest.html#str1
! HDF5 has a string datatype (H5T_C_S1) that can be modified to create a fixed length
! or variable length array of characters (a string).
!
! Using this datatype without modifying it will create a one-character string. To create
! a string longer than one character, you must get a copy or instance of the datatype,
! and modify it with the H5Tset_size API.

  CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
  CALL H5Tset_size_f(memtype, size, hdferr)

! Read initial date using the default properties.
  CALL h5dread_f(dset,memtype,init_datum,data_dims,hdferr)

  ! Close dataset and dataspace
  !
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)

end if

! Get simulation start date
  CALL h5dopen_f (file,"Simu_start_datum", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)

  ! Read block list using the default properties.
  CALL h5dread_f(dset,H5T_NATIVE_INTEGER,start_date,data_dims,hdferr)

  ! Close dataset and dataspace
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)

end subroutine



subroutine wqhdf_getradiation(ntims,f_itot,file)
! Subroutine to extract radiation
! Kai Rasmus 2015
  use hdf5
  use prec
  implicit none

! Input arguments
  integer(HID_T) :: file
  real(kind=dp), dimension(ntims) :: f_itot ! Total radiation
  integer :: ntims
  real(kind=dp), allocatable,dimension(:) :: inputbuffer
  integer :: ii ! Loop counter

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: hdferr

  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of blocks in block list
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  integer,dimension(1) :: start,count

  CALL h5dopen_f (file,"Solar_radiation", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

  allocate(inputbuffer(data_dims(1)))

  start=(/0/)
  count=(/ntims/)
  CALL h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, &
                             start,count,hdferr)

  ! Read block list using the default properties.
  CALL h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer,data_dims, hdferr)

  do ii=1,ntims
     f_itot(ii)=inputbuffer(ii)
  end do
  
  deallocate(inputbuffer)

  ! Close dataset and dataspace
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)

end subroutine

subroutine wqhdf_getdoublevector(ndim,vector,dataset,file)
! Subroutine to extract 1-dimensional double vector from file
! Kai Rasmus 2015
  use hdf5
  use prec
  implicit none

! Input arguments
  integer(HID_T) :: file
  real(kind=dp), dimension(ndim) :: vector !
  integer :: ndim
  character(len=40) :: dataset

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: hdferr

  INTEGER(HSIZE_T), dimension(2) :: blndims ! Number of blocks in block list
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  CALL h5dopen_f (file,dataset, dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

  ! Read block list using the default properties.
  CALL h5dread_f(dset,H5T_NATIVE_DOUBLE,vector,data_dims, hdferr)

  ! Close dataset and dataspace
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)

end subroutine



subroutine  wqhdf_getadvection(ntims,nifaces,nlayers,f_avl,f_avr,f_diff,s_inter_list,file)
! Subroutine to get advection
! Kai Rasmus 2015

  use hdf5
  use prec
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: ntims,nlayers,nifaces ! Dimensions
  real(kind=dp), dimension(nifaces,ntims,nlayers) :: f_avl ! Interface left advection
  real(kind=dp), dimension(nifaces,ntims,nlayers) :: f_avr ! Interface right advection
  real(kind=dp), dimension(nifaces,ntims,nlayers) :: f_diff ! Interface diffusion
  character(len=21), dimension(nifaces) :: s_inter_list ! List of interfaces in string format

! Handles
  INTEGER(HID_T) :: space,dset,space1,dset1,space2,dset2
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group,subgroup

  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
  real(kind=dp),dimension(nlayers,ntims) :: inputbuffer,inputbuffer1,inputbuffer2

  INTEGER :: II,jj,kk ! LOOP COUNTERs

  CALL h5gopen_f(file,"Interfaces", group, hdferr)
 ! Advection
  do ii=1,nifaces
        CALL h5gopen_f(group,s_inter_list(ii),subgroup,hdferr)

        CALL h5dopen_f (subgroup,"Av_l", dset, hdferr)
        CALL h5dopen_f (subgroup,"Av_r", dset1, hdferr)
        CALL h5dopen_f (subgroup,"ADd", dset2, hdferr)

        CALL h5dget_space_f(dset, space, hdferr)
        CALL h5dget_space_f(dset1, space1, hdferr)
        CALL h5dget_space_f(dset2, space2, hdferr)

!Left
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer,data_dims,hdferr)
!Right
        call h5dread_f(dset1,H5T_NATIVE_DOUBLE,inputbuffer1,data_dims,hdferr)
!Diffusion
        call h5dread_f(dset2,H5T_NATIVE_DOUBLE,inputbuffer2,data_dims,hdferr)

        do jj=1,ntims
          do kk=1,nlayers
                f_avl(ii,jj,kk)=inputbuffer(kk,jj)
                f_avr(ii,jj,kk)=inputbuffer1(kk,jj)
                f_diff(ii,jj,kk)=inputbuffer2(kk,jj)
          end do
        end do

        ! Close subgroup and dataset
        CALL h5dclose_f(dset, hdferr)
        CALL h5dclose_f(dset1, hdferr)
        CALL h5dclose_f(dset2, hdferr)

        CALL h5sclose_f(space, hdferr)
        CALL h5sclose_f(space1, hdferr)
        CALL h5sclose_f(space2, hdferr)
        CALL h5gclose_f(subgroup, hdferr)

  end do

  call h5gclose_f(group,hdferr)


end subroutine

subroutine wqhdf_getinitcond(ntims,nblocks,nlayers,svar,offset,&
           block_list,file)
! Subroutine to read initial value from input file
! Kai Rasmus 2015

  use hdf5
  use prec
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: ntims      ! Number of timesteps
  integer :: nblocks    ! Number of blocks
  integer :: nlayers    ! Number of layers
  real(kind=dp), dimension(ntims,nblocks,nlayers) :: svar ! Variable array for initial values
  integer :: offset     ! Row number of value in input file
  character(len=10),dimension(nblocks) :: block_list ! List of blocks

! Handles
  INTEGER(HID_T) :: dset
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group, subgroup

  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
! TODO: Note that the inputbuffer needs to be changed to 
! reflect the possibility of cahnging the naumber of tracers
  real(kind=dp),dimension(nlayers,8) :: inputbuffer
  real(kind=dp),dimension(nlayers*8) :: inputbuffera

  INTEGER :: ii,jj,jj2,kk               ! Loop counters
  character(len=128) :: blockname       ! String for block name

  CALL h5gopen_f(file,"Blocks", group, hdferr)

 ! Initial conditions
 ! Default value is 0.0 for everything
  svar=0.0

  blockloop: DO ii=1,nblocks
!   if(block_list(ii).gt."100000")then ! Not a boundary conditions block [COMMENTE OUT AS OBSOLETE 20190828/JR]
   !    write(blockname,'(I0)')block_list(ii)
    blockname=block_list(ii)
   !write(0,*)ii,blockname
    CALL h5gopen_f(group,blockname,subgroup,hdferr)
    CALL h5dopen_f (subgroup,"Init_inside", dset, hdferr)
    IF (hdferr.ne.0) THEN
       write(*,'(A,A12)') "No Init_inside subgroup in group. All initial values set to 0.0 for block: ", blockname
    ELSE
       ! CALL h5dget_space_f(dset, space, hdferr)
       call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffera,data_dims,hdferr)

       !TODO: what to do with -1.0 values??
       ! Set them to zero here
       kk=1
       do jj2=1,8
        do jj=1,nlayers
          inputbuffer(jj,jj2)=max(inputbuffera(kk),0.0)
          kk=kk+1
          ! inputbuffer(2,jj)=max(inputbuffer(2,jj),0.0)
        end do
       end do

       svar(1,ii,:)=inputbuffer(:,offset)

       !write(*,*)inputbuffer
       !call h5sclose_f(space,hdferr)

       CALL h5dclose_f(dset, hdferr)

    END IF

    CALL h5gclose_f(subgroup, hdferr)

!   END IF

  END DO blockloop

  call h5gclose_f(group,hdferr)
end subroutine


subroutine wqhdf_getblockproperties(ntims,nblocks,nlayers,b_volume,b_area,&
                                b_ad,b_awdw,b_awup,b_s,b_t,b_tbot,b_ubot,&
                                block_list,file)
! Subroutine to get block properties
! Kai Rasmus 2015

  use hdf5
  use prec
  use wqficos_helper
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: ntims,nblocks,nlayers ! Dimensions  
  real(kind=dp),dimension(ntims,nblocks,nlayers+1),intent(out) :: b_volume ! Volume of block
  real(kind=dp),dimension(ntims,nblocks,nlayers) :: b_area  ! Area of block
  real(kind=dp),dimension(ntims,nblocks) :: b_ad            ! Diffusion within block (TODO: vertical diffusion??)
  real(kind=dp),dimension(ntims,nblocks) :: b_awdw          ! Advection down
  real(kind=dp),dimension(ntims,nblocks) :: b_awup          ! Advection up
  real(kind=dp),dimension(ntims,nblocks,nlayers) :: b_s     ! Salinity
  real(kind=dp),dimension(ntims,nblocks,nlayers) :: b_t     ! Temperature
  real(kind=dp),dimension(ntims,nblocks) :: b_tbot          ! Bottom temperature
  real(kind=dp),dimension(ntims,nblocks) :: b_ubot          ! Bottom velocity
  character(len=10),dimension(nblocks) :: block_list ! List of blocks

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group,subgroup

  INTEGER(HSIZE_T), dimension(2) :: blndims
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
  real(kind=dp),dimension(nlayers) :: inputbuffer1
  real(kind=dp),dimension(ntims) :: inputbuffer2
  real(kind=dp),dimension(ntims*nlayers) :: inputbuffer3a
  real(kind=dp),dimension(ntims*(nlayers+1)) :: inputbuffer4a


  INTEGER :: ii,jj,kk,ll,nn                   ! Loop counters
  character(len=128) :: blockname       ! String for block name

  integer :: iblock, bindex

  CALL h5gopen_f(file,"Blocks", group, hdferr)

 ! Block properties

  do ii=1,nblocks
     blockname=block_list(ii)
     read(blockname,'(I10)')iblock

     call getblockindex(nblocks,block_list,iblock,bindex)

! Open block
     CALL h5gopen_f(group,blockname,subgroup,hdferr)
! AD
     CALL h5dopen_f (subgroup,"AD", dset, hdferr)
     call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer2,data_dims,hdferr)
     b_ad(:,bindex)=inputbuffer2
     call h5dclose_f(dset,hdferr)

! Area
     CALL h5dopen_f (subgroup,"Area", dset, hdferr)
     call h5dget_space_f(dset,space,hdferr)
     CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

     if(data_dims(2).eq.1)then
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer1,data_dims,hdferr)
     else
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer3a,(/data_dims(1),data_dims(2)/),hdferr)
     end if

     ll=1
     do jj=1,ntims
        do kk=1,nlayers
           if(data_dims(2).eq.1)then
             b_area(jj,bindex,kk)=inputbuffer1(kk)
           else
             b_area(jj,bindex,kk)=inputbuffer3a(ll)
             ll=ll+1
           end if
        end do
     end do
     call h5sclose_f(space,hdferr)
     call h5dclose_f(dset,hdferr)

! Aw_dw
     CALL h5dopen_f (subgroup,"Aw_dw", dset, hdferr)
     call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer2,data_dims,hdferr)
     b_awdw(:,bindex)=inputbuffer2
     call h5dclose_f(dset,hdferr)

! Aw_up
     CALL h5dopen_f (subgroup,"Aw_up", dset, hdferr)
     call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer2,data_dims,hdferr)
     b_awup(:,bindex)=inputbuffer2
     call h5dclose_f(dset,hdferr)

! Sl
     CALL h5dopen_f (subgroup,"Sl", dset, hdferr)
     call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer3a,data_dims,hdferr)
     ll=1
     do jj=1,ntims
       do kk=1,nlayers
         b_s(jj,bindex,kk)=inputbuffer3a(ll)
         ll=ll+1
       end do
     end do
     call h5dclose_f(dset,hdferr)

! Tp
     CALL h5dopen_f (subgroup,"Tp", dset, hdferr)
     call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer3a,data_dims,hdferr)
     ll=1
     do jj=1,ntims
        do kk=1,nlayers
          b_t(jj,bindex,kk)=inputbuffer3a(ll)
          ll=ll+1
        end do
     end do
     call h5dclose_f(dset,hdferr)

! Tp_bottom
     CALL h5dopen_f (subgroup,"Tp_bottom", dset, hdferr)
     call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer2,data_dims,hdferr)
     b_tbot(:,bindex)=inputbuffer2
     call h5dclose_f(dset,hdferr)

! V_bottom
     CALL h5dopen_f (subgroup,"V_bottom", dset, hdferr)
     call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer2,data_dims,hdferr)
     b_ubot(:,bindex)=inputbuffer2
     call h5dclose_f(dset,hdferr)

! Volume
     CALL h5dopen_f (subgroup,"Volume", dset, hdferr)
     call h5dget_space_f(dset,space,hdferr)
     CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

     if(data_dims(2).eq.1)then
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer1,data_dims,hdferr)
     else
        if(data_dims(1).eq.nlayers+1)then
           call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer4a,data_dims,hdferr)
        else 
           call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer3a,data_dims,hdferr)
           if(inputbuffer3a(1).eq.0.0)then
              write(*,*)inputbuffer3a
           end if
        end if
     end if

     nn=1
     ll=1
     do jj=1,ntims
        do kk=1,nlayers
           if(data_dims(2).eq.1)then
              b_volume(jj,bindex,kk)=inputbuffer1(kk)
           else
              if(data_dims(1).eq.nlayers+1)then
                 b_volume(jj,bindex,kk)=inputbuffer4a(nn)
                 nn=nn+1
              else
                 b_volume(jj,bindex,kk)=inputbuffer3a(ll)
                 b_volume(jj,bindex,nlayers+1)=0.0
                 ll=ll+1
              end if
           end if
        end do
        nn=nn+1
     end do

     call h5sclose_f(space,hdferr)
     call h5dclose_f(dset,hdferr) 
     CALL h5gclose_f(subgroup,hdferr)
  end do

  call h5gclose_f(group,hdferr)

  do ii=1,nblocks
     if(b_volume(1,ii,1).eq.0.0)then
        write(*,*)"VOLUME FAIL",nblocks,ii,block_list(bindex)
        stop
     end if
  end do
 
end subroutine ! Get block properties

!Not used anywhere. Commented out 2019-08-28/JR
!subroutine wqhdf_getifacearea(nifaces,nlayers,ntims,iface_area,s_inter_list,file)
!! Subroutine to get interface areas from interface list
!! Kai Rasmus 2015
!
!  use hdf5
!  implicit none
!
!! Input arguments
!  integer(HID_T) :: file
!  integer :: nifaces,nlayers,ntims ! Dimensions
!  integer,dimension(ntims,nifaces,nlayers) :: iface_area  ! Interface area
!  character(len=21), dimension(nifaces) :: s_inter_list ! List of interfaces in string format
!
!! Handles
!  INTEGER(HID_T) :: space,dset
!  INTEGER(SIZE_T) :: hdferr
!  INTEGER(HID_T)  :: group,subgroup
!
!  INTEGER(HSIZE_T), dimension(2) :: blndims
!  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
!  real(kind=dp),allocatable,dimension(:,:) :: inputbuffer
!
!  integer :: ii ! Loop counter
!
!  allocate(inputbuffer(ntims,nlayers))
!
!  CALL h5gopen_f(file,"Interfaces", group, hdferr)
!  do ii=1,nifaces
!        CALL h5gopen_f(group,s_inter_list(ii),subgroup,hdferr)
!        CALL h5dopen_f(subgroup,"A", dset, hdferr)
!        CALL h5dget_space_f(dset, space, hdferr)
!        CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)
!        call h5dread_f(dset,H5T_NATIVE_DOUBLE,inputbuffer,data_dims,hdferr)
!        call h5sclose_f(space,hdferr)
!        call h5dclose_f(dset,hdferr)
!        call h5gclose_f(subgroup,hdferr)
!  end do
!
!  call h5gclose_f(group,hdferr)
!  deallocate(inputbuffer)
!
!end subroutine


!!
subroutine wqhdf_createoutput(n_blocks,n_tims,block_list,file)
! Subroutine to create output groups and datasets
! Kai Rasmus 2015

  use hdf5
  USE ISO_C_BINDING
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: n_blocks  ! Number of blocks
  character(len=10),dimension(n_blocks) :: block_list
  integer :: n_tims    ! Number of timesteps

! Handles
  INTEGER(HID_T) :: space,space2,dset
  INTEGER(SIZE_T) :: size,hdferr
  INTEGER(HID_T)  :: group,filetype,memtype,subgroup
  INTEGER(HSIZE_T), dimension(2) :: blndims
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  integer :: ii ! Loop counters
  character(len=128) :: blockname     ! String for block name
  CHARACTER(LEN=9) :: aname   ! Attribute name
  CHARACTER(len=7), DIMENSION(1) ::  attr_data

  INTEGER(HSIZE_T), DIMENSION(1) :: dims != (/dim0, dim1/)
  integer(hsize_t) :: irank

 !dims(1)=1
 dims(1)=n_tims
 irank=1

  CALL h5gcreate_f(file, "Blocks", group, hdferr)
  call h5tvlen_create_f(H5T_FORTRAN_S1,memtype,hdferr)
  CALL h5tvlen_create_f(H5T_STD_I32LE, filetype, hdferr)
  CALL h5screate_f(H5S_NULL_F, space2, hdferr)

  CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
  size=7
  CALL H5Tset_size_f(memtype, size, hdferr)

  blndims(1)=1

  do ii=1,n_blocks
    blockname=block_list(ii)
! Create block-group
    call h5gcreate_f(group,blockname,subgroup,hdferr)

! Create space and dataset for variables
    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Cdet_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Cdet_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Cv",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Ndet_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Ndet_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Nv",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^2"
    call wqhdf_cattr(dset,"Unit"," mg/m^2",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Pdet_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Pdet_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Pfev",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^2"
    call wqhdf_cattr(dset,"Unit"," mg/m^2",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"Pv",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^2"
    call wqhdf_cattr(dset,"Unit"," mg/m^2",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"cA",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mgN/m^3"
    call wqhdf_cattr(dset,"Unit","mgN/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"cC",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mgN/m^3"
    call wqhdf_cattr(dset,"Unit","mgN/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"cDIN_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"cDIN_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"cDIP_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"cDIP_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"height_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="m     "
    call wqhdf_cattr(dset,"Unit","      m",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"height_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="m"
    call wqhdf_cattr(dset,"Unit","      m",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"totN_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"totN_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"totP_0",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"totP_1",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mg/m^3"
    call wqhdf_cattr(dset,"Unit"," mg/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)
    
    ! ---------------------------------------------------
    ! Karel Kaurila 2024 - N2 fixation
    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"rN2fixFC",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="tn N"
    call wqhdf_cattr(dset,"Unit"," tn   N",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)

    CALL h5screate_simple_f(irank,dims, space, hdferr)
    CALL h5dcreate_f(subgroup,"cN2fixFC",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
    ! Create attribute
    aname="Unit"
    attr_data(1)="mgN/m^3"
    call wqhdf_cattr(dset,"Unit","mgN/m^3",memtype)
    CALL h5dclose_f(dset , hdferr)
    call h5sclose_f(space,hdferr)
    ! -----------------------------------------------------
    
    call h5gclose_f(subgroup,hdferr)
  end do

  ! Date group
  data_dims(1)=6
  data_dims(2)=n_tims
  CALL h5screate_simple_f(2,data_dims, space, hdferr)
  CALL h5dcreate_f(file,"Dates",H5T_NATIVE_INTEGER,&
         space, dset, hdferr)
  CALL h5dclose_f(dset,hdferr)
  call h5sclose_f(space,hdferr)

  ! Time array
  CALL h5screate_simple_f(1,blndims, space, hdferr)
  CALL h5dcreate_f(file,"Time_array",H5T_NATIVE_DOUBLE,&
         space, dset, hdferr)
  CALL h5dclose_f(dset,hdferr)
  call h5sclose_f(space,hdferr)

  ! Wq model
  data_dims(1)=128
  data_dims(2)=0
  CALL h5screate_simple_f(1,data_dims, space, hdferr)
  CALL h5dcreate_f(file,"Wq_model",memtype,&
         space, dset, hdferr)
  CALL h5dclose_f(dset,hdferr)
  call h5sclose_f(space,hdferr)

  call h5sclose_f(space2,hdferr)
  CALL h5tclose_f(memtype, hdferr)
  CALL h5tclose_f(filetype, hdferr)
  call h5gclose_f(group,hdferr)

end subroutine

subroutine wqhdf_getnbndblocks(nbndblocks,file)
! Subroutine to get number of boundary blocks
! Kai Rasmus 2015

  use hdf5
  USE ISO_C_BINDING
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: nbndblocks ! Number of boundary blocks

! Handles
  INTEGER(HID_T) :: space,dset
  INTEGER(SIZE_T) :: hdferr

  INTEGER(HSIZE_T), dimension(2) :: blndims
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  CALL h5dopen_f (file,"Boundary_block_list", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

  ! Close dataset and dataspace
  CALL h5sclose_f(space, hdferr)
  CALL h5dclose_f(dset, hdferr)

  nbndblocks=blndims(1)


end subroutine

subroutine wqhdf_getbnddata(ntims,nbndblocks,nparams,bndblocks,bnddata,file)
! Subroutine to get number of boundary block data
! Kai Rasmus 2015

  use hdf5
  use prec
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: nbndblocks                       ! Number of boundary blocks
  integer,dimension(nbndblocks) :: bndblocks  ! List of boundary blocks
  integer :: nparams                          ! Number of parameters
  integer :: ntims                            ! Number of timesteps
  real(kind=dp),dimension(nbndblocks,ntims,nparams) :: bnddata ! Array for boundary condition data

  real(kind=dp),dimension(nparams,ntims) ::bufferbnddata  

! Handles
  INTEGER(HID_T) :: space,dset,space2,dset2
  INTEGER(SIZE_T) :: size,hdferr,filetype2
  INTEGER(HID_T)  :: group,filetype,memtype,subgroup

  INTEGER(HSIZE_T), dimension(2) :: blndims,blndims2
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims,data_dims2

  integer :: ii,jj,kk ! Loop counters
  character(len=128) :: bndbuffer

  character(len=10),dimension(nbndblocks) :: sbndblocks

  CALL h5dopen_f (file,"Boundary_block_list", dset, hdferr)
  CALL h5dget_space_f(dset, space, hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

  ! Get the datatype and its size.
  CALL H5Dget_type_f(dset, filetype, hdferr)
  ! Get filetype type (H5T_INTEGER_F or H5T_STRING_F)
  call h5tget_class_f(filetype,filetype2,hdferr)
  if (filetype2.eq.H5T_INTEGER_F) then
     write(*,*)"Integer-type boundary block list found"
  else
     write(*,*)"String-type boundary block list found"
  end if

  ! Read boundary block list using the default properties.
  CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
  size=10
  CALL H5Tset_size_f(memtype,size, hdferr)

! Read initial date using the default properties.
  if(filetype2.eq.H5T_INTEGER_F)then
      write(*,*)"Integer-type boundary block list found"
      CALL h5dread_f(dset,H5T_NATIVE_INTEGER,bndblocks,data_dims, hdferr)
  else
      write(*,*)"String-type boundary block list found"
      CALL h5dread_f(dset,memtype,sbndblocks,data_dims, hdferr)
      do jj=1,nbndblocks
         read(sbndblocks(jj),*)bndblocks(jj)
      end do
  end if

  call h5tclose_f(memtype,hdferr)

  CALL h5gopen_f(file,"Blocks", group, hdferr)
  do ii=1,nbndblocks
        write(bndbuffer,'(I0)')bndblocks(ii)
        CALL h5gopen_f(group,bndbuffer,subgroup,hdferr)
        CALL h5dopen_f(subgroup,"BC_vect", dset2, hdferr)
        CALL h5dget_space_f(dset2, space2, hdferr)
        CALL h5sget_simple_extent_dims_f(space2,data_dims2, blndims2, hdferr)
        call h5dread_f(dset2,H5T_NATIVE_DOUBLE,bufferbnddata,data_dims2,hdferr)
        do jj=1,ntims
          do kk=1,nparams
             bnddata(ii,jj,kk)=bufferbnddata(kk,jj) ! Note that indii are reversed
          end do
        end do
        call h5sclose_f(space2,hdferr)
        call h5dclose_f(dset2,hdferr)
        call h5gclose_f(subgroup,hdferr)
  end do

  ! Close dataset and dataspace
  CALL h5dclose_f(dset, hdferr)
  CALL h5gclose_f(group, hdferr)
  CALL h5sclose_f(space,hdferr)

end subroutine


subroutine wqhdf_loadings(ntims,nblocks,nlayers,nparams,block_list,loading,file)
! Subroutine to load loadings from file
!
! A global array of dataset names is used to store the names todether with global
! values that describe how many sources and how long their names are.
!
! Kai Rasmus 2015
!
  use hdf5
  USE ISO_C_BINDING
  USE g_iterate
  use prec
  use wqficos_helper
  implicit none

! Input arguments
  integer(HID_T) :: file
  integer :: ntims,nblocks,nlayers,nparams
  character(len=10),dimension(nblocks) :: block_list
  real(kind=dp),dimension(ntims,nblocks,nlayers,nparams) :: loading ! Point source loading

! Handles
  INTEGER(HID_T)  :: space2,dset2,layerattr,cyclicattr
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group,subgroup,blockgroup

  INTEGER(HSIZE_T), dimension(2) :: blndims2
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims2

  integer :: ii,jj,kk,ll,mm,hh ! Loop counters
  character(len=128) :: s_block   ! String for block names

  INTEGER :: nlinks     ! Number of links in group
  INTEGER :: max_corder ! Current maximum creation order value for group

  integer :: storage_type

  integer :: iblock

  TYPE(C_FUNPTR) :: funptr ! Variables related to the callback function
  TYPE(C_PTR) :: ptr
  INTEGER(hsize_t) :: idx
  INTEGER :: ret_value

  INTEGER,dimension(1),target :: rdata,cdata ! Read buffer for layer attribute and cyclical atrribute
  TYPE(C_PTR) :: f_ptr

  real(kind=dp),allocatable,dimension(:,:) :: databuffer
  real(kind=dp),allocatable,dimension(:) :: databuffera

  integer :: bindex ! Block index

! Zero loadings
  loading=0.0

  idx = 0
  funptr = C_FUNLOC(op_func) ! call back function
  ptr    = C_NULL_PTR

 ! Open block-list
  CALL h5gopen_f(file,"Blocks", group, hdferr)

  do ii=1,nblocks

    s_block=block_list(ii)
    read(s_block,'(I10)')iblock

    if(s_block(1:1).ne."8".and.iblock.gt.15000)then ! Skip boundary blocks
       call h5gopen_f(group,s_block,blockgroup,hdferr)

       call h5gopen_f(blockgroup,"Loads",subgroup,hdferr)

       call h5gget_info_f(subgroup,storage_type,nlinks,max_corder,hdferr)
       write(*,*)
       write(*,'(A)')"--------------------------------------------------------"
       write(*,'(A7A10)')"Block: ",s_block(1:10)
       write(*,'(A17I3)')"Number of loads: ",nlinks
       write(*,*)

       idx=0

       loading_stations=""
       numstations=0
       call h5literate_f(subgroup,H5_INDEX_NAME_F,H5_ITER_NATIVE_F,&
                         idx,funptr,ptr,ret_value,hdferr)

       do jj=1,numstations
          write(*,'(A16A)')"Loading source: ",&
                         loading_stations(jj)(1:stationlens(jj))
! Open station
          CALL h5dopen_f(subgroup,loading_stations(jj)(1:stationlens(jj)+1),&
                         dset2, hdferr)

! Get Cyclic attribute
          CALL h5aopen_f(dset2,'Cyclical',cyclicattr, hdferr)
          f_ptr = C_LOC(cdata(1))
          CALL h5aread_f(cyclicattr, H5T_NATIVE_INTEGER,f_ptr, hdferr)
          write(*,'(AI1)')"Cyclical: ",cdata(1)

! Get Layer attribute
          CALL h5aopen_f(dset2,'Layer',layerattr, hdferr)
          f_ptr = C_LOC(rdata(1))
          CALL h5aread_f(layerattr, H5T_NATIVE_INTEGER,f_ptr, hdferr)
          write(*,'(AI1)')"Layer: ",rdata(1)

! Get data
          CALL h5dget_space_f(dset2, space2, hdferr)
          CALL h5sget_simple_extent_dims_f(space2,data_dims2, blndims2, hdferr)

          allocate(databuffer(data_dims2(1),data_dims2(2)))
          allocate(databuffera(data_dims2(1)*data_dims2(2)))
          call h5dread_f(dset2,H5T_NATIVE_DOUBLE,databuffera,data_dims2,hdferr)

! Put data in correct arrays
          call getblockindex(nblocks,block_list,iblock,bindex)
          hh=1
          mm=1
          write(*,'(A)')"Point source loading"
          if (cdata(1).eq.1) then
             ! cyclic loading
             do kk=1,ntims
              do ll=1,nparams
                 loading(kk,bindex,rdata+1,ll)=loading(kk,bindex,rdata+1,ll)&
                          +databuffera(hh)!databuffer(ll,mm)
                 mm=mm+1
                 hh=hh+1
                 if (hh.eq.data_dims2(1)*data_dims2(2)-1) then
                    mm=1
                    hh=1
                 end if
              end do
             end do
          else
             do kk=1,ntims
              do ll=1,nparams
                 loading(kk,bindex,rdata+1,ll)=loading(kk,bindex,rdata+1,ll)&
                          +databuffera(hh)!databuffer(ll,kk)
                 hh=hh+1
              end do
             end do
          end if

          write(*,*)
          deallocate(databuffer)
          deallocate(databuffera)

! Release resources
          call h5aclose_f(layerattr,hdferr)
          call h5aclose_f(cyclicattr,hdferr)
          call h5sclose_f(space2,hdferr)
          call h5dclose_f(dset2,hdferr)
       end do

       call h5gclose_f(subgroup,hdferr)
       call h5gclose_f(blockgroup,hdferr)

    end if

  end do

  call h5gclose_f(group,hdferr)

end subroutine


subroutine hdfwq_writedate(model_date,tt,outfile)
! Subroutine to write a date to the output file
! 
! Called by main program
!
! Kai Rasmus 2015
  use hdf5
  USE ISO_C_BINDING
  USE g_iterate
  implicit none

! Input arguments
  integer(HID_T) :: outfile             ! Outputfile
  integer :: tt                         ! Timestep index, number of timesteps
  integer,dimension(6) :: model_date    ! Date vector

! Handles
  INTEGER(HID_T) :: space,dset,memspace
  INTEGER(SIZE_T) :: hdferr

! Dimensions
  INTEGER(HSIZE_T), dimension(2)   :: blndims
  INTEGER(HSIZE_T), DIMENSION(2)   :: data_dims
  INTEGER(HSIZE_T), DIMENSION(1:2) :: start, count

! Date group
  data_dims(1)=1
  data_dims(2)=6


  call h5dopen_f(outfile,"Dates",dset,hdferr)
  CALL h5dget_space_f(dset,space,hdferr)
  CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)
  CALL h5screate_simple_f(1,data_dims, memspace,hdferr)
  start=(/0,tt-1/)
  count=(/6,1/)
  CALL h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, &
                             start,count,hdferr)
  CALL h5dwrite_f(dset,H5T_NATIVE_INTEGER,model_date,(/1,6/),&
                hdferr,memspace,file_space_id=space)
  call h5sclose_f(memspace,hdferr)  
  call h5sclose_f(space,hdferr)
  call h5dclose_f(dset,hdferr)

end subroutine

subroutine hdfwq_writeoutput(varname,s_cc,nblocks,block_list,tt,outfile)
! Subroutine to write a vector the output file
!
! Called by the main program
!
! Kai Rasmus 2015

! Note to self:
! The Fortran HDF interface sucks!!

  use hdf5
  USE ISO_C_BINDING
  USE g_iterate
  USE prec
  implicit none

! Input arguments
  integer(HID_T) :: outfile             ! Outputfile
  integer :: nblocks                    ! Number of blocks 
  integer :: tt                         ! Timestep index
  character(LEN=*) :: varname           ! Name of parametr to write
  character(LEN=10),dimension(nblocks) :: block_list
  real(kind=dp),dimension(nblocks) :: s_cc     ! Vector to write to file

! Handles
  INTEGER(HID_T) :: space,dset,memspace
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group,blockgroup

! Dimensions
  INTEGER(HSIZE_T), dimension(2) :: blndims
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  integer :: ii ! Loop counters
  INTEGER(HSIZE_T), DIMENSION(1)   :: start, count

  real(kind=dp) :: datapoint

  character(len=128) :: s_block ! Character string for block name

  CALL h5gopen_f(outfile,"Blocks", group, hdferr)

  do ii=1,nblocks
        s_block=block_list(ii)
        call h5gopen_f(group,s_block,blockgroup,hdferr)
        call h5dopen_f(blockgroup,varname,dset,hdferr)
        CALL h5dget_space_f(dset,space,hdferr)
        CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

        CALL h5screate_simple_f(1,(/1/), memspace,hdferr)
        start=(/tt-1/)
        count=(/1/)
        CALL h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, &
                                   start,count,hdferr)

        datapoint=s_cc(ii)

        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,datapoint,(/1/),&
                        hdferr,memspace,file_space_id=space)

        call h5sclose_f(memspace,hdferr)
        call h5sclose_f(space,hdferr)
        call h5dclose_f(dset,hdferr)
        call h5gclose_f(blockgroup,hdferr)
  end do
  call h5gclose_f(group,hdferr)
end subroutine


subroutine wqhdf_cattr(dset,aname,attr_data,memtype)

  use hdf5
  implicit none

! Input arguments
  CHARACTER(LEN=4) :: aname   ! Attribute name
  CHARACTER(len=7), DIMENSION(1) ::  attr_data  
  INTEGER(HID_T) :: dset,memtype

! Handles
  INTEGER(HID_T) :: aspace,attr
  INTEGER(SIZE_T) :: hdferr

    CALL h5screate_simple_f(1,(/1,1/), aspace, hdferr)
    CALL h5acreate_f(dset,aname,memtype, aspace, attr, hdferr)
    CALL h5awrite_f(attr, memtype,attr_data,(/1,7/),hdferr)
    CALL h5aclose_f(attr , hdferr)
    call h5sclose_f(aspace,hdferr)

end subroutine

!!
subroutine hdfwq_writethickness(nblocks,b_thick_top,b_thick_bottom,&
                                block_list,tt,outfile)
! Subroutine to write thicknesses to HDF-file
! Kai Rasmus, 2015

  use hdf5
  use prec
  implicit none

! Input argumetns
integer :: tt    ! time step index
integer(HID_T) :: outfile ! File identifier for output file
integer :: nblocks ! Number of blocks
character(len=10),dimension(nblocks) :: block_list
real(kind=dp),dimension(nblocks) :: b_thick_top,b_thick_bottom ! Thickness of layers

! Handles
  INTEGER(HID_T) :: space,dset,memspace
  INTEGER(SIZE_T) :: hdferr
  INTEGER(HID_T)  :: group,blockgroup

! Local variables
  integer :: ii ! Loop counters
  INTEGER(HSIZE_T), DIMENSION(1)   :: start, count
  real(kind=dp) :: datapoint

  character(len=128) :: s_block ! Character string for block name

! Dimensions
  INTEGER(HSIZE_T), dimension(2) :: blndims
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims

  CALL h5gopen_f(outfile,"Blocks", group, hdferr)

  do ii=1,nblocks
        s_block=block_list(ii)
        call h5gopen_f(group,s_block,blockgroup,hdferr)
        call h5dopen_f(blockgroup,"height_0",dset,hdferr)
        CALL h5dget_space_f(dset,space,hdferr)
        CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

        CALL h5screate_simple_f(1,(/1/), memspace,hdferr)
        start=(/tt-1/)
        count=(/1/)
        CALL h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, &
                start,count,hdferr)

        datapoint=b_thick_top(ii)

        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,datapoint,(/1/),&
                hdferr,memspace,file_space_id=space)

        call h5sclose_f(memspace,hdferr)
        call h5sclose_f(space,hdferr)
        call h5dclose_f(dset,hdferr)


        call h5dopen_f(blockgroup,"height_1",dset,hdferr)
        CALL h5dget_space_f(dset,space,hdferr)
        CALL h5sget_simple_extent_dims_f(space,data_dims, blndims, hdferr)

        CALL h5screate_simple_f(1,(/1/), memspace,hdferr)
        start=(/tt-1/)
        count=(/1/)
        CALL h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, &
                start,count,hdferr)

        datapoint=b_thick_bottom(ii)

        CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,datapoint,(/1/),&
                hdferr,memspace,file_space_id=space)

        call h5sclose_f(memspace,hdferr)
        call h5sclose_f(space,hdferr)
        call h5dclose_f(dset,hdferr)

        call h5gclose_f(blockgroup,hdferr)
  end do

 call h5gclose_f(group,hdferr)

end subroutine

END MODULE hdfio
