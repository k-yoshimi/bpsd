! bpsd_species.f90

module bpsd_species

  use bpsd_kinds
  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_put_species,bpsd_get_species, &
         bpsd_save_species,bpsd_load_species, &
         bpsd_get_species_kdata
  private

  ! Number of scalar fields packed per species in the flat speciesx
  ! buffer. Single source of truth for the put/get/setup_kdata layout;
  ! changing bpsd_species_data must be matched by updating this and
  ! the kid/kunit assignments below.
  integer(ikind), parameter :: nfields = 3   ! pa, pz, npa

  logical, save :: bpsd_speciesx_init_flag = .TRUE.
  type(bpsd_data0Dx_type), save :: speciesx

contains

!-----------------------------------------------------------------------
  subroutine bpsd_init_speciesx
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    speciesx%status=0
    speciesx%dataName='species'
    speciesx%ndmax=0

    bpsd_speciesx_init_flag = .FALSE.

    return
  end subroutine bpsd_init_speciesx

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_setup_species_kdata
!-----------------------------------------------------------------------
    IMPLICIT NONE
    integer(ikind):: nd

    ! Single source of truth for the layout is the `nfields` parameter
    ! at module top. Earlier revisions had ndmax = nsmax*5 with stride
    ! 3, which left 2 trailing slots per species uninitialized and
    ! visible to bpsd_save_data0Dx.
    do nd=0,speciesx%ndmax-nfields,nfields
       speciesx%kid(nd+1)='species%pa'
       speciesx%kid(nd+2)='species%pz'
       speciesx%kid(nd+3)='species%npa'
       speciesx%kunit(nd+1)=' '
       speciesx%kunit(nd+2)=' '
       speciesx%kunit(nd+3)=' '
    enddo
    RETURN
  END SUBROUTINE bpsd_setup_species_kdata

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_species_data(data,n1)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(bpsd_species_data),DIMENSION(:),ALLOCATABLE:: data
    INTEGER(ikind),INTENT(IN):: n1

    IF(ALLOCATED(data)) THEN
       IF(n1.LE.0) THEN
          DEALLOCATE(data)
       ELSE IF(n1.NE.SIZE(data,1)) THEN
          DEALLOCATE(data)
          ALLOCATE(data(n1))
       END IF
    ELSE
       IF(n1.GT.0)  THEN
          ALLOCATE(data(n1))
       ENDIF
    ENDIF
  END SUBROUTINE bpsd_adjust_species_data

!-----------------------------------------------------------------------
  subroutine bpsd_put_species(species_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_species_type),intent(in):: species_in
    integer,intent(out) :: ierr
    integer :: ns,nd

    if(bpsd_speciesx_init_flag) call bpsd_init_speciesx

    speciesx%ndmax=species_in%nsmax*nfields
    CALL bpsd_adjust_karray(speciesx%kid,speciesx%ndmax)
    CALL bpsd_adjust_karray(speciesx%kunit,speciesx%ndmax)
    CALL bpsd_adjust_array1D(speciesx%data,speciesx%ndmax)

    CALL bpsd_setup_species_kdata

    do ns=1,species_in%nsmax
       nd=nfields*(ns-1)
       speciesx%data(nd+1)=species_in%data(ns)%pa
       speciesx%data(nd+2)=species_in%data(ns)%pz
       speciesx%data(nd+3)=species_in%data(ns)%npa
    enddo
    CALL DATE_AND_TIME(speciesx%created_date, &
                       speciesx%created_time, &
                       speciesx%created_timezone)
    speciesx%status = 2
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_put_species'
       do nd=1,speciesx%ndmax
          write(6,'(A32,1PE12.4)') &
                           speciesx%kid(nd),speciesx%data(nd)
       enddo
    endif
    return
  end subroutine bpsd_put_species

!-----------------------------------------------------------------------
  subroutine bpsd_get_species(species_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_kinds
    use bpsd_subs
    implicit none
    type(bpsd_species_type),intent(out) :: species_out
    integer,intent(out) :: ierr
    integer :: nr, nd, ns
    real(dp) :: s
    real(dp), dimension(:), ALLOCATABLE :: v

    if(bpsd_speciesx_init_flag) call bpsd_init_speciesx

    if(speciesx%status.eq.0) then
       write(6,*) 'XX bpsd_get_species: no space allocated to species%data'
       ierr=1
       return
    endif

    if(speciesx%status.eq.1) then
       write(6,*) 'XX bpsd_get_species: no data in species%data'
       ierr=2
       return
    endif

    species_out%nsmax=speciesx%ndmax/nfields

    CALL bpsd_adjust_species_data(species_out%data,species_out%nsmax)

    do ns=1,species_out%nsmax
       nd=nfields*(ns-1)
       species_out%data(ns)%pa =speciesx%data(nd+1)
       species_out%data(ns)%pz =speciesx%data(nd+2)
       species_out%data(ns)%npa=NINT(speciesx%data(nd+3))
    enddo
    speciesx%status=2
    ierr=0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_species'
       do nd=1,speciesx%ndmax
          write(6,'(A32,1PE12.4)') speciesx%kid(nd),speciesx%data(nd)
       enddo
    endif
    return
  end subroutine bpsd_get_species

!-----------------------------------------------------------------------
  subroutine bpsd_save_species(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_speciesx_init_flag) call bpsd_init_speciesx

    if(speciesx%status.gt.1) &
               call bpsd_save_data0Dx(fid,speciesx,ierr)
    return

  end subroutine bpsd_save_species

!-----------------------------------------------------------------------
  subroutine bpsd_load_species(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data0Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: ns,nd

    if(bpsd_speciesx_init_flag) call bpsd_init_speciesx

    speciesx%dataName = datax%dataName
    speciesx%time = datax%time
    speciesx%ndmax=datax%ndmax
    CALL bpsd_adjust_karray(speciesx%kid,speciesx%ndmax)
    CALL bpsd_adjust_karray(speciesx%kunit,speciesx%ndmax)
    CALL bpsd_adjust_array1D(speciesx%data,speciesx%ndmax)
    do nd=1,speciesx%ndmax
       speciesx%data(nd) = datax%data(nd)
    enddo
    do nd=1,speciesx%ndmax
       speciesx%kid(nd)=datax%kid(nd)
       speciesx%kunit(nd)=datax%kunit(nd)
    enddo
    speciesx%created_date = datax%created_date
    speciesx%created_time = datax%created_time
    speciesx%created_timezone = datax%created_timezone
    speciesx%status=2
    ierr=0
    return

  end subroutine bpsd_load_species

!-----------------------------------------------------------------------
  subroutine bpsd_get_species_kdata(ndmax_out,kid_out,kunit_out,ierr)
!-----------------------------------------------------------------------
! Read-only accessor for the internal speciesx kid/kunit metadata.
! Mirrors bpsd_get_trmatrix_kdata so tests can verify the on-array
! layout directly without parsing the unformatted save file.

    use bpsd_subs
    implicit none
    integer,intent(out) :: ndmax_out
    character(len=32),dimension(:),allocatable,intent(out) :: kid_out
    character(len=32),dimension(:),allocatable,intent(out) :: kunit_out
    integer,intent(out) :: ierr
    integer :: nd

    if(bpsd_speciesx_init_flag) call bpsd_init_speciesx

    if(speciesx%status.lt.2) then
       ndmax_out = 0
       ierr = 1
       return
    endif

    ndmax_out = speciesx%ndmax
    if(allocated(kid_out))   deallocate(kid_out)
    if(allocated(kunit_out)) deallocate(kunit_out)
    allocate(kid_out(ndmax_out))
    allocate(kunit_out(ndmax_out))
    do nd = 1, ndmax_out
       kid_out(nd)   = speciesx%kid(nd)
       kunit_out(nd) = speciesx%kunit(nd)
    end do
    ierr = 0
  end subroutine bpsd_get_species_kdata

end module bpsd_species
