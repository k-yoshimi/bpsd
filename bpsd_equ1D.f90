! bpsd_equ1D.f90

module bpsd_equ1D

  use bpsd_kinds
  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_put_equ1D,bpsd_get_equ1D, &
         bpsd_save_equ1D,bpsd_load_equ1D
  private

  logical, save :: bpsd_equ1Dx_init_flag = .TRUE.
  type(bpsd_data1Dx_type), save :: equ1Dx

contains

!-----------------------------------------------------------------------
  subroutine bpsd_init_equ1Dx
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    equ1Dx%status=0
    equ1Dx%dataName='equ1D'
    equ1Dx%ndmax=6

    allocate(equ1Dx%kid(6))
    allocate(equ1Dx%kunit(6))
    equ1Dx%kid(1)='equ1D%psit'
    equ1Dx%kid(2)='equ1D%psip'
    equ1Dx%kid(3)='equ1D%ppp'
    equ1Dx%kid(4)='equ1D%piq'
    equ1Dx%kid(5)='equ1D%pip'
    equ1Dx%kid(6)='equ1D%pit'
    equ1Dx%kunit(1)='Wb'
    equ1Dx%kunit(2)='Wb'
    equ1Dx%kunit(3)='Pa'
    equ1Dx%kunit(4)=' '
    equ1Dx%kunit(5)='A'
    equ1Dx%kunit(6)='A'

    bpsd_equ1Dx_init_flag = .FALSE.

    return
  end subroutine bpsd_init_equ1Dx

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_equ1D_data(data,n1)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(bpsd_equ1D_data),DIMENSION(:),ALLOCATABLE:: data
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
  END SUBROUTINE bpsd_adjust_equ1D_data

!-----------------------------------------------------------------------
  subroutine bpsd_put_equ1D(equ1D_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_equ1D_type):: equ1D_in
    integer :: ierr
    integer :: nr,nd

    if(bpsd_equ1Dx_init_flag) call bpsd_init_equ1Dx

    equ1Dx%nrmax=equ1D_in%nrmax
    CALL bpsd_adjust_array1D(equ1Dx%rho,equ1Dx%nrmax)
    CALL bpsd_adjust_array2D(equ1Dx%data,equ1Dx%nrmax,equ1Dx%ndmax)

    equ1Dx%time=equ1D_in%time
    do nr=1,equ1Dx%nrmax
       equ1Dx%rho(nr)    = equ1D_in%rho(nr)
       equ1Dx%data(nr,1) = equ1D_in%data(nr)%psit
       equ1Dx%data(nr,2) = equ1D_in%data(nr)%psip
       equ1Dx%data(nr,3) = equ1D_in%data(nr)%ppp
       equ1Dx%data(nr,4) = equ1D_in%data(nr)%piq
       equ1Dx%data(nr,5) = equ1D_in%data(nr)%pip
       equ1Dx%data(nr,6) = equ1D_in%data(nr)%pit
    enddo
    CALL DATE_AND_TIME(equ1Dx%created_date, &
                       equ1Dx%created_time, &
                       equ1Dx%created_timezone)

    if(equ1Dx%status.ge.3) then
       equ1Dx%status=3
    else
       equ1Dx%status=2
    endif
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_put_equ1D'
       write(6,*) '---- equ1Dx%s'
       write(6,'(1P5E12.4)') (equ1Dx%rho(nr),nr=1,equ1Dx%nrmax)
       do nd=1,equ1Dx%ndmax
          write(6,*) '---- ',equ1Dx%kid(nd)
          write(6,'(1P5E12.4)') (equ1Dx%data(nr,nd),nr=1,equ1Dx%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_put_equ1D

!-----------------------------------------------------------------------
  subroutine bpsd_get_equ1D(equ1D_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_equ1D_type),intent(inout) :: equ1D_out
    integer,intent(out) :: ierr
    integer :: nr, nd, mode
    real(dp), dimension(6) :: v

    if(bpsd_equ1Dx_init_flag) call bpsd_init_equ1Dx

    if(equ1Dx%status.eq.0) then
       write(6,*) 'XX bpsd_get_equ1D: no space allocated to equ1Dx%data'
       ierr=1
       return
    endif

    if(equ1Dx%status.eq.1) then
       write(6,*) 'XX bpsd_get_equ1D: no data in equ1Dx%data'
       ierr=2
       return
    endif

    if(equ1D_out%nrmax.eq.0) then
       mode=0
       equ1D_out%nrmax = equ1Dx%nrmax
    else
       mode=1
    endif

    CALL bpsd_adjust_array1D(equ1D_out%rho,equ1D_out%nrmax)
    CALL bpsd_adjust_equ1D_data(equ1D_out%data,equ1D_out%nrmax)

    if(mode.eq.0) then
       equ1D_out%time  = equ1Dx%time
       do nr=1,equ1D_out%nrmax
          equ1D_out%rho(nr)       = equ1Dx%rho(nr)
          equ1D_out%data(nr)%psit = equ1Dx%data(nr,1)
          equ1D_out%data(nr)%psip = equ1Dx%data(nr,2)
          equ1D_out%data(nr)%ppp  = equ1Dx%data(nr,3)
          equ1D_out%data(nr)%piq  = equ1Dx%data(nr,4)
          equ1D_out%data(nr)%pip  = equ1Dx%data(nr,5)
          equ1D_out%data(nr)%pit  = equ1Dx%data(nr,6)
       enddo
       ierr=0
       return
    endif

!   A put on a different grid after a get leaves status 3 with the spline
!   of the earlier grid: adjust it here too (a no-op when the size is the same).
    if(equ1Dx%status.eq.2 .or. equ1Dx%status.eq.3) then
       CALL bpsd_adjust_array3D(equ1Dx%spline,4,equ1Dx%nrmax,equ1Dx%ndmax)
       equ1Dx%status=3
    endif

    if(equ1Dx%status.eq.3) then
       CALL bpsd_adjust_array1D(equ1Dx%s,equ1Dx%nrmax)
       do nr=1,equ1Dx%nrmax
          equ1Dx%s(nr)=equ1Dx%rho(nr)**2
       enddo
       do nd=1,6
          call bpsd_spl1D(equ1Dx,nd,ierr)
       enddo
       equ1Dx%status=4
    endif

    do nr=1,equ1D_out%nrmax
       do nd=1,6
          call bpsd_spl1DF(equ1D_out%rho(nr),v(nd),equ1Dx,nd,ierr)
       enddo
       equ1D_out%data(nr)%psit = v(1)
       equ1D_out%data(nr)%psip = v(2)
       equ1D_out%data(nr)%ppp  = v(3)
       equ1D_out%data(nr)%piq  = v(4)
       equ1D_out%data(nr)%pip  = v(5)
       equ1D_out%data(nr)%pit  = v(6)
    enddo
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_equ1D'
       write(6,*) '---- equ1Dx%rho'
       write(6,'(1P5E12.4)') &
                     (equ1D_out%rho(nr),nr=1,equ1D_out%nrmax)
       write(6,*) '---- equ1Dx%psit'
       write(6,'(1P5E12.4)') &
                     (equ1D_out%data(nr)%psit,nr=1,equ1D_out%nrmax)
       write(6,*) '---- equ1Dx%psip'
       write(6,'(1P5E12.4)') &
                     (equ1D_out%data(nr)%psip,nr=1,equ1D_out%nrmax)
       write(6,*) '---- equ1Dx%ppp'
       write(6,'(1P5E12.4)') &
                     (equ1D_out%data(nr)%ppp,nr=1,equ1D_out%nrmax)
       write(6,*) '---- equ1Dx%piq'
       write(6,'(1P5E12.4)') &
                     (equ1D_out%data(nr)%piq,nr=1,equ1D_out%nrmax)
       write(6,*) '---- equ1Dx%pip'
       write(6,'(1P5E12.4)') &
                     (equ1D_out%data(nr)%pip,nr=1,equ1D_out%nrmax)
       write(6,*) '---- equ1Dx%pit'
       write(6,'(1P5E12.4)') &
                     (equ1D_out%data(nr)%pit,nr=1,equ1D_out%nrmax)
    endif
    return
  end subroutine bpsd_get_equ1D

!-----------------------------------------------------------------------
  subroutine bpsd_save_equ1D(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_equ1Dx_init_flag) call bpsd_init_equ1Dx

    if(equ1Dx%status.gt.1) call bpsd_save_data1Dx(fid,equ1Dx,ierr)
    return

  end subroutine bpsd_save_equ1D

!-----------------------------------------------------------------------
  subroutine bpsd_load_equ1D(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data1Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: ns,nr,nd

    if(bpsd_equ1Dx_init_flag) call bpsd_init_equ1Dx

    equ1Dx%dataName=datax%dataName
    equ1Dx%time = datax%time
    equ1Dx%nrmax=datax%nrmax
    equ1Dx%ndmax=datax%ndmax
    CALL bpsd_adjust_karray(equ1Dx%kid,equ1Dx%ndmax)
    CALL bpsd_adjust_karray(equ1Dx%kunit,equ1Dx%ndmax)
    CALL bpsd_adjust_array1D(equ1Dx%rho,equ1Dx%nrmax)
    CALL bpsd_adjust_array2D(equ1Dx%data,equ1Dx%nrmax,equ1Dx%ndmax)

    do nr=1,equ1Dx%nrmax
       equ1Dx%rho(nr) = datax%rho(nr)
       do nd=1,equ1Dx%ndmax
          equ1Dx%data(nr,nd) = datax%data(nr,nd)
       enddo
    enddo

    if(equ1Dx%status.ge.3) then
       equ1Dx%status=3
    else
       equ1Dx%status=2
    endif
    ierr=0
    return

  end subroutine bpsd_load_equ1D

end module bpsd_equ1D
