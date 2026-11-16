! bpsd_metric1D.f90

module bpsd_metric1D

  use bpsd_kinds
  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_put_metric1D,bpsd_get_metric1D, &
         bpsd_save_metric1D,bpsd_load_metric1D
  private

  logical, save :: bpsd_metric1Dx_init_flag = .TRUE.
  type(bpsd_data1Dx_type), save :: metric1Dx

contains

!-----------------------------------------------------------------------
  subroutine bpsd_init_metric1Dx
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    metric1Dx%status=0
    metric1Dx%dataName='metric1D'
    metric1Dx%ndmax=20

    allocate(metric1Dx%kid(metric1Dx%ndmax))
    allocate(metric1Dx%kunit(metric1Dx%ndmax))
    metric1Dx%kid( 1)='metric1D%pvol'
    metric1Dx%kid( 2)='metric1D%psur'
    metric1Dx%kid( 3)='metric1D%dvpsit'
    metric1Dx%kid( 4)='metric1D%dvpsip'
    metric1Dx%kid( 5)='metric1D%aver2'
    metric1Dx%kid( 6)='metric1D%aver2i'
    metric1Dx%kid( 7)='metric1D%aveb2'
    metric1Dx%kid( 8)='metric1D%aveb2i'
    metric1Dx%kid( 9)='metric1D%avegv'
    metric1Dx%kid(10)='metric1D%avegv2'
    metric1Dx%kid(11)='metric1D%avegvr2'
    metric1Dx%kid(12)='metric1D%avegr'
    metric1Dx%kid(13)='metric1D%avegr2'
    metric1Dx%kid(14)='metric1D%avegrr2'
    metric1Dx%kid(15)='metric1D%avegpp2'
    metric1Dx%kid(16)='metric1D%rr'
    metric1Dx%kid(17)='metric1D%rs'
    metric1Dx%kid(18)='metric1D%elip'
    metric1Dx%kid(19)='metric1D%trig'
    metric1Dx%kid(20)='metric1D%aveb'
    metric1Dx%kunit( 1)='m^3 '
    metric1Dx%kunit( 2)='m^2 '
    metric1Dx%kunit( 3)=' '
    metric1Dx%kunit( 4)=' '
    metric1Dx%kunit( 5)='m^2 '
    metric1Dx%kunit( 6)='m^-2 '
    metric1Dx%kunit( 7)='T^2 '
    metric1Dx%kunit( 8)='T^-2'
    metric1Dx%kunit( 9)=' '
    metric1Dx%kunit(10)=' '
    metric1Dx%kunit(11)=' '
    metric1Dx%kunit(12)=' '
    metric1Dx%kunit(13)=' '
    metric1Dx%kunit(14)=' '
    metric1Dx%kunit(15)=' '
    metric1Dx%kunit(16)='m'
    metric1Dx%kunit(17)='m'
    metric1Dx%kunit(18)=' '
    metric1Dx%kunit(19)=' '
    metric1Dx%kunit(20)='T '

    bpsd_metric1Dx_init_flag = .FALSE.

    return
  end subroutine bpsd_init_metric1Dx


!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_metric1D_data(data,n1)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(bpsd_metric1D_data),DIMENSION(:),ALLOCATABLE:: data
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
  END SUBROUTINE bpsd_adjust_metric1D_data

!-----------------------------------------------------------------------
  subroutine bpsd_put_metric1D(metric1D_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_metric1D_type):: metric1D_in
    integer :: ierr
    integer :: nr,nd

    if(bpsd_metric1Dx_init_flag) call bpsd_init_metric1Dx

    metric1Dx%nrmax=metric1D_in%nrmax
    CALL bpsd_adjust_array1D(metric1Dx%rho,metric1Dx%nrmax)
    CALL bpsd_adjust_array2D(metric1Dx%data,metric1Dx%nrmax,metric1Dx%ndmax)

    metric1Dx%time=metric1D_in%time
    do nr=1,metric1Dx%nrmax
       metric1Dx%rho(nr)     = metric1D_in%rho(nr)
       metric1Dx%data(nr, 1) = metric1D_in%data(nr)%pvol
       metric1Dx%data(nr, 2) = metric1D_in%data(nr)%psur
       metric1Dx%data(nr, 3) = metric1D_in%data(nr)%dvpsit
       metric1Dx%data(nr, 4) = metric1D_in%data(nr)%dvpsip
       metric1Dx%data(nr, 5) = metric1D_in%data(nr)%aver2
       metric1Dx%data(nr, 6) = metric1D_in%data(nr)%aver2i
       metric1Dx%data(nr, 7) = metric1D_in%data(nr)%aveb2
       metric1Dx%data(nr, 8) = metric1D_in%data(nr)%aveb2i
       metric1Dx%data(nr, 9) = metric1D_in%data(nr)%avegv
       metric1Dx%data(nr,10) = metric1D_in%data(nr)%avegv2
       metric1Dx%data(nr,11) = metric1D_in%data(nr)%avegvr2
       metric1Dx%data(nr,12) = metric1D_in%data(nr)%avegr
       metric1Dx%data(nr,13) = metric1D_in%data(nr)%avegr2
       metric1Dx%data(nr,14) = metric1D_in%data(nr)%avegrr2
       metric1Dx%data(nr,15) = metric1D_in%data(nr)%avegpp2
       metric1Dx%data(nr,16) = metric1D_in%data(nr)%rr
       metric1Dx%data(nr,17) = metric1D_in%data(nr)%rs
       metric1Dx%data(nr,18) = metric1D_in%data(nr)%elip
       metric1Dx%data(nr,19) = metric1D_in%data(nr)%trig
       metric1Dx%data(nr,20) = metric1D_in%data(nr)%aveb
    enddo
    CALL DATE_AND_TIME(metric1Dx%created_date, &
                       metric1Dx%created_time, &
                       metric1Dx%created_timezone)

    if(metric1Dx%status.ge.3) then
       metric1Dx%status=3
    else
       metric1Dx%status=2
    endif
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_put_metric1D'
       write(6,*) '---- metric1Dx%s'
       write(6,'(1P5E12.4)') &
                     (metric1Dx%rho(nr),nr=1,metric1Dx%nrmax)
       do nd=1,metric1Dx%ndmax
          write(6,*) '---- ',metric1Dx%kid(nd)
          write(6,'(1P5E12.4)') &
                           (metric1Dx%data(nr,nd),nr=1,metric1Dx%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_put_metric1D

!-----------------------------------------------------------------------
  subroutine bpsd_get_metric1D(metric1D_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_metric1D_type),intent(inout) :: metric1D_out
    integer,intent(out) :: ierr
    integer :: nr, nd, mode
    real(dp), dimension(20) :: v

    if(bpsd_metric1Dx_init_flag) call bpsd_init_metric1Dx

    if(metric1Dx%status.eq.0) then
       write(6,*) &
              'XX bpsd_get_metric1D: no space allocated to metric1Dx%data'
       ierr=1
       return
    endif

    if(metric1Dx%status.eq.1) then
       write(6,*) 'XX bpsd_get_metric1D: no data in metric1Dx%data'
       ierr=2
       return
    endif

    if(metric1D_out%nrmax.eq.0) then
       mode=0
       metric1D_out%nrmax = metric1Dx%nrmax
    else
       mode=1
    endif

    CALL bpsd_adjust_array1D(metric1D_out%rho,metric1D_out%nrmax)
    CALL bpsd_adjust_metric1D_data(metric1D_out%data,metric1D_out%nrmax)

    if(mode.eq.0) then
       metric1D_out%time  = metric1Dx%time
       do nr=1,metric1D_out%nrmax
          metric1D_out%rho(nr)           = metric1Dx%rho(nr)
          metric1D_out%data(nr)%pvol     = metric1Dx%data(nr, 1)
          metric1D_out%data(nr)%psur     = metric1Dx%data(nr, 2)
          metric1D_out%data(nr)%dvpsit   = metric1Dx%data(nr, 3)
          metric1D_out%data(nr)%dvpsip   = metric1Dx%data(nr, 4)
          metric1D_out%data(nr)%aver2    = metric1Dx%data(nr, 5)
          metric1D_out%data(nr)%aver2i   = metric1Dx%data(nr, 6)
          metric1D_out%data(nr)%aveb2    = metric1Dx%data(nr, 7)
          metric1D_out%data(nr)%aveb2i   = metric1Dx%data(nr, 8)
          metric1D_out%data(nr)%avegv    = metric1Dx%data(nr, 9)
          metric1D_out%data(nr)%avegv2   = metric1Dx%data(nr,10)
          metric1D_out%data(nr)%avegvr2  = metric1Dx%data(nr,11)
          metric1D_out%data(nr)%avegr    = metric1Dx%data(nr,12)
          metric1D_out%data(nr)%avegr2   = metric1Dx%data(nr,13)
          metric1D_out%data(nr)%avegrr2  = metric1Dx%data(nr,14)
          metric1D_out%data(nr)%avegpp2  = metric1Dx%data(nr,15)
          metric1D_out%data(nr)%rr       = metric1Dx%data(nr,16)
          metric1D_out%data(nr)%rs       = metric1Dx%data(nr,17)
          metric1D_out%data(nr)%elip     = metric1Dx%data(nr,18)
          metric1D_out%data(nr)%trig     = metric1Dx%data(nr,19)
          metric1D_out%data(nr)%aveb     = metric1Dx%data(nr,20)
       enddo
       ierr=0
       return
    endif

!   A put on a different grid after a get leaves status 3 with the spline
!   of the earlier grid: adjust it here too (a no-op when the size is the same).
    if(metric1Dx%status.eq.2 .or. metric1Dx%status.eq.3) then
       CALL bpsd_adjust_array3D(metric1Dx%spline,4,metric1Dx%nrmax, &
                                                   metric1Dx%ndmax)
       metric1Dx%status=3
    endif

    if(metric1Dx%status.eq.3) then
       CALL bpsd_adjust_array1D(metric1Dx%s,metric1Dx%nrmax)
       do nr=1,metric1Dx%nrmax
          metric1Dx%s(nr)=metric1Dx%rho(nr)**2
       enddo
       do nd=1,metric1Dx%ndmax
          call bpsd_spl1D(metric1Dx,nd,ierr)
       enddo
       metric1Dx%status=4
    endif

    do nr=1,metric1D_out%nrmax
       do nd=1,metric1Dx%ndmax
          call bpsd_spl1DF(metric1D_out%rho(nr),v(nd),metric1Dx,nd,ierr)
       enddo
       metric1D_out%data(nr)%pvol     = v( 1)
       metric1D_out%data(nr)%psur     = v( 2)
       metric1D_out%data(nr)%dvpsit   = v( 3)
       metric1D_out%data(nr)%dvpsip   = v( 4)
       metric1D_out%data(nr)%aver2    = v( 5)
       metric1D_out%data(nr)%aver2i   = v( 6)
       metric1D_out%data(nr)%aveb2    = v( 7)
       metric1D_out%data(nr)%aveb2i   = v( 8)
       metric1D_out%data(nr)%avegv    = v( 9)
       metric1D_out%data(nr)%avegv2   = v(10)
       metric1D_out%data(nr)%avegvr2  = v(11)
       metric1D_out%data(nr)%avegr    = v(12)
       metric1D_out%data(nr)%avegr2   = v(13)
       metric1D_out%data(nr)%avegrr2  = v(14)
       metric1D_out%data(nr)%avegpp2  = v(15)
       metric1D_out%data(nr)%rr       = v(16)
       metric1D_out%data(nr)%rs       = v(17)
       metric1D_out%data(nr)%elip     = v(18)
       metric1D_out%data(nr)%trig     = v(19)
       metric1D_out%data(nr)%aveb     = v(20)
    enddo
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_metric1D'
       write(6,*) '---- metric1Dx%rho'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%rho(nr),nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%pvol'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%pvol,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%psur'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%psur,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%dvpssit'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%dvpsit,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%dvpssit'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%dvpsit,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%aver2'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%aver2,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%aver2i'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%aver2i,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%aveb2'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%aveb2,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%aveb2i'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%aveb2i,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%avegv'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%avegv,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%avegv2'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%avegv2,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%avegvr2'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%avegvr2,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%avegr'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%avegr,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%avegr2'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%avegr2,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%avegrr2'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%avegrr2,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%avegpp2'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%avegpp2,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%rr'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%rr,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%rs'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%rs,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%elip'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%elip,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%trig'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%trig,nr=1,metric1D_out%nrmax)
       write(6,*) '---- metric1Dx%aveb'
       write(6,'(1P5E12.4)') &
                     (metric1D_out%data(nr)%aveb,nr=1,metric1D_out%nrmax)
    endif
    return
  end subroutine bpsd_get_metric1D

!-----------------------------------------------------------------------
  subroutine bpsd_save_metric1D(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_metric1Dx_init_flag) call bpsd_init_metric1Dx

    if(metric1Dx%status.gt.1) call bpsd_save_data1Dx(fid,metric1Dx,ierr)
    return

  end subroutine bpsd_save_metric1D

!-----------------------------------------------------------------------
  subroutine bpsd_load_metric1D(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data1Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: ns,nr,nd

    if(bpsd_metric1Dx_init_flag) call bpsd_init_metric1Dx

    metric1Dx%dataName=datax%dataName
    metric1Dx%time = datax%time
    metric1Dx%nrmax=datax%nrmax
    metric1Dx%ndmax=datax%ndmax
    CALL bpsd_adjust_karray(metric1Dx%kid,metric1Dx%ndmax)
    CALL bpsd_adjust_karray(metric1Dx%kunit,metric1Dx%ndmax)
    CALL bpsd_adjust_array1D(metric1Dx%rho,metric1Dx%nrmax)
    CALL bpsd_adjust_array2D(metric1Dx%data,metric1Dx%nrmax, &
                                            metric1Dx%ndmax)

    do nr=1,metric1Dx%nrmax
       metric1Dx%rho(nr) = datax%rho(nr)
       do nd=1,metric1Dx%ndmax
          metric1Dx%data(nr,nd) = datax%data(nr,nd)
       enddo
    enddo

    if(metric1Dx%status.ge.3) then
       metric1Dx%status=3
    else
       metric1Dx%status=2
    endif
    ierr=0
    return

  end subroutine bpsd_load_metric1D

end module bpsd_metric1D

