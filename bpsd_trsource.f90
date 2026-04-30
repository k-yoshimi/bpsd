! bpsd_trsource.f90

module bpsd_trsource

  use bpsd_kinds
  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_put_trsource,bpsd_get_trsource, &
         bpsd_save_trsource,bpsd_load_trsource
  private

  ! Number of scalar fields packed per species in the flat trsourcex
  ! buffer. Single source of truth for the put/get/setup_kdata layout.
  integer(ikind), parameter :: nfields = 10
  ! nip, nim, ncx, Pec, Plh, Pic, Pbr, Pcy, Plr, Poh

  logical, save :: bpsd_trsourcex_init_flag = .TRUE.
  type(bpsd_data1Dx_type), save :: trsourcex

contains

!-----------------------------------------------------------------------
  subroutine bpsd_init_trsourcex
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    trsourcex%status=0
    trsourcex%dataName='trsource'
    trsourcex%ndmax=0

    bpsd_trsourcex_init_flag = .FALSE.

    return
  end subroutine bpsd_init_trsourcex

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_trsource_data(data,n1,n2)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(bpsd_trsource_data),DIMENSION(:,:),ALLOCATABLE:: data
    INTEGER(ikind),INTENT(IN):: n1,n2

    IF(ALLOCATED(data)) THEN
       IF(n1.LE.0.or.n2.LE.0) THEN
          DEALLOCATE(data)
       ELSE IF(n1.NE.SIZE(data,1).OR.n2.NE.SIZE(data,2)) THEN
          DEALLOCATE(data)
          ALLOCATE(data(n1,n2))
       END IF
    ELSE
       IF(n1.GT.0.AND.n2.GT.0)  THEN
          ALLOCATE(data(n1,n2))
       ENDIF
    ENDIF
  END SUBROUTINE bpsd_adjust_trsource_data

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_setup_trsource_kdata
!-----------------------------------------------------------------------
    IMPLICIT NONE
    integer(ikind):: nd

    do nd=0,trsourcex%ndmax-nfields,nfields
       trsourcex%kid(nd+1)='trsource%nip'
       trsourcex%kid(nd+2)='trsource%nim'
       trsourcex%kid(nd+3)='trsource%ncx'
       trsourcex%kid(nd+4)='trsource%Pec'
       trsourcex%kid(nd+5)='trsource%Plh'
       trsourcex%kid(nd+6)='trsource%Pic'
       trsourcex%kid(nd+7)='trsource%Pbr'
       trsourcex%kid(nd+8)='trsource%Pcy'
       trsourcex%kid(nd+9)='trsource%Plr'
       trsourcex%kid(nd+10)='trsource%Poh'
       trsourcex%kunit(nd+1)='1/(m^3 s)'
       trsourcex%kunit(nd+2)='1/(m^3 s)'
       trsourcex%kunit(nd+3)='1/(m^3 s)'
       trsourcex%kunit(nd+4)='W/m^3'
       trsourcex%kunit(nd+5)='W/m^3'
       trsourcex%kunit(nd+6)='W/m^3'
       trsourcex%kunit(nd+7)='W/m^3'
       trsourcex%kunit(nd+8)='W/m^3'
       trsourcex%kunit(nd+9)='W/m^3'
       trsourcex%kunit(nd+10)='W/m^3'
    enddo
    RETURN
  END SUBROUTINE bpsd_setup_trsource_kdata

!-----------------------------------------------------------------------
  subroutine bpsd_put_trsource(trsource_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_trsource_type),intent(in):: trsource_in
    integer,intent(out) :: ierr
    integer :: ns,nr,nd

    if(bpsd_trsourcex_init_flag) call bpsd_init_trsourcex

    trsourcex%nrmax=trsource_in%nrmax
    trsourcex%ndmax=trsource_in%nsmax*nfields
    CALL bpsd_adjust_karray(trsourcex%kid,trsourcex%ndmax)
    CALL bpsd_adjust_karray(trsourcex%kunit,trsourcex%ndmax)
    CALL bpsd_adjust_array1D(trsourcex%rho,trsourcex%nrmax)
    CALL bpsd_adjust_array2D(trsourcex%data,trsourcex%nrmax,trsourcex%ndmax)

    CALL bpsd_setup_trsource_kdata

    trsourcex%time = trsource_in%time
    do nr=1,trsource_in%nrmax
       trsourcex%rho(nr) = trsource_in%rho(nr)
       do ns=1,trsource_in%nsmax
          nd=nfields*(ns-1)
          trsourcex%data(nr,nd+1) = trsource_in%data(nr,ns)%nip
          trsourcex%data(nr,nd+2) = trsource_in%data(nr,ns)%nim
          trsourcex%data(nr,nd+3) = trsource_in%data(nr,ns)%ncx
          trsourcex%data(nr,nd+4) = trsource_in%data(nr,ns)%Pec
          trsourcex%data(nr,nd+5) = trsource_in%data(nr,ns)%Plh
          trsourcex%data(nr,nd+6) = trsource_in%data(nr,ns)%Pic
          trsourcex%data(nr,nd+7) = trsource_in%data(nr,ns)%Pbr
          trsourcex%data(nr,nd+8) = trsource_in%data(nr,ns)%Pcy
          trsourcex%data(nr,nd+9) = trsource_in%data(nr,ns)%Plr
          trsourcex%data(nr,nd+10)= trsource_in%data(nr,ns)%Poh
       enddo
    enddo
    CALL DATE_AND_TIME(trsourcex%created_date, &
                       trsourcex%created_time, &
                       trsourcex%created_timezone)

    if(trsourcex%status.ge.3) then
       trsourcex%status=3
    else
       trsourcex%status=2
    endif

    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_put_trsource'
       write(6,*) '---- trsourcex%rho'
       write(6,'(1P5E12.4)') &
                     (trsourcex%rho(nr),nr=1,trsourcex%nrmax)
       do nd=1,trsourcex%ndmax
          write(6,*) '---- ',trsourcex%kid(nd)
          write(6,'(1P5E12.4)') &
                           (trsourcex%data(nr,nd),nr=1,trsourcex%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_put_trsource

!-----------------------------------------------------------------------
  subroutine bpsd_get_trsource(trsource_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_trsource_type),intent(out) :: trsource_out
    integer,intent(out) :: ierr
    integer :: nr, nd, ns, mode
    real(dp), dimension(:), ALLOCATABLE :: v

    ! INTENT(OUT) on a derived type with allocatable components leaves
    ! scalar fields UNDEFINED on entry per Fortran 2003+. The body
    ! below tests trsource_out%nrmax to choose mode; without this reset
    ! that test reads heap garbage. Mirrors the plasmaf fix in
    ! bpsd_get_plasmaf.
    trsource_out%nrmax = 0
    trsource_out%nsmax = 0
    trsource_out%time  = 0.0_dp

    if(bpsd_trsourcex_init_flag) call bpsd_init_trsourcex

    if(trsourcex%status.eq.0) then
       write(6,*) 'XX bpsd_get_trsource: no space allocated to trsourcex%data'
       ierr=1
       return
    endif

    if(trsourcex%status.eq.1) then
       write(6,*) 'XX bpsd_get_trsource: no data in trsourcex%data'
       ierr=2
       return
    endif

    if(trsource_out%nrmax.eq.0) then
       mode=0
       trsource_out%nrmax = trsourcex%nrmax
    else
       mode=1
    endif
    trsource_out%nsmax = trsourcex%ndmax/nfields

    CALL bpsd_adjust_array1D(trsource_out%rho,trsource_out%nrmax)
    CALL bpsd_adjust_trsource_data(trsource_out%data,trsource_out%nrmax, &
                                                     trsource_out%nsmax)

    if(mode.eq.0) then
       trsource_out%time  = trsourcex%time
       do nr=1,trsourcex%nrmax
          trsource_out%rho(nr)=trsourcex%rho(nr)
          do ns=1,trsource_out%nsmax
             nd=nfields*(ns-1)
             trsource_out%data(nr,ns)%nip =trsourcex%data(nr,nd+1)
             trsource_out%data(nr,ns)%nim =trsourcex%data(nr,nd+2)
             trsource_out%data(nr,ns)%ncx =trsourcex%data(nr,nd+3)
             trsource_out%data(nr,ns)%Pec =trsourcex%data(nr,nd+4)
             trsource_out%data(nr,ns)%Plh =trsourcex%data(nr,nd+5)
             trsource_out%data(nr,ns)%Pic =trsourcex%data(nr,nd+6)
             trsource_out%data(nr,ns)%Pbr =trsourcex%data(nr,nd+7)
             trsource_out%data(nr,ns)%Pcy =trsourcex%data(nr,nd+8)
             trsource_out%data(nr,ns)%Plr =trsourcex%data(nr,nd+9)
             trsource_out%data(nr,ns)%Poh =trsourcex%data(nr,nd+10)
          enddo
       enddo
       ierr=0
       return
    endif

    if(trsourcex%status.eq.2) then
       ! spline first-dim is the cubic-spline coefficient count from
       ! spl1D's U(4, NXMAX) interface; unrelated to nfields.
       CALL bpsd_adjust_array3D(trsourcex%spline,4,trsourcex%nrmax, &
                                                    trsourcex%ndmax)
       trsourcex%status=3
    endif

    if(trsourcex%status.eq.3) then
       CALL bpsd_adjust_array1D(trsourcex%s,trsourcex%nrmax)
       do nr=1,trsourcex%nrmax
          trsourcex%s(nr)=trsourcex%rho(nr)**2
       enddo
       do nd=1,trsourcex%ndmax
          call bpsd_spl1D(trsourcex,nd,ierr)
       enddo
       trsourcex%status=4
    endif

    allocate(v(trsourcex%ndmax))

    do nr=1,trsource_out%nrmax
       do nd=1,trsourcex%ndmax
          call bpsd_spl1DF(trsource_out%rho(nr),v(nd),trsourcex,nd,ierr)
       enddo
       do ns=1,trsource_out%nsmax
          nd=nfields*(ns-1)
          trsource_out%data(nr,ns)%nip  = v(nd+1)
          trsource_out%data(nr,ns)%nim  = v(nd+2)
          trsource_out%data(nr,ns)%ncx  = v(nd+3)
          trsource_out%data(nr,ns)%Pec  = v(nd+4)
          trsource_out%data(nr,ns)%Plh  = v(nd+5)
          trsource_out%data(nr,ns)%Pic  = v(nd+6)
          trsource_out%data(nr,ns)%Pbr  = v(nd+7)
          trsource_out%data(nr,ns)%Pcy  = v(nd+8)
          trsource_out%data(nr,ns)%Plr  = v(nd+9)
          trsource_out%data(nr,ns)%Poh  = v(nd+10)
       enddo
    enddo
    deallocate(v)
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_trsource'
       write(6,*) '---- trsourcex%rho'
       write(6,'(1P5E12.4)') &
                     (trsource_out%rho(nr),nr=1,trsource_out%nrmax)
       do ns=1,trsource_out%nsmax
          write(6,*) '---- trsourcex%nip(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%nip,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%nim(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%nim,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%ncx(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%ncx,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%Pec(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%Pec,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%Plh(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%Plh,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%Pic(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%Pic,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%Pbr(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%Pbr,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%Pcy(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%Pcy,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%Plr(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%Plr,nr=1,trsource_out%nrmax)
          write(6,*) '---- trsourcex%Poh(',ns,')'
          write(6,'(1P5E12.4)') &
                     (trsource_out%data(nr,ns)%Poh,nr=1,trsource_out%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_get_trsource

!-----------------------------------------------------------------------
  subroutine bpsd_save_trsource(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_trsourcex_init_flag) call bpsd_init_trsourcex

    if(trsourcex%status.gt.1) call bpsd_save_data1Dx(fid,trsourcex,ierr)
    return

  end subroutine bpsd_save_trsource

!-----------------------------------------------------------------------
  subroutine bpsd_load_trsource(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data1Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: ns,nr,nd

    if(bpsd_trsourcex_init_flag) call bpsd_init_trsourcex

    trsourcex%dataName=datax%dataName
    trsourcex%time = datax%time
    trsourcex%nrmax=datax%nrmax
    trsourcex%ndmax=datax%ndmax
    CALL bpsd_adjust_karray(trsourcex%kid,trsourcex%ndmax)
    CALL bpsd_adjust_karray(trsourcex%kunit,trsourcex%ndmax)
    CALL bpsd_adjust_array1D(trsourcex%rho,trsourcex%nrmax)
    CALL bpsd_adjust_array2D(trsourcex%data,trsourcex%nrmax, &
                                            trsourcex%ndmax)

    do nr=1,trsourcex%nrmax
       trsourcex%rho(nr) = datax%rho(nr)
       do nd=1,trsourcex%ndmax
          trsourcex%data(nr,nd) = datax%data(nr,nd)
       enddo
    enddo
    do nd=1,trsourcex%ndmax
       trsourcex%kid(nd)=datax%kid(nd)
       trsourcex%kunit(nd)=datax%kunit(nd)
    enddo
    trsourcex%created_date = datax%created_date
    trsourcex%created_time = datax%created_time
    trsourcex%created_timezone = datax%created_timezone

    if(trsourcex%status.ge.3) then
       trsourcex%status=3
    else
       trsourcex%status=2
    endif

    ierr=0
    return

  end subroutine bpsd_load_trsource

end module bpsd_trsource
