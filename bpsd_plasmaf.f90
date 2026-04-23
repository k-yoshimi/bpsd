! bpsd_plasma.f90

module bpsd_plasmaf

  use bpsd_kinds
  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal

  PRIVATE
  PUBLIC bpsd_put_plasmaf
  PUBLIC bpsd_get_plasmaf
  PUBLIC bpsd_save_plasmaf
  PUBLIC bpsd_load_plasmaf

  logical, save :: bpsd_plasmafx_init_flag = .TRUE.
  type(bpsd_data1Dx_type), save :: plasmafx

contains

!-----------------------------------------------------------------------
  subroutine bpsd_init_plasmafx
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    plasmafx%status=0
    plasmafx%dataName='plasmaf'
    plasmafx%ndmax=0

    bpsd_plasmafx_init_flag = .FALSE.

    return
  end subroutine bpsd_init_plasmafx

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_plasmaf_data(data,n1,n2)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(bpsd_plasmaf_data),DIMENSION(:,:),ALLOCATABLE:: data
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
  END SUBROUTINE bpsd_adjust_plasmaf_data

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_setup_plasmaf_kdata
!-----------------------------------------------------------------------
    IMPLICIT NONE
    integer(ikind):: nd

    do nd=0,plasmafx%ndmax-2,12
       plasmafx%kid(nd+1)='plasmaf%density'
       plasmafx%kid(nd+2)='plasmaf%temperature'
       plasmafx%kid(nd+3)='plasmaf%temperature_para'
       plasmafx%kid(nd+4)='plasmaf%temperature_perp'
       plasmafx%kid(nd+5)='plasmaf%velocity_tor'
       plasmafx%kid(nd+6)='plasmaf%velocity_pol'
       plasmafx%kid(nd+7)='plasmaf%velocity_para'
       plasmafx%kid(nd+8)='plasmaf%velocity_perp'
       plasmafx%kid(nd+9)='plasmaf%zave'
       plasmafx%kid(nd+10)='plasmaf%z2ave'
       plasmafx%kid(nd+11)='plasmaf%density_fastion'
       plasmafx%kid(nd+12)='plasmaf%energy_fastion'
       plasmafx%kunit(nd+1)='10^20/m^3'
       plasmafx%kunit(nd+2)='eV'
       plasmafx%kunit(nd+3)='eV'
       plasmafx%kunit(nd+4)='eV'
       plasmafx%kunit(nd+5)='m/s'
       plasmafx%kunit(nd+6)='m/s'
       plasmafx%kunit(nd+7)='m/s'
       plasmafx%kunit(nd+8)='m/s'
       plasmafx%kunit(nd+9)=' '
       plasmafx%kunit(nd+10)=' '
       plasmafx%kunit(nd+11)='10^20/m^3'
       plasmafx%kunit(nd+12)='eV'
    enddo
    plasmafx%kid(plasmafx%ndmax)='plasmaf%qinv'
    plasmafx%kunit(plasmafx%ndmax)=' '
    RETURN
  END SUBROUTINE bpsd_setup_plasmaf_kdata

!-----------------------------------------------------------------------
  subroutine bpsd_put_plasmaf(plasmaf_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_plasmaf_type),intent(in):: plasmaf_in
    integer,intent(out) :: ierr
    integer :: ns,nr,nd

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    plasmafx%nrmax=plasmaf_in%nrmax
    plasmafx%ndmax=plasmaf_in%nsmax*12+1
    CALL bpsd_adjust_karray(plasmafx%kid,plasmafx%ndmax)
    CALL bpsd_adjust_karray(plasmafx%kunit,plasmafx%ndmax)
    CALL bpsd_adjust_array1D(plasmafx%rho,plasmafx%nrmax)
    CALL bpsd_adjust_array2D(plasmafx%data,plasmafx%nrmax,plasmafx%ndmax)

    CALL bpsd_setup_plasmaf_kdata

    plasmafx%time = plasmaf_in%time
    do nr=1,plasmaf_in%nrmax
       plasmafx%rho(nr) = plasmaf_in%rho(nr)
       do ns=1,plasmaf_in%nsmax
          nd=12*(ns-1)
          plasmafx%data(nr,nd+1) = plasmaf_in%data(nr,ns)%density
          plasmafx%data(nr,nd+2) = plasmaf_in%data(nr,ns)%temperature
          plasmafx%data(nr,nd+3) = plasmaf_in%data(nr,ns)%temperature_para
          plasmafx%data(nr,nd+4) = plasmaf_in%data(nr,ns)%temperature_perp
          plasmafx%data(nr,nd+5) = plasmaf_in%data(nr,ns)%velocity_tor
          plasmafx%data(nr,nd+6) = plasmaf_in%data(nr,ns)%velocity_pol
          plasmafx%data(nr,nd+7) = plasmaf_in%data(nr,ns)%velocity_para
          plasmafx%data(nr,nd+8) = plasmaf_in%data(nr,ns)%velocity_perp
          plasmafx%data(nr,nd+9) = plasmaf_in%data(nr,ns)%zave
          plasmafx%data(nr,nd+10)= plasmaf_in%data(nr,ns)%z2ave
          plasmafx%data(nr,nd+11)= plasmaf_in%data(nr,ns)%density_fastion
          plasmafx%data(nr,nd+12)= plasmaf_in%data(nr,ns)%energy_fastion
       enddo
       plasmafx%data(nr,plasmafx%ndmax) = plasmaf_in%qinv(nr)
    enddo
    CALL DATE_AND_TIME(plasmafx%created_date, &
                       plasmafx%created_time, &
                       plasmafx%created_timezone)

    if(plasmafx%status.ge.3) then
       plasmafx%status=3
    else
       plasmafx%status=2
    endif

    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_put_plasmaf'
       write(6,*) '---- plasmafx%rho'
       write(6,'(1P5E12.4)') &
                     (plasmafx%rho(nr),nr=1,plasmafx%nrmax)
       do nd=1,plasmafx%ndmax
          write(6,*) '---- ',plasmafx%kid(nd)
          write(6,'(1P5E12.4)') &
                           (plasmafx%data(nr,nd),nr=1,plasmafx%nrmax)
       enddo
    endif
    return
  end subroutine bpsd_put_plasmaf

!-----------------------------------------------------------------------
  subroutine bpsd_get_plasmaf(plasmaf_out,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_plasmaf_type),intent(out) :: plasmaf_out
    integer,intent(out) :: ierr
    integer :: nr, nd, ns, mode
    real(dp), dimension(:), ALLOCATABLE :: v

    ! Defensive zero-init of the scalar fields. INTENT(OUT) on a
    ! derived type with allocatable components leaves scalar fields
    ! UNDEFINED on entry per Fortran 2003+; the body below tests
    ! `plasmaf_out%nrmax.eq.0` at line 186 to choose mode 0 vs 1.
    ! Without this reset that test reads heap garbage and dispatches
    ! to mode=1 with stale nrmax → bpsd_adjust_array1D leaves rho's
    ! descriptor in a corrupt state → "Index '1' below lower bound of
    ! <garbage>" SIGABRT in callers (libtotapi.so via tr_loop ->
    ! tr_bpsd_get). Confirmed via file-based dump showing
    ! `ALLOC(rho)=F` AND `out%nrmax=<stale>` at entry. Tracked in
    ! task-private issue #148.
    plasmaf_out%nrmax = 0
    plasmaf_out%nsmax = 0
    plasmaf_out%time  = 0.0_dp

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    if(plasmafx%status.eq.0) then
       write(6,*) 'XX bpsd_get_plasmaf: no space allocated to plasmafx%data'
       ierr=1
       return
    endif

    if(plasmafx%status.eq.1) then
       write(6,*) 'XX bpsd_get_plasmaf: no data in plasmafx%data'
       ierr=2
       return
    endif

    if(plasmaf_out%nrmax.eq.0) then
       mode=0
       plasmaf_out%nrmax = plasmafx%nrmax
    else
       mode=1
    endif
    plasmaf_out%nsmax = (plasmafx%ndmax-1)/12

    CALL bpsd_adjust_array1D(plasmaf_out%rho,plasmaf_out%nrmax)
    CALL bpsd_adjust_array1D(plasmaf_out%qinv,plasmaf_out%nrmax)
    CALL bpsd_adjust_plasmaf_data(plasmaf_out%data,plasmaf_out%nrmax, &
                                                   plasmaf_out%nsmax)

!    if(mode.eq.0) then
       plasmaf_out%time  = plasmafx%time
       do nr=1,plasmafx%nrmax
          plasmaf_out%rho(nr)=plasmafx%rho(nr)
          do ns=1,plasmaf_out%nsmax
             nd=12*(ns-1)
             plasmaf_out%data(nr,ns)%density         =plasmafx%data(nr,nd+1)
             plasmaf_out%data(nr,ns)%temperature     =plasmafx%data(nr,nd+2)
             plasmaf_out%data(nr,ns)%temperature_para=plasmafx%data(nr,nd+3)
             plasmaf_out%data(nr,ns)%temperature_perp=plasmafx%data(nr,nd+4)
             plasmaf_out%data(nr,ns)%velocity_tor    =plasmafx%data(nr,nd+5)
             plasmaf_out%data(nr,ns)%velocity_pol    =plasmafx%data(nr,nd+6)
             plasmaf_out%data(nr,ns)%velocity_para   =plasmafx%data(nr,nd+7)
             plasmaf_out%data(nr,ns)%velocity_perp   =plasmafx%data(nr,nd+8)
             plasmaf_out%data(nr,ns)%zave            =plasmafx%data(nr,nd+9)
             plasmaf_out%data(nr,ns)%z2ave           =plasmafx%data(nr,nd+10)
             plasmaf_out%data(nr,ns)%density_fastion =plasmafx%data(nr,nd+11)
             plasmaf_out%data(nr,ns)%energy_fastion  =plasmafx%data(nr,nd+12)
          enddo
          plasmaf_out%qinv(nr)=plasmafx%data(nr,plasmafx%ndmax)
       enddo
       ierr=0
       return
!    endif

    if(plasmafx%status.eq.2) then
       CALL bpsd_adjust_array3D(plasmafx%spline,4,plasmafx%nrmax, &
                                                  plasmafx%ndmax)
       plasmafx%status=3
    endif

    if(plasmafx%status.eq.3) then
       CALL bpsd_adjust_array1D(plasmafx%s,plasmafx%nrmax)
       do nr=1,plasmafx%nrmax
          plasmafx%s(nr)=plasmafx%rho(nr)**2
       enddo
       do nd=1,plasmafx%ndmax
          call bpsd_spl1D(plasmafx,nd,ierr)
       enddo
       plasmafx%status=4
    endif

    allocate(v(plasmafx%ndmax))

    do nr=1,plasmaf_out%nrmax
       do nd=1,plasmafx%ndmax
          call bpsd_spl1DF(plasmaf_out%rho(nr),v(nd),plasmafx,nd,ierr)
       enddo
       do ns=1,plasmaf_out%nsmax
          nd=12*(ns-1)
          plasmaf_out%data(nr,ns)%density          = v(nd+1)
          plasmaf_out%data(nr,ns)%temperature      = v(nd+2)
          plasmaf_out%data(nr,ns)%temperature_para = v(nd+3)
          plasmaf_out%data(nr,ns)%temperature_perp = v(nd+4)
          plasmaf_out%data(nr,ns)%velocity_tor     = v(nd+5)
          plasmaf_out%data(nr,ns)%velocity_pol     = v(nd+6)
          plasmaf_out%data(nr,ns)%velocity_para    = v(nd+7)
          plasmaf_out%data(nr,ns)%velocity_perp    = v(nd+8)
          plasmaf_out%data(nr,ns)%zave             = v(nd+9)
          plasmaf_out%data(nr,ns)%z2ave            = v(nd+10)
          plasmaf_out%data(nr,ns)%density_fastion  = v(nd+11)
          plasmaf_out%data(nr,ns)%energy_fastion   = v(nd+12)
       enddo
       plasmaf_out%qinv(nr)  = v(plasmafx%ndmax)
    enddo
    deallocate(v)
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_plasmaf'
       write(6,*) '---- plasmafx%rho'
       write(6,'(1P5E12.4)') &
                     (plasmaf_out%rho(nr), &
                     nr=1,plasmaf_out%nrmax)
       do ns=1,plasmaf_out%nsmax
          write(6,*) '---- plasmafx%density(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%density, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%temperature(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%temperature, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%temperature_para(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%temperature_para, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%temperature_perp(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%temperature_perp, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%velocity_tor(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%velocity_tor, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%velocity_pol(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%velocity_pol, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%velocity_para(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%velocity_para, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%velocity_perp(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%velocity_perp, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%zave(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%zave, &
                     nr=1,plasmaf_out%nrmax)
          write(6,*) '---- plasmafx%z2ave(',ns,')'
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%z2ave, &
                     nr=1,plasmaf_out%nrmax)
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%density_fastion, &
                     nr=1,plasmaf_out%nrmax)
          write(6,'(1P5E12.4)') &
                     (plasmaf_out%data(nr,ns)%energy_fastion, &
                     nr=1,plasmaf_out%nrmax)
       enddo
       write(6,*) '---- plasmafx%qinv'
       write(6,'(1P5E12.4)') &
                     (plasmaf_out%qinv(nr),nr=1,plasmaf_out%nrmax)
    endif
    return
  end subroutine bpsd_get_plasmaf

!-----------------------------------------------------------------------
  subroutine bpsd_save_plasmaf(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    if(plasmafx%status.gt.1) call bpsd_save_data1Dx(fid,plasmafx,ierr)
    return

  end subroutine bpsd_save_plasmaf

!-----------------------------------------------------------------------
  subroutine bpsd_load_plasmaf(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data1Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: ns,nr,nd

    if(bpsd_plasmafx_init_flag) call bpsd_init_plasmafx

    plasmafx%dataName=datax%dataName
    plasmafx%time = datax%time
    plasmafx%nrmax=datax%nrmax
    plasmafx%ndmax=datax%ndmax
    CALL bpsd_adjust_karray(plasmafx%kid,plasmafx%ndmax)
    CALL bpsd_adjust_karray(plasmafx%kunit,plasmafx%ndmax)
    CALL bpsd_adjust_array1D(plasmafx%rho,plasmafx%nrmax)
    CALL bpsd_adjust_array2D(plasmafx%data,plasmafx%nrmax, &
                                           plasmafx%ndmax)

    do nr=1,plasmafx%nrmax
       plasmafx%rho(nr) = datax%rho(nr)
       do nd=1,plasmafx%ndmax
          plasmafx%data(nr,nd) = datax%data(nr,nd)
       enddo
    enddo
    do nd=1,plasmafx%ndmax
       plasmafx%kid(nd)=datax%kid(nd)
       plasmafx%kunit(nd)=datax%kunit(nd)
    enddo
    plasmafx%created_date = datax%created_date
    plasmafx%created_time = datax%created_time
    plasmafx%created_timezone = datax%created_timezone

    if(plasmafx%status.ge.3) then
       plasmafx%status=3
    else
       plasmafx%status=2
    endif

    ierr=0
    return

  end subroutine bpsd_load_plasmaf

end module bpsd_plasmaf
