! bpsd_device.f90

module bpsd_device

  use bpsd_flags
  use bpsd_types
  use bpsd_types_internal
  public bpsd_put_device,bpsd_get_device, &
         bpsd_save_device,bpsd_load_device
  private

  logical, save :: bpsd_devicex_init_flag = .TRUE.
  type(bpsd_data0Dx_type), save :: devicex

contains

!-----------------------------------------------------------------------
  subroutine bpsd_devicex_init
!-----------------------------------------------------------------------
    use bpsd_subs
    implicit none

    devicex%status=0
    devicex%dataName='device'
    devicex%ndmax=8
    allocate(devicex%kid(8))
    allocate(devicex%kunit(8))
    allocate(devicex%data(8))
    devicex%kid(1)='device%rr'
    devicex%kid(2)='device%zz'
    devicex%kid(3)='device%ra'
    devicex%kid(4)='device%rb'
    devicex%kid(5)='device%bb'
    devicex%kid(6)='device%ip'
    devicex%kid(7)='device%elip'
    devicex%kid(8)='device%trig'
    devicex%kunit(1)='m'
    devicex%kunit(2)='m'
    devicex%kunit(3)='m'
    devicex%kunit(4)='m'
    devicex%kunit(5)='T'
    devicex%kunit(6)='A'
    devicex%kunit(7)=' '
    devicex%kunit(8)=' '

    bpsd_devicex_init_flag = .FALSE.

    return
  end subroutine bpsd_devicex_init

!-----------------------------------------------------------------------
  subroutine bpsd_put_device(device_in,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_device_type) :: device_in
    integer :: ierr, nd

    if(bpsd_devicex_init_flag) call bpsd_devicex_init

    devicex%dataName = 'device'
    devicex%time = 0.D0
    devicex%data(1) = device_in%rr
    devicex%data(2) = device_in%zz
    devicex%data(3) = device_in%ra
    devicex%data(4) = device_in%rb
    devicex%data(5) = device_in%bb
    devicex%data(6) = device_in%ip
    devicex%data(7) = device_in%elip
    devicex%data(8) = device_in%trig
    CALL DATE_AND_TIME(devicex%created_date, &
                       devicex%created_time, &
                       devicex%created_timezone)
    devicex%status = 2
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_put_device'
       do nd=1,devicex%ndmax
          write(6,'(A32,1PE12.4,A)') devicex%kid(nd),devicex%data(nd),TRIM(devicex%kunit(nd))
       enddo
    endif
    return
  end subroutine bpsd_put_device

!-----------------------------------------------------------------------
  subroutine bpsd_get_device(device_out,ierr)
!-----------------------------------------------------------------------c
    use bpsd_subs
    implicit none
    type(bpsd_device_type) :: device_out
    integer :: ierr, nd

    if(bpsd_devicex_init_flag) call bpsd_devicex_init

    ! status meanings:
    !   0 = devicex_init done; data(8) allocated but never put
    !   2 = bpsd_put_device has populated data(8)
    ! Reject status<2 (covers the never-allocated and the
    ! allocated-but-never-put cases). The earlier check (status==1)
    ! never fired -- devicex_init sets status=0, so a fresh init
    ! followed by a get would silently read uninitialised data(1..8)
    ! (ALLOCATE without zero-init) and return ierr=0. Match the lt-2
    ! pattern used in bpsd_get_species_kunit (bpsd_species.f90:240).
    if(devicex%status.lt.2) then
       write(6,*) 'XX bpsd_get_device: no data in device (status<2)'
       ierr=2
       return
    endif

    device_out%rr    = devicex%data(1)
    device_out%zz    = devicex%data(2)
    device_out%ra    = devicex%data(3)
    device_out%rb    = devicex%data(4)
    device_out%bb    = devicex%data(5)
    device_out%ip    = devicex%data(6)
    device_out%elip  = devicex%data(7)
    device_out%trig  = devicex%data(8)
    ierr = 0

    if(bpsd_debug_flag) then
       write(6,*) '-- bpsd_get_device'
       do nd=1,devicex%ndmax
          write(6,'(A32,1PE12.4)') &
                           devicex%kid(nd),devicex%data(nd)
       enddo
    endif
    return
  end subroutine bpsd_get_device

!-----------------------------------------------------------------------
  subroutine bpsd_save_device(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    integer,intent(in) :: fid
    integer,intent(out) :: ierr

    if(bpsd_devicex_init_flag) call bpsd_devicex_init

    if(devicex%status.gt.1) &
               call bpsd_save_data0Dx(fid,devicex,ierr)
    return

  end subroutine bpsd_save_device

!-----------------------------------------------------------------------
  subroutine bpsd_load_device(datax,ierr)
!-----------------------------------------------------------------------

    use bpsd_subs
    implicit none
    type(bpsd_data0Dx_type),intent(in) :: datax
    integer,intent(out) :: ierr
    integer:: nd

    if(bpsd_devicex_init_flag) call bpsd_devicex_init

    devicex%dataName=datax%dataName
    devicex%ndmax=datax%ndmax
    devicex%time = datax%time
    do nd=1,devicex%ndmax
       devicex%data(nd) = datax%data(nd)
    enddo
    do nd=1,devicex%ndmax
       devicex%kid(nd)=datax%kid(nd)
       devicex%kunit(nd)=datax%kunit(nd)
    enddo
    devicex%created_date = datax%created_date
    devicex%created_time = datax%created_time
    devicex%created_timezone = datax%created_timezone
    devicex%status=2
    ierr=0
    return

  end subroutine bpsd_load_device

end module bpsd_device
