!     $Id$
!=======================================================================
module bpsd

  use bpsd_kinds
  use bpsd_constants
  use bpsd_flags
  use bpsd_types
  use bpsd_shot
  use bpsd_device
  use bpsd_species
  use bpsd_equ1D
  use bpsd_metric1D
  use bpsd_plasmaf
  use bpsd_trmatrix
  use bpsd_trsource

  interface bpsd_put_data
     module procedure bpsd_put_shot, &
                      bpsd_put_device, &
                      bpsd_put_species, &
                      bpsd_put_equ1D, &
                      bpsd_put_metric1D, &
                      bpsd_put_plasmaf, &
                      bpsd_put_trmatrix, &
                      bpsd_put_trsource
  end interface

  interface bpsd_get_data
     module procedure bpsd_get_shot, &
                      bpsd_get_device, &
                      bpsd_get_species, &
                      bpsd_get_equ1D, &
                      bpsd_get_metric1D, &
                      bpsd_get_plasmaf, &
                      bpsd_get_trmatrix, &
                      bpsd_get_trsource
  end interface

contains

!-----------------------------------------------------------------------
  subroutine bpsd_save(ierr)
!-----------------------------------------------------------------------

    use bpsd_types_internal
    use bpsd_libfio
    implicit none
    integer,intent(out):: ierr
    integer:: fid

    fid=31

    call fwopen(fid,'bpsd.data',0,0,'bpsd',ierr)
    if(ierr.ne.0) return

    write(fid) '##bspd.data version 0.01'

    call bpsd_save_shot(fid,ierr)
    if(ierr.ne.0) return
    call bpsd_save_device(fid,ierr)
    if(ierr.ne.0) return
    call bpsd_save_equ1D(fid,ierr)
    if(ierr.ne.0) return
    call bpsd_save_metric1D(fid,ierr)
    if(ierr.ne.0) return
    call bpsd_save_species(fid,ierr)
    if(ierr.ne.0) return
    call bpsd_save_plasmaf(fid,ierr)
    if(ierr.ne.0) return
    call bpsd_save_trmatrix(fid,ierr)
    if(ierr.ne.0) return
    call bpsd_save_trsource(fid,ierr)
    if(ierr.ne.0) return

    write(fid) 'data:end'

    close(fid)
      
    return
  end subroutine bpsd_save

!-----------------------------------------------------------------------
  subroutine bpsd_load(ierr)
!-----------------------------------------------------------------------

    use bpsd_types_internal
    use bpsd_libfio
    implicit none
    integer,intent(out):: ierr
    integer:: fid
    character(len=24):: line
    character(len=8):: datadim

    fid=32

    call fropen(fid,'bpsd.data',0,0,'bpsd',ierr)
    if(ierr.ne.0) return

    read(fid,IOSTAT=ierr,ERR=8,END=9) line
    if(line(1:12).ne.'##bspd.data ') then
       ierr=1
       goto 9000
    endif
    if(line(21:24).ne.'0.01') then
       ierr=2
       goto 9000
    endif

1   continue
    read(fid,IOSTAT=ierr,ERR=8,END=9) datadim

    if(datadim.eq.'data:str') then
       call bpsd_load_shotx(fid,ierr)
       if(ierr.ne.0) then
          ierr=2000+ierr
          goto 9000
       endif
    endif

    if(datadim.eq.'data::0D') then
       call bpsd_load_data0Dx(fid,ierr)
       if(ierr.ne.0) then
          ierr=3000+ierr
          goto 9000
       endif
    endif

    if(datadim.eq.'data::1D') then
       call bpsd_load_data1Dx(fid,ierr)
       if(ierr.ne.0) then
          ierr=4000+ierr
          goto 9000
       endif
    endif

    if(datadim.eq.'data::2D') then
       call bpsd_load_data2Dx(fid,ierr)
       if(ierr.ne.0) then
          ierr=5000+ierr
          goto 9000
       endif
    endif

    if(datadim.eq.'data::3D') then
       call bpsd_load_data3Dx(fid,ierr)
       if(ierr.ne.0) then
          ierr=6000+ierr
          goto 9000
       endif
    endif

    if(datadim.eq.'data:end') then
       ierr=0
       goto 9000
    endif
    go to 1

8   ierr=8000
    goto 9000
9   ierr=9000

9000 close(fid)
    return
  end subroutine bpsd_load

!-----------------------------------------------------------------------
  subroutine bpsd_load_shotx(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_shotx_type):: datax
    integer,intent(out):: ierr

    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%deviceID
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%shotID,datax%modelID

    if(datax%dataName(1:4).eq.'shot') call bpsd_load_shot(datax,ierr)

    write(6,*) '+++ ',datax%dataName,': loaded'
    ierr=0
    return

8   ierr=800+ierr
    return
9   ierr=900+ierr
    return
  end subroutine bpsd_load_shotx

!-----------------------------------------------------------------------
  subroutine bpsd_load_data0Dx(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data0Dx_type):: datax
    integer,intent(out):: ierr

    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%ndmax
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%created_date, &
                                      datax%created_time, &
                                      datax%created_timezone 
    allocate(datax%kid(datax%ndmax))
    allocate(datax%kunit(datax%ndmax))
    allocate(datax%data(datax%ndmax))
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kunit
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

    if(datax%dataName(1:6).eq.'device') &
               call bpsd_load_device(datax,ierr)
    if(datax%dataName(1:7).eq.'species') &
               call bpsd_load_species(datax,ierr)

    deallocate(datax%data)
    deallocate(datax%kid)

    write(6,*) '+++ ',datax%dataName,': loaded'
    ierr=0
    return

8   ierr=800+ierr
    return
9   ierr=900+ierr
    return
  end subroutine bpsd_load_data0Dx

!-----------------------------------------------------------------------
  subroutine bpsd_load_data1Dx(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data1Dx_type):: datax
    integer,intent(out):: ierr

    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%nrmax,datax%ndmax
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%created_date, &
                                      datax%created_time, &
                                      datax%created_timezone
    allocate(datax%kid(datax%ndmax))
    allocate(datax%kunit(datax%ndmax))
    allocate(datax%rho(datax%nrmax))
    allocate(datax%data(datax%nrmax,datax%ndmax))
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kunit
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%rho
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

    if(datax%dataName(1:5).eq.'equ1D') &
               call bpsd_load_equ1D(datax,ierr)
    if(datax%dataName(1:8).eq.'metric1D') &
               call bpsd_load_metric1D(datax,ierr)
    if(datax%dataName(1:7).eq.'plasmaf')  &
               call bpsd_load_plasmaf(datax,ierr)
    if(datax%dataName(1:8).eq.'trmatrix') &
               call bpsd_load_trmatrix(datax,ierr)
    if(datax%dataName(1:8).eq.'trsource') &
               call bpsd_load_trsource(datax,ierr)

    deallocate(datax%data)
    deallocate(datax%rho)
    deallocate(datax%kid)

    write(6,*) '+++ ',datax%dataName,': loaded'
    ierr=0
    return

8   ierr=800+ierr
    return
9   ierr=900+ierr
    return
  end subroutine bpsd_load_data1Dx

!-----------------------------------------------------------------------
  subroutine bpsd_load_data2Dx(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data2Dx_type):: datax
    integer,intent(out):: ierr

    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
    read(fid,IOSTAT=ierr,ERR=8,END=9) &
               datax%nrmax,datax%nthmax,datax%ndmax
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%created_date, &
                                      datax%created_time, &
                                      datax%created_timezone
    allocate(datax%kid(datax%ndmax))
    allocate(datax%kunit(datax%ndmax))
    allocate(datax%rho(datax%nrmax))
    allocate(datax%th(datax%nthmax))
    allocate(datax%data(datax%nrmax,datax%nthmax,datax%ndmax))
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kunit
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%rho
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%th
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

!      if(datax%dataName(1:7).eq.'plasmaf') &
!           call bpsd_load_plasmaf(datax,ierr)

    deallocate(datax%data)
    deallocate(datax%th)
    deallocate(datax%rho)
    deallocate(datax%kid)

    write(6,*) '+++ ',datax%dataName,': loaded'
    ierr=0
    return

8   ierr=800+ierr
    return
9   ierr=900+ierr
    return
  end subroutine bpsd_load_data2Dx

!-----------------------------------------------------------------------
  subroutine bpsd_load_data3Dx(fid,ierr)
!-----------------------------------------------------------------------

    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data3Dx_type):: datax
    integer,intent(out):: ierr

    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%dataName
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%time
    read(fid,IOSTAT=ierr,ERR=8,END=9) &
               datax%nrmax,datax%nthmax,datax%nphmax,datax%ndmax
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%created_date, &
                                      datax%created_time, &
                                      datax%created_timezone
    allocate(datax%kid(datax%ndmax))
    allocate(datax%kunit(datax%ndmax))
    allocate(datax%rho(datax%nrmax))
    allocate(datax%th(datax%nthmax))
    allocate(datax%ph(datax%nphmax))
    allocate(datax%data(datax%nrmax,datax%nthmax,datax%nphmax, &
                   datax%ndmax))
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kid
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%kunit
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%rho
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%th
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%ph
    read(fid,IOSTAT=ierr,ERR=8,END=9) datax%data

!      if(datax%dataName(1:7).eq.'plasmaf') 
!           call bpsd_load_plasmaf(datax,ierr)

    deallocate(datax%data)
    deallocate(datax%ph)
    deallocate(datax%th)
    deallocate(datax%rho)
    deallocate(datax%kid)

    write(6,*) '+++ ',datax%dataName,': loaded'
    ierr=0
    return

8   ierr=800+ierr
    return
9   ierr=900+ierr
    return
  end subroutine bpsd_load_data3Dx

end module bpsd
