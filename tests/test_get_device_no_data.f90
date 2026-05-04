! tests/test_get_device_no_data.f90
! Regression test for the bpsd_get_device fresh-init garbage bug.
!
! Pre-fix behaviour: bpsd_devicex_init allocates devicex%data(8) but
! does NOT zero-init, and sets status=0. bpsd_get_device only checked
! status==1, so a fresh init followed by bpsd_get_data(device, ierr)
! would silently return ierr=0 with whatever bytes happened to be in
! the just-allocated data array — i.e. garbage with a "success" code.
!
! Post-fix behaviour: bpsd_get_device rejects status<2 (covers both
! the never-allocated and the allocated-but-never-put cases),
! returning ierr=2.
!
! This test exercises only the no-data path (status=0). The normal
! put -> get round trip is implicitly covered by every other
! test_put_get_*.f90 that exercises bpsd_put/get pairs.

program test_get_device_no_data
  use bpsd_kinds
  use bpsd_types
  use bpsd
  implicit none

  type(bpsd_device_type) :: out_data
  integer :: ierr

  ! No bpsd_put_device call — fresh BPSD state with status=0.
  call bpsd_get_data(out_data, ierr)

  if (ierr == 0) then
     write(*,*) 'FAIL: bpsd_get_device returned ierr=0 on fresh init '// &
                '(should reject status<2; pre-fix bug)'
     error stop 1
  end if

  write(*,*) 'PASS: test_get_device_no_data (ierr=', ierr, ')'
end program test_get_device_no_data
