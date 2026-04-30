! tests/test_put_get_metric1D.f90
! Round-trip test: put -> get for bpsd_metric1D_type.

program test_put_get_metric1D
  use bpsd_kinds
  use bpsd_types
  use bpsd
  implicit none

  integer, parameter :: nrmax = 5
  type(bpsd_metric1D_type) :: in_data, out_data
  integer :: ierr, nr, nfail
  real(rkind), parameter :: tol = 1.0e-12_rkind

  in_data%nrmax = nrmax
  in_data%time  = 0.5_rkind
  allocate(in_data%rho(nrmax))
  allocate(in_data%data(nrmax))
  do nr = 1, nrmax
     in_data%rho(nr)            = (nr - 1) / real(nrmax - 1, rkind)
     in_data%data(nr)%pvol      = 1.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%psur      = 2.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%dvpsit    = 3.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%dvpsip    = 4.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%aver2     = 5.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%aver2i    = 6.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%aveb2     = 7.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%aveb2i    = 8.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%avegv     = 9.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%avegv2    = 10.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%avegvr2   = 11.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%avegr     = 12.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%avegr2    = 13.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%avegrr2   = 14.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%avegpp2   = 15.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%rr        = 16.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%rs        = 17.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%elip      = 18.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%trig      = 19.0_rkind + 0.1_rkind  * nr
     in_data%data(nr)%aveb      = 20.0_rkind + 0.1_rkind  * nr
  end do

  call bpsd_put_data(in_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: put ierr=', ierr; error stop 1
  end if

  out_data%nrmax = 0
  call bpsd_get_data(out_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: get ierr=', ierr; error stop 2
  end if

  nfail = 0
  if (out_data%nrmax /= nrmax) nfail = nfail + 1
  if (abs(out_data%time - in_data%time) > tol) nfail = nfail + 1
  do nr = 1, nrmax
     if (abs(out_data%rho(nr)          - in_data%rho(nr))          > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%pvol    - in_data%data(nr)%pvol)    > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%psur    - in_data%data(nr)%psur)    > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%dvpsit  - in_data%data(nr)%dvpsit)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%dvpsip  - in_data%data(nr)%dvpsip)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%aver2   - in_data%data(nr)%aver2)   > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%aver2i  - in_data%data(nr)%aver2i)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%aveb2   - in_data%data(nr)%aveb2)   > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%aveb2i  - in_data%data(nr)%aveb2i)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%avegv   - in_data%data(nr)%avegv)   > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%avegv2  - in_data%data(nr)%avegv2)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%avegvr2 - in_data%data(nr)%avegvr2) > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%avegr   - in_data%data(nr)%avegr)   > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%avegr2  - in_data%data(nr)%avegr2)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%avegrr2 - in_data%data(nr)%avegrr2) > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%avegpp2 - in_data%data(nr)%avegpp2) > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%rr      - in_data%data(nr)%rr)      > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%rs      - in_data%data(nr)%rs)      > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%elip    - in_data%data(nr)%elip)    > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%trig    - in_data%data(nr)%trig)    > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%aveb    - in_data%data(nr)%aveb)    > tol) nfail = nfail + 1
  end do

  if (nfail /= 0) then
     write(*,*) 'FAIL: test_put_get_metric1D mismatches=', nfail; error stop 99
  end if
  write(*,*) 'PASS: test_put_get_metric1D'
end program test_put_get_metric1D
