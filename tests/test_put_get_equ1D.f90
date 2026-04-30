! tests/test_put_get_equ1D.f90
! Round-trip test: put -> get for bpsd_equ1D_type.

program test_put_get_equ1D
  use bpsd_kinds
  use bpsd_types
  use bpsd
  implicit none

  integer, parameter :: nrmax = 7
  type(bpsd_equ1D_type) :: in_data, out_data
  integer :: ierr, nr, nfail
  real(rkind), parameter :: tol = 1.0e-12_rkind

  in_data%nrmax = nrmax
  in_data%time  = 1.25_rkind
  allocate(in_data%rho(nrmax))
  allocate(in_data%data(nrmax))
  do nr = 1, nrmax
     in_data%rho(nr)        = (nr - 1) / real(nrmax - 1, rkind)
     in_data%data(nr)%psit  = 1.0_rkind  + 0.1_rkind  * nr
     in_data%data(nr)%psip  = 2.0_rkind  + 0.2_rkind  * nr
     in_data%data(nr)%ppp   = 100.0_rkind + 5.0_rkind * nr
     in_data%data(nr)%piq   = 0.5_rkind  + 0.05_rkind * nr
     in_data%data(nr)%pip   = 1.0e6_rkind  + 1.0e3_rkind * nr
     in_data%data(nr)%pit   = 2.0e6_rkind  + 2.0e3_rkind * nr
  end do

  call bpsd_put_data(in_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: bpsd_put_data ierr=', ierr
     error stop 1
  end if

  out_data%nrmax = 0
  call bpsd_get_data(out_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: bpsd_get_data ierr=', ierr
     error stop 2
  end if

  nfail = 0
  if (out_data%nrmax /= nrmax) then
     write(*,*) 'FAIL: nrmax', out_data%nrmax, '/=', nrmax
     nfail = nfail + 1
  end if
  if (abs(out_data%time - in_data%time) > tol) then
     write(*,*) 'FAIL: time mismatch'
     nfail = nfail + 1
  end if
  do nr = 1, nrmax
     if (abs(out_data%rho(nr)        - in_data%rho(nr))        > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%psit  - in_data%data(nr)%psit)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%psip  - in_data%data(nr)%psip)  > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%ppp   - in_data%data(nr)%ppp)   > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%piq   - in_data%data(nr)%piq)   > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%pip   - in_data%data(nr)%pip)   > tol) nfail = nfail + 1
     if (abs(out_data%data(nr)%pit   - in_data%data(nr)%pit)   > tol) nfail = nfail + 1
  end do

  if (nfail /= 0) then
     write(*,*) 'FAIL: test_put_get_equ1D mismatches=', nfail
     error stop 99
  end if
  write(*,*) 'PASS: test_put_get_equ1D'
end program test_put_get_equ1D
