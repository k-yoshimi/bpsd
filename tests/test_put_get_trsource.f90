! tests/test_put_get_trsource.f90
! Round-trip test: put -> get for bpsd_trsource_type.

program test_put_get_trsource
  use bpsd_kinds
  use bpsd_types
  use bpsd
  implicit none

  integer, parameter :: nrmax = 5, nsmax = 2
  type(bpsd_trsource_type) :: in_data, out_data
  integer :: ierr, nr, ns, nfail
  real(rkind), parameter :: tol = 1.0e-12_rkind

  in_data%nrmax = nrmax
  in_data%nsmax = nsmax
  in_data%time  = 4.5_rkind
  allocate(in_data%rho(nrmax))
  allocate(in_data%data(nrmax, nsmax))
  do nr = 1, nrmax
     in_data%rho(nr) = (nr - 1) / real(nrmax - 1, rkind)
     do ns = 1, nsmax
        in_data%data(nr,ns)%nip = 1.0e18_rkind + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%nim = 2.0e18_rkind + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%ncx = 3.0e18_rkind + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%Pec = 1.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%Plh = 2.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%Pic = 3.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%Pbr = 4.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%Pcy = 5.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%Plr = 6.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns
        in_data%data(nr,ns)%Poh = 7.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns
     end do
  end do

  call bpsd_put_data(in_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: put ierr=', ierr; error stop 1
  end if

  call bpsd_get_data(out_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: get ierr=', ierr; error stop 2
  end if

  nfail = 0
  if (out_data%nrmax /= nrmax) nfail = nfail + 1
  if (out_data%nsmax /= nsmax) nfail = nfail + 1
  if (abs(out_data%time - in_data%time) > tol) nfail = nfail + 1
  do nr = 1, nrmax
     if (abs(out_data%rho(nr) - in_data%rho(nr)) > tol) nfail = nfail + 1
     do ns = 1, nsmax
        if (abs(out_data%data(nr,ns)%nip - in_data%data(nr,ns)%nip) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%nim - in_data%data(nr,ns)%nim) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%ncx - in_data%data(nr,ns)%ncx) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Pec - in_data%data(nr,ns)%Pec) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Plh - in_data%data(nr,ns)%Plh) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Pic - in_data%data(nr,ns)%Pic) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Pbr - in_data%data(nr,ns)%Pbr) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Pcy - in_data%data(nr,ns)%Pcy) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Plr - in_data%data(nr,ns)%Plr) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Poh - in_data%data(nr,ns)%Poh) > tol) nfail = nfail + 1
     end do
  end do

  if (nfail /= 0) then
     write(*,*) 'FAIL: test_put_get_trsource mismatches=', nfail; error stop 99
  end if
  write(*,*) 'PASS: test_put_get_trsource'
end program test_put_get_trsource
