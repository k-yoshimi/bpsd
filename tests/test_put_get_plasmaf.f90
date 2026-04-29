! tests/test_put_get_plasmaf.f90
! Round-trip test: put -> get for bpsd_plasmaf_type.

program test_put_get_plasmaf
  use bpsd_kinds
  use bpsd_types
  use bpsd
  implicit none

  integer, parameter :: nrmax = 6, nsmax = 3
  type(bpsd_plasmaf_type) :: in_data, out_data
  integer :: ierr, nr, ns, nfail
  real(rkind), parameter :: tol = 1.0e-12_rkind

  in_data%nrmax = nrmax
  in_data%nsmax = nsmax
  in_data%time  = 2.0_rkind
  allocate(in_data%rho(nrmax))
  allocate(in_data%qinv(nrmax))
  allocate(in_data%data(nrmax, nsmax))
  do nr = 1, nrmax
     in_data%rho(nr)  = (nr - 1) / real(nrmax - 1, rkind)
     in_data%qinv(nr) = 0.5_rkind + 0.05_rkind * nr
     do ns = 1, nsmax
        in_data%data(nr,ns)%density          = 1.0e19_rkind * nr * ns
        in_data%data(nr,ns)%temperature      = 1000.0_rkind * (nr + ns)
        in_data%data(nr,ns)%temperature_para = 1100.0_rkind * (nr + ns)
        in_data%data(nr,ns)%temperature_perp = 1200.0_rkind * (nr + ns)
        in_data%data(nr,ns)%velocity_tor     = 100.0_rkind  * (nr + ns)
        in_data%data(nr,ns)%velocity_pol     = 200.0_rkind  * (nr + ns)
        in_data%data(nr,ns)%velocity_para    = 300.0_rkind  * (nr + ns)
        in_data%data(nr,ns)%velocity_perp    = 400.0_rkind  * (nr + ns)
        in_data%data(nr,ns)%zave             = 1.0_rkind    * ns
        in_data%data(nr,ns)%z2ave            = 1.0_rkind    * ns * ns
        in_data%data(nr,ns)%density_fastion  = 1.0e17_rkind * nr * ns
        in_data%data(nr,ns)%energy_fastion   = 5000.0_rkind * (nr + ns)
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
     if (abs(out_data%rho(nr)  - in_data%rho(nr))  > tol) nfail = nfail + 1
     if (abs(out_data%qinv(nr) - in_data%qinv(nr)) > tol) nfail = nfail + 1
     do ns = 1, nsmax
        if (abs(out_data%data(nr,ns)%density          - in_data%data(nr,ns)%density)          > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%temperature      - in_data%data(nr,ns)%temperature)      > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%temperature_para - in_data%data(nr,ns)%temperature_para) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%temperature_perp - in_data%data(nr,ns)%temperature_perp) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%velocity_tor     - in_data%data(nr,ns)%velocity_tor)     > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%velocity_pol     - in_data%data(nr,ns)%velocity_pol)     > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%velocity_para    - in_data%data(nr,ns)%velocity_para)    > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%velocity_perp    - in_data%data(nr,ns)%velocity_perp)    > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%zave             - in_data%data(nr,ns)%zave)             > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%z2ave            - in_data%data(nr,ns)%z2ave)            > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%density_fastion  - in_data%data(nr,ns)%density_fastion)  > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%energy_fastion   - in_data%data(nr,ns)%energy_fastion)   > tol) nfail = nfail + 1
     end do
  end do

  if (nfail /= 0) then
     write(*,*) 'FAIL: test_put_get_plasmaf mismatches=', nfail; error stop 99
  end if
  write(*,*) 'PASS: test_put_get_plasmaf'
end program test_put_get_plasmaf
