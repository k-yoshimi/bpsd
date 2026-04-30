! tests/test_save_load_roundtrip.f90
! Save/load round-trip test for plasmaf, trmatrix, trsource.
! Pattern: put v1 -> save -> put v2 -> load -> get must equal v1.

program test_save_load_roundtrip
  use bpsd_kinds
  use bpsd_types
  use bpsd
  implicit none

  integer, parameter :: nrmax = 4, nsmax = 2
  type(bpsd_plasmaf_type)  :: pf1, pf2, pfo
  type(bpsd_trmatrix_type) :: tm1, tm2, tmo
  type(bpsd_trsource_type) :: ts1, ts2, tso
  integer :: ierr, nfail
  real(rkind), parameter :: tol = 1.0e-12_rkind

  call build_plasmaf(pf1, 1.0_rkind)
  call build_plasmaf(pf2, 7.0_rkind)
  call build_trmatrix(tm1, 1.0_rkind)
  call build_trmatrix(tm2, 7.0_rkind)
  call build_trsource(ts1, 1.0_rkind)
  call build_trsource(ts2, 7.0_rkind)

  ! Put v1 and save.
  call bpsd_put_data(pf1, ierr); if (ierr /= 0) call die('put pf1', ierr)
  call bpsd_put_data(tm1, ierr); if (ierr /= 0) call die('put tm1', ierr)
  call bpsd_put_data(ts1, ierr); if (ierr /= 0) call die('put ts1', ierr)
  call bpsd_save(ierr);          if (ierr /= 0) call die('bpsd_save', ierr)

  ! Overwrite internal state with v2.
  call bpsd_put_data(pf2, ierr); if (ierr /= 0) call die('put pf2', ierr)
  call bpsd_put_data(tm2, ierr); if (ierr /= 0) call die('put tm2', ierr)
  call bpsd_put_data(ts2, ierr); if (ierr /= 0) call die('put ts2', ierr)

  ! Load file -> internal state should be v1 again.
  call bpsd_load(ierr);          if (ierr /= 0) call die('bpsd_load', ierr)

  call bpsd_get_data(pfo, ierr); if (ierr /= 0) call die('get pfo', ierr)
  tmo%nrmax = 0
  call bpsd_get_data(tmo, ierr); if (ierr /= 0) call die('get tmo', ierr)
  call bpsd_get_data(tso, ierr); if (ierr /= 0) call die('get tso', ierr)

  nfail = 0
  call cmp_plasmaf (pf1, pfo, nfail)
  call cmp_trmatrix(tm1, tmo, nfail)
  call cmp_trsource(ts1, tso, nfail)

  if (nfail /= 0) then
     write(*,*) 'FAIL: test_save_load_roundtrip mismatches=', nfail
     error stop 99
  end if
  write(*,*) 'PASS: test_save_load_roundtrip'

contains

  subroutine die(msg, ierr_val)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: ierr_val
    write(*,*) 'FAIL:', msg, ' ierr=', ierr_val
    error stop 1
  end subroutine die

  subroutine build_plasmaf(p, scale)
    type(bpsd_plasmaf_type), intent(out) :: p
    real(rkind), intent(in) :: scale
    integer :: nr, ns
    p%nrmax = nrmax; p%nsmax = nsmax; p%time = scale
    allocate(p%rho(nrmax), p%qinv(nrmax), p%data(nrmax, nsmax))
    do nr = 1, nrmax
       p%rho(nr)  = (nr - 1) / real(nrmax - 1, rkind)
       p%qinv(nr) = 0.5_rkind + scale * 0.01_rkind * nr
       do ns = 1, nsmax
          p%data(nr,ns)%density          = scale * 1.0e19_rkind * nr * ns
          p%data(nr,ns)%temperature      = scale * 1000.0_rkind * (nr + ns)
          p%data(nr,ns)%temperature_para = scale * 1100.0_rkind * (nr + ns)
          p%data(nr,ns)%temperature_perp = scale * 1200.0_rkind * (nr + ns)
          p%data(nr,ns)%velocity_tor     = scale * 100.0_rkind  * (nr + ns)
          p%data(nr,ns)%velocity_pol     = scale * 200.0_rkind  * (nr + ns)
          p%data(nr,ns)%velocity_para    = scale * 300.0_rkind  * (nr + ns)
          p%data(nr,ns)%velocity_perp    = scale * 400.0_rkind  * (nr + ns)
          p%data(nr,ns)%zave             = scale * 1.0_rkind    * ns
          p%data(nr,ns)%z2ave            = scale * 1.0_rkind    * ns * ns
          p%data(nr,ns)%density_fastion  = scale * 1.0e17_rkind * nr * ns
          p%data(nr,ns)%energy_fastion   = scale * 5000.0_rkind * (nr + ns)
       end do
    end do
  end subroutine build_plasmaf

  subroutine build_trmatrix(t, scale)
    type(bpsd_trmatrix_type), intent(out) :: t
    real(rkind), intent(in) :: scale
    integer :: nr, ns
    t%nrmax = nrmax; t%nsmax = nsmax; t%time = scale
    allocate(t%rho(nrmax), t%data(nrmax, nsmax))
    do nr = 1, nrmax
       t%rho(nr) = (nr - 1) / real(nrmax - 1, rkind)
       do ns = 1, nsmax
          t%data(nr,ns)%Dn = scale * (1.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns)
          t%data(nr,ns)%Dp = scale * (2.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns)
          t%data(nr,ns)%DT = scale * (3.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns)
          t%data(nr,ns)%un = scale * (4.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns)
          t%data(nr,ns)%up = scale * (5.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns)
          t%data(nr,ns)%uT = scale * (6.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns)
       end do
    end do
  end subroutine build_trmatrix

  subroutine build_trsource(s, scale)
    type(bpsd_trsource_type), intent(out) :: s
    real(rkind), intent(in) :: scale
    integer :: nr, ns
    s%nrmax = nrmax; s%nsmax = nsmax; s%time = scale
    allocate(s%rho(nrmax), s%data(nrmax, nsmax))
    do nr = 1, nrmax
       s%rho(nr) = (nr - 1) / real(nrmax - 1, rkind)
       do ns = 1, nsmax
          s%data(nr,ns)%nip = scale * (1.0e18_rkind + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%nim = scale * (2.0e18_rkind + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%ncx = scale * (3.0e18_rkind + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%Pec = scale * (1.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%Plh = scale * (2.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%Pic = scale * (3.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%Pbr = scale * (4.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%Pcy = scale * (5.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%Plr = scale * (6.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns)
          s%data(nr,ns)%Poh = scale * (7.0e6_rkind  + 0.1_rkind * nr + 0.01_rkind * ns)
       end do
    end do
  end subroutine build_trsource

  subroutine cmp_plasmaf(a, b, fc)
    type(bpsd_plasmaf_type), intent(in) :: a, b
    integer, intent(inout) :: fc
    integer :: nr, ns
    if (b%nrmax /= a%nrmax) fc = fc + 1
    if (b%nsmax /= a%nsmax) fc = fc + 1
    if (abs(b%time - a%time) > tol) fc = fc + 1
    do nr = 1, a%nrmax
       if (abs(b%rho(nr)  - a%rho(nr))  > tol) fc = fc + 1
       if (abs(b%qinv(nr) - a%qinv(nr)) > tol) fc = fc + 1
       do ns = 1, a%nsmax
          if (abs(b%data(nr,ns)%density          - a%data(nr,ns)%density)          > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%temperature      - a%data(nr,ns)%temperature)      > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%temperature_para - a%data(nr,ns)%temperature_para) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%temperature_perp - a%data(nr,ns)%temperature_perp) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%velocity_tor     - a%data(nr,ns)%velocity_tor)     > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%velocity_pol     - a%data(nr,ns)%velocity_pol)     > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%velocity_para    - a%data(nr,ns)%velocity_para)    > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%velocity_perp    - a%data(nr,ns)%velocity_perp)    > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%zave             - a%data(nr,ns)%zave)             > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%z2ave            - a%data(nr,ns)%z2ave)            > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%density_fastion  - a%data(nr,ns)%density_fastion)  > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%energy_fastion   - a%data(nr,ns)%energy_fastion)   > tol) fc = fc + 1
       end do
    end do
  end subroutine cmp_plasmaf

  subroutine cmp_trmatrix(a, b, fc)
    type(bpsd_trmatrix_type), intent(in) :: a, b
    integer, intent(inout) :: fc
    integer :: nr, ns
    if (b%nrmax /= a%nrmax) fc = fc + 1
    if (b%nsmax /= a%nsmax) fc = fc + 1
    if (abs(b%time - a%time) > tol) fc = fc + 1
    do nr = 1, a%nrmax
       if (abs(b%rho(nr) - a%rho(nr)) > tol) fc = fc + 1
       do ns = 1, a%nsmax
          if (abs(b%data(nr,ns)%Dn - a%data(nr,ns)%Dn) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Dp - a%data(nr,ns)%Dp) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%DT - a%data(nr,ns)%DT) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%un - a%data(nr,ns)%un) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%up - a%data(nr,ns)%up) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%uT - a%data(nr,ns)%uT) > tol) fc = fc + 1
       end do
    end do
  end subroutine cmp_trmatrix

  subroutine cmp_trsource(a, b, fc)
    type(bpsd_trsource_type), intent(in) :: a, b
    integer, intent(inout) :: fc
    integer :: nr, ns
    if (b%nrmax /= a%nrmax) fc = fc + 1
    if (b%nsmax /= a%nsmax) fc = fc + 1
    if (abs(b%time - a%time) > tol) fc = fc + 1
    do nr = 1, a%nrmax
       if (abs(b%rho(nr) - a%rho(nr)) > tol) fc = fc + 1
       do ns = 1, a%nsmax
          if (abs(b%data(nr,ns)%nip - a%data(nr,ns)%nip) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%nim - a%data(nr,ns)%nim) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%ncx - a%data(nr,ns)%ncx) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Pec - a%data(nr,ns)%Pec) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Plh - a%data(nr,ns)%Plh) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Pic - a%data(nr,ns)%Pic) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Pbr - a%data(nr,ns)%Pbr) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Pcy - a%data(nr,ns)%Pcy) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Plr - a%data(nr,ns)%Plr) > tol) fc = fc + 1
          if (abs(b%data(nr,ns)%Poh - a%data(nr,ns)%Poh) > tol) fc = fc + 1
       end do
    end do
  end subroutine cmp_trsource

end program test_save_load_roundtrip
