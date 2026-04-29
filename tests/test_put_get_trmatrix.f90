! tests/test_put_get_trmatrix.f90
! Round-trip test: put -> get for bpsd_trmatrix_type.
! Also asserts kid/kunit reflect the 6 trmatrix variables (Dn,Dp,DT,un,up,uT)
! using the bpsd_get_trmatrix_kdata accessor (no on-disk binary parsing).

program test_put_get_trmatrix
  use bpsd_kinds
  use bpsd_types
  use bpsd_trmatrix, only: bpsd_get_trmatrix_kdata
  use bpsd
  implicit none

  integer, parameter :: nrmax = 5, nsmax = 2
  type(bpsd_trmatrix_type) :: in_data, out_data
  integer :: ierr, nr, ns, nfail
  real(rkind), parameter :: tol = 1.0e-12_rkind

  in_data%nrmax = nrmax
  in_data%nsmax = nsmax
  in_data%time  = 3.5_rkind
  allocate(in_data%rho(nrmax))
  allocate(in_data%data(nrmax, nsmax))
  do nr = 1, nrmax
     in_data%rho(nr) = (nr - 1) / real(nrmax - 1, rkind)
     do ns = 1, nsmax
        in_data%data(nr,ns)%Dn = 1.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns
        in_data%data(nr,ns)%Dp = 2.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns
        in_data%data(nr,ns)%DT = 3.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns
        in_data%data(nr,ns)%un = 4.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns
        in_data%data(nr,ns)%up = 5.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns
        in_data%data(nr,ns)%uT = 6.0_rkind + 0.01_rkind * nr + 0.1_rkind * ns
     end do
  end do

  call bpsd_put_data(in_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: put ierr=', ierr
     error stop 1
  end if

  call bpsd_get_data(out_data, ierr)
  if (ierr /= 0) then
     write(*,*) 'FAIL: get ierr=', ierr
     error stop 2
  end if

  nfail = 0
  if (out_data%nrmax /= nrmax) then
     write(*,*) 'FAIL: nrmax', out_data%nrmax, '/=', nrmax; nfail = nfail + 1
  end if
  if (out_data%nsmax /= nsmax) then
     write(*,*) 'FAIL: nsmax', out_data%nsmax, '/=', nsmax; nfail = nfail + 1
  end if
  if (abs(out_data%time - in_data%time) > tol) nfail = nfail + 1
  do nr = 1, nrmax
     if (abs(out_data%rho(nr) - in_data%rho(nr)) > tol) nfail = nfail + 1
     do ns = 1, nsmax
        if (abs(out_data%data(nr,ns)%Dn - in_data%data(nr,ns)%Dn) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%Dp - in_data%data(nr,ns)%Dp) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%DT - in_data%data(nr,ns)%DT) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%un - in_data%data(nr,ns)%un) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%up - in_data%data(nr,ns)%up) > tol) nfail = nfail + 1
        if (abs(out_data%data(nr,ns)%uT - in_data%data(nr,ns)%uT) > tol) nfail = nfail + 1
     end do
  end do

  call check_kid_kunit(nfail)

  if (nfail /= 0) then
     write(*,*) 'FAIL: test_put_get_trmatrix mismatches=', nfail
     error stop 99
  end if
  write(*,*) 'PASS: test_put_get_trmatrix'
contains

  subroutine check_kid_kunit(fail_count)
    integer, intent(inout) :: fail_count
    integer :: ndmax_out, ios, nd, k, ns
    character(len=32), allocatable :: kid_out(:), kunit_out(:)
    character(len=32), parameter :: expected_kid(6) = [ &
       character(len=32) ::                                  &
       'trmatrix%Dn', 'trmatrix%Dp', 'trmatrix%DT',          &
       'trmatrix%un', 'trmatrix%up', 'trmatrix%uT'           ]
    character(len=32), parameter :: expected_kunit(6) = [ &
       character(len=32) ::                                  &
       'm^2/s', 'm^2/s', 'm^2/s',                            &
       'm/s',   'm/s',   'm/s'                               ]

    call bpsd_get_trmatrix_kdata(ndmax_out, kid_out, kunit_out, ios)
    if (ios /= 0) then
       write(*,*) 'FAIL: bpsd_get_trmatrix_kdata ierr=', ios
       fail_count = fail_count + 1
       return
    end if
    if (ndmax_out /= nsmax * 6) then
       write(*,*) 'FAIL: trmatrix ndmax', ndmax_out, '/=', nsmax * 6
       fail_count = fail_count + 1
    end if
    do ns = 1, nsmax
       do k = 1, 6
          nd = (ns - 1) * 6 + k
          if (trim(kid_out(nd)) /= trim(expected_kid(k))) then
             write(*,*) 'FAIL: kid', nd, '=[', trim(kid_out(nd)), &
                        '] expected [', trim(expected_kid(k)), ']'
             fail_count = fail_count + 1
          end if
          if (trim(kunit_out(nd)) /= trim(expected_kunit(k))) then
             write(*,*) 'FAIL: kunit', nd, '=[', trim(kunit_out(nd)), &
                        '] expected [', trim(expected_kunit(k)), ']'
             fail_count = fail_count + 1
          end if
       end do
    end do
    if (allocated(kid_out))   deallocate(kid_out)
    if (allocated(kunit_out)) deallocate(kunit_out)
  end subroutine check_kid_kunit

end program test_put_get_trmatrix
