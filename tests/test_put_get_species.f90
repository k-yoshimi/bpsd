! tests/test_put_get_species.f90
! Round-trip test: put -> get for bpsd_species_type.
! Also asserts that the internal flat-array layout is exactly nsmax*3
! (3 fields per species: pa, pz, npa). The previous layout was
! nsmax*5 with stride 3, which left 2 trailing slots per species
! uninitialized and caused garbage in saved bpsd.data files.

program test_put_get_species
  use bpsd_kinds
  use bpsd_types
  use bpsd_species, only: bpsd_get_species_kdata
  use bpsd
  implicit none

  integer, parameter :: nsmax = 3
  type(bpsd_species_type) :: in_data, out_data
  integer :: ierr, ns, nfail
  real(rkind), parameter :: tol = 1.0e-12_rkind

  in_data%nsmax = nsmax
  allocate(in_data%data(nsmax))
  do ns = 1, nsmax
     in_data%data(ns)%pa  = 1.0_rkind  * ns + 0.5_rkind
     in_data%data(ns)%pz  = 1.0_rkind  * ns - 0.25_rkind
     in_data%data(ns)%npa = ns
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
  if (out_data%nsmax /= nsmax) then
     write(*,*) 'FAIL: nsmax', out_data%nsmax, '/=', nsmax
     nfail = nfail + 1
  end if
  do ns = 1, nsmax
     if (abs(out_data%data(ns)%pa - in_data%data(ns)%pa) > tol) nfail = nfail + 1
     if (abs(out_data%data(ns)%pz - in_data%data(ns)%pz) > tol) nfail = nfail + 1
     if (out_data%data(ns)%npa /= in_data%data(ns)%npa) nfail = nfail + 1
  end do

  call check_layout(nfail)

  if (nfail /= 0) then
     write(*,*) 'FAIL: test_put_get_species mismatches=', nfail
     error stop 99
  end if
  write(*,*) 'PASS: test_put_get_species'

contains

  subroutine check_layout(fail_count)
    integer, intent(inout) :: fail_count
    integer :: ndmax_out, ios, nd, k, ns
    character(len=32), allocatable :: kid_out(:), kunit_out(:)
    character(len=32), parameter :: expected_kid(3) = [ &
       character(len=32) ::                                  &
       'species%pa', 'species%pz', 'species%npa'             ]

    call bpsd_get_species_kdata(ndmax_out, kid_out, kunit_out, ios)
    if (ios /= 0) then
       write(*,*) 'FAIL: bpsd_get_species_kdata ierr=', ios
       fail_count = fail_count + 1
       return
    end if
    if (ndmax_out /= nsmax * 3) then
       write(*,*) 'FAIL: species ndmax', ndmax_out, '/=', nsmax * 3
       fail_count = fail_count + 1
    end if
    do ns = 1, nsmax
       do k = 1, 3
          nd = (ns - 1) * 3 + k
          if (nd > ndmax_out) cycle
          if (trim(kid_out(nd)) /= trim(expected_kid(k))) then
             write(*,*) 'FAIL: kid', nd, '=[', trim(kid_out(nd)), &
                        '] expected [', trim(expected_kid(k)), ']'
             fail_count = fail_count + 1
          end if
          if (trim(kunit_out(nd)) /= '') then
             write(*,*) 'FAIL: kunit', nd, '=[', trim(kunit_out(nd)), &
                        '] expected blank'
             fail_count = fail_count + 1
          end if
       end do
    end do
    if (allocated(kid_out))   deallocate(kid_out)
    if (allocated(kunit_out)) deallocate(kunit_out)
  end subroutine check_layout

end program test_put_get_species
