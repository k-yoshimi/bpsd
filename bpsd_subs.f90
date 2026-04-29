! bpsd_subs.f90

MODULE bpsd_subs

CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_karray(data,n1)
!-----------------------------------------------------------------------
    USE bpsd_kinds
    IMPLICIT NONE
    CHARACTER(LEN=*),DIMENSION(:),ALLOCATABLE:: data
    INTEGER(ikind),INTENT(IN):: n1

    IF(ALLOCATED(data)) THEN
       IF(n1.LE.0) THEN
          DEALLOCATE(data)
       ELSE IF(n1.NE.SIZE(data,1)) THEN
          DEALLOCATE(data)
          ALLOCATE(data(n1))
       END IF
    ELSE
       IF(n1.GT.0)  THEN
          ALLOCATE(data(n1))
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE bpsd_adjust_karray

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_array1D(data,n1)
!-----------------------------------------------------------------------
    USE bpsd_kinds
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: data
    INTEGER(ikind),INTENT(IN):: n1

    IF(ALLOCATED(data)) THEN
       IF(n1.LE.0) THEN
          DEALLOCATE(data)
       ELSE IF(n1.NE.SIZE(data,1)) THEN
          DEALLOCATE(data)
          ALLOCATE(data(n1))
       END IF
    ELSE
       IF(n1.GT.0)  THEN
          ALLOCATE(data(n1))
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE bpsd_adjust_array1D

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_array2D(data,n1,n2)
!-----------------------------------------------------------------------
    USE bpsd_kinds
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: data
    INTEGER(ikind),INTENT(IN):: n1,n2

    IF(ALLOCATED(data)) THEN
       IF(n1.LE.0.or.n2.LE.0) THEN
          DEALLOCATE(data)
       ELSE IF(n1.NE.SIZE(data,1).OR.n2.NE.SIZE(data,2)) THEN
          DEALLOCATE(data)
          ALLOCATE(data(n1,n2))
       END IF
    ELSE
       IF(n1.GT.0.AND.n2.GT.0)  THEN
          ALLOCATE(data(n1,n2))
       ENDIF
    ENDIF
  END SUBROUTINE bpsd_adjust_array2D

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_adjust_array3D(data,n1,n2,n3)
!-----------------------------------------------------------------------
    USE bpsd_kinds
    IMPLICIT NONE
    REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: data
    INTEGER(ikind),INTENT(IN):: n1,n2,n3

    IF(ALLOCATED(data)) THEN
       IF(n1.LE.0.OR.n2.LE.0.OR.n3.LE.0) THEN
          DEALLOCATE(data)
       ELSE IF(n1.NE.SIZE(data,1).OR.n2.NE.SIZE(data,2).OR. &
               n3.NE.SIZE(data,3)) THEN
          DEALLOCATE(data)
          ALLOCATE(data(n1,n2,n3))
       END IF
    ELSE
       IF(n1.GT.0.AND.n2.GT.0.AND.n3.GT.0)  THEN
          ALLOCATE(data(n1,n2,n3))
       ENDIF
    ENDIF
  END SUBROUTINE bpsd_adjust_array3D

!-----------------------------------------------------------------------
  SUBROUTINE bpsd_spl1D(data1D,nd,ierr)
!-----------------------------------------------------------------------
    USE bpsd_types_internal
    USE bpsd_libspl
    IMPLICIT NONE
    TYPE(bpsd_data1Dx_type) :: data1D
    INTEGER(ikind) :: nd      ! position of dependent variable
    INTEGER(ikind) :: ierr    ! error indicator
    REAL(rkind), DIMENSION(:), ALLOCATABLE :: deriv

    allocate(deriv(data1D%nrmax))
    call spl1D(data1D%s,data1D%data(:,nd),deriv,data1D%spline(:,:,nd), &
               data1D%nrmax,0,ierr)
    if(ierr.ne.0) &
               write(6,*) 'XX bpsd_spl1D : '//data1D%kid(nd)//': ierr=',ierr
    deallocate(deriv)
    return
  end subroutine bpsd_spl1D

!-----------------------------------------------------------------------
  subroutine bpsd_spl1DF(pos,val,data1D,nd,ierr)
!-----------------------------------------------------------------------
    use bpsd_types_internal
    USE bpsd_libspl
    implicit none
    real(rkind) :: pos     ! value of independent variable
    real(rkind) :: val     ! value of dependent variable
    type(bpsd_data1Dx_type) :: data1D
    integer :: nd      ! position of dependent variable
    integer :: ierr    ! error indicator

    call spl1DF(pos**2,val,data1D%s,data1D%spline(:,:,nd),data1D%nrmax,ierr)
    if(ierr.ne.0) then
       write(6,*) 'XX bpsd_spl1DF : '//data1D%kid(nd)//': ierr=',ierr
       write(6,'(1P3E12.4)')  &
                     pos**2,data1D%s(1),data1D%s(data1D%nrmax)
    endif
  end subroutine bpsd_spl1DF

!-----------------------------------------------------------------------
  subroutine bpsd_save_shotx(fid,datax,ierr)
!-----------------------------------------------------------------------
    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_shotx_type),intent(in):: datax
    integer,intent(out):: ierr

    write(fid,IOSTAT=ierr,ERR=8) 'data:str'
    write(fid,IOSTAT=ierr,ERR=8) datax%dataName
    write(fid,IOSTAT=ierr,ERR=8) datax%deviceID
    write(fid,IOSTAT=ierr,ERR=8) datax%shotID,datax%modelID

    write(6,*) '+++ ',datax%dataName,': saved'
    ierr=0
    return

8   continue
    return
  end subroutine bpsd_save_shotx

!-----------------------------------------------------------------------
  subroutine bpsd_save_data0Dx(fid,datax,ierr)
!-----------------------------------------------------------------------
    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data0Dx_type),intent(in):: datax
    integer,intent(out):: ierr

    write(fid,IOSTAT=ierr,ERR=8) 'data::0D'
    write(fid,IOSTAT=ierr,ERR=8) datax%dataName
    write(fid,IOSTAT=ierr,ERR=8) datax%time
    write(fid,IOSTAT=ierr,ERR=8) datax%ndmax
    write(fid,IOSTAT=ierr,ERR=8) datax%created_date, &
                                 datax%created_time, &
                                 datax%created_timezone
    write(fid,IOSTAT=ierr,ERR=8) datax%kid
    write(fid,IOSTAT=ierr,ERR=8) datax%kunit
    write(fid,IOSTAT=ierr,ERR=8) datax%data

    write(6,*) '+++ ',datax%dataName,': saved'
    ierr=0
    return

8   continue
    write(6,*) 'XXX', ierr
    return
  end subroutine bpsd_save_data0Dx

!-----------------------------------------------------------------------
  subroutine bpsd_save_data1Dx(fid,datax,ierr)
!-----------------------------------------------------------------------
    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data1Dx_type),intent(in):: datax
    integer,intent(out):: ierr

    write(fid,IOSTAT=ierr,ERR=8) 'data::1D'
    write(fid,IOSTAT=ierr,ERR=8) datax%dataName
    write(fid,IOSTAT=ierr,ERR=8) datax%time
    write(fid,IOSTAT=ierr,ERR=8) datax%nrmax,datax%ndmax
    write(fid,IOSTAT=ierr,ERR=8) datax%created_date, &
                                 datax%created_time, &
                                 datax%created_timezone
    write(fid,IOSTAT=ierr,ERR=8) datax%kid
    write(fid,IOSTAT=ierr,ERR=8) datax%kunit
    write(fid,IOSTAT=ierr,ERR=8) datax%rho
    write(fid,IOSTAT=ierr,ERR=8) datax%data

    write(6,*) '+++ ',datax%dataName,': saved'
    ierr=0
    return

8   continue
    return
  end subroutine bpsd_save_data1Dx

!-----------------------------------------------------------------------
  subroutine bpsd_save_data2Dx(fid,datax,ierr)
!-----------------------------------------------------------------------
    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data2Dx_type),intent(in):: datax
    integer,intent(out):: ierr

    write(fid,IOSTAT=ierr,ERR=8) 'data::2D'
    write(fid,IOSTAT=ierr,ERR=8) datax%dataName
    write(fid,IOSTAT=ierr,ERR=8) datax%time
    write(fid,IOSTAT=ierr,ERR=8) datax%nrmax,datax%nthmax,datax%ndmax
    write(fid,IOSTAT=ierr,ERR=8) datax%created_date, &
                                 datax%created_time, &
                                 datax%created_timezone
    write(fid,IOSTAT=ierr,ERR=8) datax%kid
    write(fid,IOSTAT=ierr,ERR=8) datax%rho
    write(fid,IOSTAT=ierr,ERR=8) datax%th
    write(fid,IOSTAT=ierr,ERR=8) datax%data

    write(6,*) '+++ ',datax%dataName,': saved'
    ierr=0
    return

8   continue
    return
  end subroutine bpsd_save_data2Dx

!-----------------------------------------------------------------------
  subroutine bpsd_save_data3Dx(fid,datax,ierr)
!-----------------------------------------------------------------------
    use bpsd_types_internal
    implicit none
    integer,intent(in):: fid
    type(bpsd_data3Dx_type),intent(in):: datax
    integer,intent(out):: ierr

    write(fid,IOSTAT=ierr,ERR=8) 'data::3D'
    write(fid,IOSTAT=ierr,ERR=8) datax%dataName
    write(fid,IOSTAT=ierr,ERR=8) datax%time
    write(fid,IOSTAT=ierr,ERR=8) &
               datax%nrmax,datax%nthmax,datax%nphmax,datax%ndmax
    write(fid,IOSTAT=ierr,ERR=8) datax%created_date, &
                                 datax%created_time, &
                                 datax%created_timezone
    write(fid,IOSTAT=ierr,ERR=8) datax%kid
    write(fid,IOSTAT=ierr,ERR=8) datax%rho
    write(fid,IOSTAT=ierr,ERR=8) datax%th
    write(fid,IOSTAT=ierr,ERR=8) datax%ph
    write(fid,IOSTAT=ierr,ERR=8) datax%data

    write(6,*) '+++ ',datax%dataName,': saved'
    ierr=0
    return

8   continue
    return
  end subroutine bpsd_save_data3Dx
end module bpsd_subs
