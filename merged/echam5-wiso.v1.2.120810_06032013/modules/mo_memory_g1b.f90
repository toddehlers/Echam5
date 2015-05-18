MODULE mo_memory_g1b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add_stream_element

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g1b ! construct the g1b table
  PUBLIC :: destruct_g1b  ! destruct  the g1b table

  PUBLIC :: g1b           ! the g1b table

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: vof(:,:,:)
  REAL(dp), POINTER, PUBLIC :: df(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: alpsf(:,:)
  REAL(dp), POINTER, PUBLIC :: qf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xlf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xif(:,:,:)
!---wiso-code
  REAL(dp), POINTER, PUBLIC :: wisoqf(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxlf(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxif(:,:,:,:)
!---wiso-code-end
  REAL(dp), POINTER, PUBLIC :: xtf(:,:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: g1b

CONTAINS

!---wiso-code - add lnwiso and nwiso to call of subroutine
  SUBROUTINE construct_g1b (lnlon, lnlev, lntrac, lnwiso, lngl, &
                             nlon,  nlev,  ntrac,  nwiso,  ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lntrac, lnwiso, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ntrac,  nwiso,  ngl
!---wiso-code-end

    INTEGER :: nlevp1
    INTEGER :: dim1(3), dim1p(3)
    INTEGER :: dim2(2), dim2p(2)
    INTEGER :: dim3(4), dim3p(4)
    INTEGER :: dim4(4), dim4p(4) !---wiso-code

    ! construct the g1b table
    !
    ! all information specific to this table is set in this subroutine

    nlevp1 = nlev  + 1

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/ lnlon, lnlev, lngl  /)
    dim1  = (/  nlon,  nlev,  ngl  /)

    dim2p = (/ lnlon, lngl  /)
    dim2  = (/  nlon,  ngl  /)

    dim3p = (/ lnlon, lnlev,  lntrac, lngl /)
    dim3  = (/  nlon,  nlev,   ntrac,  ngl /)

!---wiso-code
    dim4p = (/ lnlon, lnlev,  lnwiso, lngl /)
    dim4  = (/  nlon,  nlev,   nwiso,  ngl /)
!---wiso-code-end

    CALL add_stream_element (g1b, 'vof',   vof,   dim1p, dim1)
    CALL add_stream_element (g1b, 'df',    df,    dim1p, dim1)
    CALL add_stream_element (g1b, 'tf',    tf,    dim1p, dim1)
    CALL add_stream_element (g1b, 'alpsf', alpsf, dim2p, dim2)
    CALL add_stream_element (g1b, 'qf',    qf,    dim1p, dim1)
    CALL add_stream_element (g1b, 'xlf',   xlf,   dim1p, dim1)
    CALL add_stream_element (g1b, 'xif',   xif,   dim1p, dim1)
!---wiso-code
    CALL add_stream_element (g1b, 'wisoqf', wisoqf,dim4p, dim4)
    CALL add_stream_element (g1b, 'wisoxlf',wisoxlf,dim4p,dim4)
    CALL add_stream_element (g1b, 'wisoxif',wisoxif,dim4p,dim4)
!---wiso-code-end
    CALL add_stream_element (g1b, 'xtf',   xtf,   dim3p, dim3)

  END SUBROUTINE construct_g1b

  SUBROUTINE destruct_g1b

    CALL delete_stream (g1b)

  END SUBROUTINE destruct_g1b

END MODULE mo_memory_g1b
