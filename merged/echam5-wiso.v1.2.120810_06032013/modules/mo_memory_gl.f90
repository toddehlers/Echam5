MODULE mo_memory_gl
!------------------------------------------------------------------------------
  !
  ! Modules used
  !

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: new_stream, delete_stream, add_stream_element, &
                            default_stream_setting, add_stream_reference
  USE mo_netCDF,      ONLY: max_dim_name
  USE mo_tracdef,     ONLY: trlist
  USE mo_filename,    ONLY: trac_filetype
!---wiso-code
  USE mo_wiso,        ONLY: wisoq_names, wisoxl_names, wisoxi_names
!---wiso-code-end

  IMPLICIT NONE
!------------------------------------------------------------------------------
  !
  ! Public entities
  !

  PRIVATE

  PUBLIC :: construct_gl ! construct the gl table
  PUBLIC :: destruct_gl  ! destruct  the gl table

  PUBLIC :: gl           ! the gl table
  PUBLIC :: tracer       ! the tracer table
!------------------------------------------------------------------------------

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: q(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xl(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xi(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xt(:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: lammp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: phimp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sigmp(:,:,:)
!---wiso-code
  REAL(dp), POINTER, PUBLIC :: wisoq  (:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxl (:,:,:,:)
  REAL(dp), POINTER, PUBLIC :: wisoxi (:,:,:,:)
!---wiso-code-end


  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: gl
  TYPE (t_stream), POINTER :: tracer

  ! private storage for tracer fields

!---wiso-code
  REAL(dp), POINTER :: pwisoq  (:,:,:,:,:)
  REAL(dp), POINTER :: pwisoxl (:,:,:,:,:)
  REAL(dp), POINTER :: pwisoxi (:,:,:,:,:)
!---wiso-code-end
  REAL(dp), POINTER :: pxt (:,:,:,:,:)

CONTAINS
!------------------------------------------------------------------------------
!---wiso-code
  SUBROUTINE construct_gl (lnlon, lnlev, lntrac, lnwiso, lngl, &
                            nlon,  nlev,  ntrac,  nwiso,  ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lntrac, lnwiso, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ntrac,  nwiso,  ngl
!---wiso-code-end

    INTEGER                      :: dim1(3), dim1p(3)
    INTEGER                      :: dim2(4), dim2p(4)
    INTEGER                      :: dim3(4), dim3p(4)  !---wiso-code
    CHARACTER (len=max_dim_name) :: dim1n(3), dim2n(4)
    CHARACTER (len=max_dim_name) :: dim3n(4)           !---wiso-code
    INTEGER                      :: i
    INTEGER                      :: zi                 !---wiso-code
    REAL(dp), POINTER            :: p3(:,:,:), p4(:,:,:,:)

    ! construct the gl table
    !
    ! all information specific to this table is set in this subroutine


    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/  lnlon,  lnlev, lngl  /)
    dim1  = (/   nlon,   nlev,  ngl  /)
    dim1n = (/   "lon ","lev ","lat "/)

    dim2p = (/ lnlon,    lnlev,   lntrac,  lngl    /)
    dim2  = (/  nlon,     nlev,    ntrac,   ngl    /)
    dim2n = (/  "lon   ","lev   ","ntrac ","lat   "/)

!---wiso-code
    dim3p = (/ lnlon,    lnlev,   lnwiso,  lngl    /)
    dim3  = (/  nlon,     nlev,    nwiso,   ngl    /)
    dim3n = (/  "lon   ","lev   ","nwiso ","lat   "/)
!---wiso-code-end

    CALL default_stream_setting (gl ,dimnames=dim1n ,lrerun=.TRUE. ,table=128)

    CALL add_stream_element (gl,'q',    q,    code=133 ,longname='specific humidity'        ,units='kg/kg')
    CALL add_stream_element (gl,'xl',   xl,   code=153 ,longname='cloud water'              ,units='kg/kg')
    CALL add_stream_element (gl,'xi',   xi,   code=154 ,longname='cloud ice'                ,units='kg/kg')
    CALL add_stream_element (gl,'lammp',lammp,lpost=.FALSE., &
         lrerun=.FALSE.,contnorest=.TRUE.) 
    CALL add_stream_element (gl,'phimp',phimp,lpost=.FALSE., &
         lrerun=.FALSE.,contnorest=.TRUE.) 
    CALL add_stream_element (gl,'sigmp',sigmp,lpost=.FALSE., &
         lrerun=.FALSE.,contnorest=.TRUE.) 

    !
    ! Special handling for tracers
    !

!---wiso-code

    IF (nwiso > 0) THEN
      !
      ! Allocate a 5d-array (with dummy index 5) for tracers. This array 
      ! is referenced by the 4-d array 'XT' and by the 3-d arrays of 
      ! individual tracers.
      !
      ALLOCATE (pwisoq  (lnlon, lnlev, lnwiso, lngl, 1))
      ALLOCATE (pwisoxl (lnlon, lnlev, lnwiso, lngl, 1))
      ALLOCATE (pwisoxi (lnlon, lnlev, lnwiso, lngl, 1))
      !
      ! Set meta information on water isotope arrays.
      ! Set restart flag, obtain reference to memory info entry.
      !
      p4 => pwisoq(:,:,:,:,1)
      CALL add_stream_element (gl, 'wisoq', wisoq, ldims=dim3p, gdims=dim3,dimnames=dim3n, &
                               longname='specific humidity - water isotopes', p4=p4,lrerun=.FALSE.,lpost=.FALSE.,table=128)
      !
      p4 => pwisoxl(:,:,:,:,1)
      CALL add_stream_element (gl, 'wisoxl', wisoxl, ldims=dim3p, gdims=dim3,dimnames=dim3n,               &
                               longname='cloud water - water isotopes', p4=p4,lrerun=.TRUE.,table=128)
      !
      p4 => pwisoxi(:,:,:,:,1)
      CALL add_stream_element (gl, 'wisoxi', wisoxi, ldims=dim3p, gdims=dim3,dimnames=dim3n,               &
                               longname='cloud ice - water isotopes', p4=p4,lrerun=.TRUE.,table=128)
      !
      ! provide additional meta-information for individual tracers.
      !
      zi = 240     ! Water isotope tracer fields start at code #241
      !
      DO i = 1, lnwiso
        p4 => pwisoq(:,:,i,:,:)
        zi = zi+1
        CALL add_stream_element (gl, wisoq_names(i), p3, ldims=dim1p, gdims=dim1, dimnames=dim1n,          &
                                 longname='specific humidity - water isotopes',units='kg/kg',lrerun=.TRUE.,  &
                                 lpost=.TRUE.,code = zi,p4=p4,table=128)
        !
        p4 => pwisoxl(:,:,i,:,:)
        zi = zi+1
        CALL add_stream_element (gl, wisoxl_names(i), p3, ldims=dim1p, gdims=dim1, dimnames=dim1n,         &
                                 longname='cloud water - water isotopes',units='kg/kg',lrerun=.TRUE.,        &
                                 lpost=.TRUE.,code = zi,p4=p4,table=128)
        !
        p4 => pwisoxi(:,:,i,:,:)
        zi = zi+1
        CALL add_stream_element (gl, wisoxi_names(i), p3, ldims=dim1p, gdims=dim1, dimnames=dim1n,         &
                                 longname='cloud ice - water isotopes',units='kg/kg',lrerun=.TRUE.,          &
                                 lpost=.TRUE.,code = zi,p4=p4,table=128)
      END DO
    ELSE
      NULLIFY (pwisoq)
      NULLIFY (pwisoxl)
      NULLIFY (pwisoxi)
      CALL add_stream_element (gl,'wisoq', wisoq, dim3p, dim3, &
                               lpost   =.FALSE.,         &
                               lrerun  =.FALSE.)
      CALL add_stream_element (gl,'wisoxl', wisoxl, dim3p, dim3, &
                               lpost   =.FALSE.,         &
                               lrerun  =.FALSE.)
      CALL add_stream_element (gl,'wisoxi', wisoxi, dim3p, dim3, &
                               lpost   =.FALSE.,         &
                               lrerun  =.FALSE.)
    END IF

!---wiso-code-end

    IF (ntrac > 0) THEN
      !
      ! Allocate a 5d-array (with dummy index 5) for tracers. This array 
      ! is referenced by the 4-d array 'XT' and by the 3-d arrays of 
      ! individual tracers.
      !
      ALLOCATE (pxt(lnlon, lnlev, lntrac, lngl, 1))
      !
      ! Set meta information on XT array.
      ! Set restart flag, obtain reference to memory info entry.
      !
      p4 => pxt(:,:,:,:,1)
      CALL add_stream_element (gl, 'xt', xt, dim2p, dim2,       &
                               dimnames   = dim2n,              &
                               lrerun     = trlist% oldrestart, &
                               contnorest = .TRUE.,             &   
                               mem_info   = trlist% mixt,       &
                               lpost      = .FALSE.,            &
                               p4=p4)
      !
      ! setup a new output stream
      !
      CALL new_stream (tracer ,'tracer',trac_filetype)
      !
      ! add entries for geopotential, log surface pressure, grid-box area
      !
      CALL add_stream_reference (tracer, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (tracer, 'lsp'     ,'sp'    ,lpost=.TRUE.)
      CALL add_stream_reference (tracer, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
      CALL add_stream_reference (tracer, 'gboxarea','geoloc',lpost=.TRUE.)
      !
      ! provide additional meta-information for individual tracers.
      !
      DO i = 1, lntrac
        p4 => pxt(:,:,i,:,:)

        CALL add_stream_element (tracer, trlist% ti(i)% fullname, p3,      &
                                 dim1p, dim1,                              &
                                 dimnames   = dim1n,                       &
                                 units      = trlist% ti(i)% units,        &
                                 longname   = trlist% ti(i)% longname,     &
                                 lrerun     = trlist% ti(i)% nrerun==1,    &
                                 contnorest = .TRUE.,                      &
                                 mem_info   = trlist% mi(i)% xt,           &
                                 lpost      = trlist% ti(i)% nwrite==1,    &
                                 code       = trlist% ti(i)% code,         &
                                 table      = trlist% ti(i)% table,        &
                                 bits       = trlist% ti(i)% gribbits,     &
                                 tracidx    = i,                           &
                                 p4         = p4)
      END DO
    ELSE
      NULLIFY (tracer)
      NULLIFY (pxt)
      CALL add_stream_element (gl,'xt', xt, dim2p, dim2, &
                               lpost   =.FALSE.,         &
                               lrerun  =.FALSE.,         &
                               mem_info=trlist% mixt)
    END IF

  END SUBROUTINE construct_gl
!------------------------------------------------------------------------------
  SUBROUTINE destruct_gl

    CALL delete_stream (gl)
    CALL delete_stream (tracer)

!---wiso-code
    IF(ASSOCIATED (pwisoq)) DEALLOCATE (pwisoq)
    IF(ASSOCIATED (pwisoxl)) DEALLOCATE (pwisoxl)
    IF(ASSOCIATED (pwisoxi)) DEALLOCATE (pwisoxi)
!---wiso-code-end
    IF(ASSOCIATED (pxt)) DEALLOCATE (pxt)

  END SUBROUTINE destruct_gl
!------------------------------------------------------------------------------
END MODULE mo_memory_gl
