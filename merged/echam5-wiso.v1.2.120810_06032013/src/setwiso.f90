  SUBROUTINE setwiso
  !
  ! preset water isotope physical constants according to namelist 
  !
  ! M. Werner, AWI Bremerhaven, 2010
  !
  ! This subroutine is called by initialize.f90
  !

    USE mo_doctor,        ONLY: nout
    USE mo_kind,          ONLY: dp
    USE mo_mpi,           ONLY: p_io, p_parallel, p_parallel_io, p_bcast
    USE mo_namelist,      ONLY: position_nml, nnml, POSITIONED 
    USE mo_constants,     ONLY: rhoh2o
    USE mo_exception,     ONLY: finish
    
    USE mo_wiso,          ONLY: wisoctl, nwiso, nwisotyp, nwisomax,                      &
                                wisoq_names,   wisoxl_names,   wisoxi_names,             &
                                wisoqm1_names, wisoxlm1_names, wisoxim1_names,           &
                                twisoatm, tnat, toce, twisorhoh2o,                       &
                                twisosur1, twisosur2,                                    &
                                talphal1, talphal2, talphal3,                            &
                                talphas1, talphas2, talphas3,                            &
                                tkinsl, tkinfa1, tkinfa2, tdifrel
                                

    IMPLICIT NONE

    CHARACTER(LEN=7) :: q_name,   xl_name,   xi_name,   &
                        qm1_name, xlm1_name, xim1_name
    
    REAL(dp):: &
            zwisoatm,znat,zoce,zwisorhoh2o,             &
            zwisosur1,zwisosur2,                        &
            zalphal1,zalphal2,zalphal3,                 &
            zalphas1,zalphas2,zalphas3,                 &
            zkinsl,zkinfa1,zkinfa2,zdifrel
    
    INTEGER     :: ierr        ! error return value from position_nml
    INTEGER     :: jt


    ! set default values of NAMELIST variables
    nwiso      = 0

    nwisotyp(1)=1     ! H2-16O
    nwisotyp(2)=2     ! H2-18O
    nwisotyp(3)=3     ! HDO

    ! read namelist wisoctl
    IF (p_parallel_io) THEN
       CALL position_nml ('WISOCTL', status=ierr)
       SELECT CASE (ierr)
          CASE (POSITIONED)
             READ (nnml, wisoctl)
       END SELECT
    ENDIF
    IF (p_parallel) THEN
     CALL p_bcast (nwiso, p_io)
    ENDIF

    ! check if nwiso .LE. nwisomax
    IF (p_parallel_io) THEN
      IF (nwiso .GT. nwisomax) THEN
        CALL finish('setwiso','Specified number of tracers NWISO too large (> 3)!')
      ENDIF
    ENDIF

    ! check if nwisotyp is known
    IF (p_parallel_io) THEN
      DO jt=1,nwiso
        IF (nwisotyp(jt).LT.1 .OR. nwisotyp(jt).GT.3) THEN
          CALL finish('setwiso','Specified type of tracer NWISOTYP unknown!')
        ENDIF
      END DO
    ENDIF

    ! set isotope-dependent coefficients
    ! (tracer types: 1 = H2-16O, 2 = H2-18O, 3 = HD-16O)
    DO jt=1,nwiso
       IF (nwisotyp(jt).eq.1) THEN ! O16
          q_name   ='q16o   '
          qm1_name ='q16om1 '
          xl_name  ='xl16o  '
          xlm1_name='xl16om1'
          xi_name  ='xi16o  '
          xim1_name='xi16om1'
          zwisoatm=0._dp
          znat=1._dp
          zoce=1._dp
          zwisosur1=0._dp
          zwisosur2=0._dp
          zwisorhoh2o=rhoh2o
          zdifrel=1._dp
          zalphal1=0._dp
          zalphal2=0._dp
          zalphal3=0._dp
          zalphas1=0._dp
          zalphas2=0._dp
          zalphas3=0._dp
          zkinsl=0._dp
          zkinfa1=0._dp
          zkinfa2=0._dp
       ELSEIF (nwisotyp(jt).eq.2) THEN ! O18
          q_name   ='q18o   '
          qm1_name ='q18om1 '
          xl_name  ='xl18o  '
          xlm1_name='xl18om1'
          xi_name  ='xi18o  '
          xim1_name='xi18om1'
          zwisoatm=-20._dp/1000._dp
          znat=2005.2e-6_dp
          zoce=znat
          zwisosur1=0.69_dp/1000._dp
          zwisosur2=13.6_dp/1000._dp
          zwisorhoh2o=rhoh2o
          zdifrel=1._dp/0.9723_dp
          zalphal1=1137._dp
          zalphal2=-0.4156_dp
          zalphal3=-2.0667e-3_dp

! default fractionation factors for vapour-ice
          zalphas1=0._dp
          zalphas2=11.839_dp
          zalphas3=-0.028244_dp
! new fractionation factors for vapour-ice by M. Ellehoj
!          zalphas1=8312.5_dp
!          zalphas2=-49.192_dp
!          zalphas3=0.0831_dp

          zkinsl=0.006_dp
          zkinfa1=0.000285_dp
          zkinfa2=0.00082_dp
       ELSEIF (nwisotyp(jt).eq.3) THEN ! HDO
          q_name   ='qhdo   '
          qm1_name ='qhdom1 '
          xl_name  ='xlhdo  '
          xlm1_name='xlhdom1'
          xi_name  ='xihdo  '
          xim1_name='xihdom1'
          zwisoatm=-150._dp/1000._dp
          znat=155.76e-6_dp
          zoce=znat
          zwisosur1=5.6_dp/1000._dp
          zwisosur2=100._dp/1000._dp
          zwisorhoh2o=rhoh2o
          zdifrel=1._dp/0.9755_dp
          zalphal1=24844._dp
          zalphal2=-76.248_dp
          zalphal3=52.612e-3_dp

! default fractionation factors for vapour-ice
          zalphas1=16288._dp
          zalphas2=0._dp
          zalphas3=-0.0934_dp
! new fractionation factors for vapour-ice by M. Ellehoj
!          zalphas1=48888._dp
!          zalphas2=-203.10_dp
!          zalphas3=0.2133_dp

          zkinsl=0.00528_dp
          zkinfa1=0.0002508_dp
          zkinfa2=0.0007216_dp
       ENDIF
       ! set variable names
       wisoq_names(jt)    = q_name
       wisoqm1_names(jt)  = qm1_name
       wisoxl_names(jt)   = xl_name
       wisoxlm1_names(jt) = xlm1_name
       wisoxi_names(jt)   = xi_name
       wisoxim1_names(jt) = xim1_name
       ! initialize deviation from SMOW in the atmosphere (if not set via wiso namelist)
       !  (default: 18o:-20 permill, hdo:-150 permill)
!       if (twisoatm(jt).EQ. -9.e9) twisoatm(jt)=zwisoatm
       twisoatm(jt)=zwisoatm
       ! tnat = natural isotope distribution to V-SMOW,1:T,2:D,3:O18,4:H2O^16
       tnat(jt)=znat
       ! toce = mean concentration in ocean surface water
       toce(jt)=zoce
       ! txtsur = initial surface concentration constants
       twisosur1(jt)=zwisosur1
       twisosur2(jt)=zwisosur2
       ! wisorhoh2o=density of the isotopic water
       twisorhoh2o(jt)=zwisorhoh2o
       ! the relation between the diffusivity of "water" and the water isotope
       tdifrel(jt)=zdifrel
       ! some constants for the calculation of the fractionation factor
       talphal1(jt)=zalphal1
       talphal2(jt)=zalphal2
       talphal3(jt)=zalphal3
       ! some constants for the calculation of the fractionation factor over ice
       talphas1(jt)=zalphas1
       talphas2(jt)=zalphas2
       talphas3(jt)=zalphas3
       ! some constants for the kinetic fractionation over the ocean
       tkinsl(jt)=zkinsl
       tkinfa1(jt)=zkinfa1
       tkinfa2(jt)=zkinfa2
    ENDDO
  
    ! write namelist wisoctl to stdout    
    IF (p_parallel_io) THEN
      IF (nwiso .gt. 0) WRITE(nout,wisoctl)
    END IF
       
  END SUBROUTINE setwiso
