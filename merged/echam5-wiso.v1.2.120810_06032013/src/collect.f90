SUBROUTINE collect (kproma                                             &
            , pahflw,  pahfsw,  pahfice                                &
            , ptrflw,  psoflw                                          &
            , pqres,   pevapw,  pevapi                                 &
            , pustrw,  pvstrw,  pustri,  pvstri                        &
            , palake,  pslf,    pseaice                                &
            , pwind10w                                                 &
            , pawhea,  pawsol,  pawfre,  pawust                        &
            , pawvst,  paicon,  paiqre,  paifre                        &
            , paiust,  paivst,  pawsta                                 &
            , prsfc,   pssfc,   prsfl,   pssfl          )
!
!  ---------------------------------------------------------------------
!
!  Collects surface values for input to the ocean model
!
!  *collect* is called from *physc*
!
!     Authors:
!
!     R. Voss, DKRZ,  August 1997, origianl source
!     R. Voss, U. Schlese, DKRZ, December 1999, modified for echam5
!     S. Legutke, DKRZ, July 2000, modifications for coupling with C-HOPE
!     S. Legutke, DKRZ, Nov  2000, exchange of 4 fluxes (#ifdef PFLUXES4)
!                                  for ocean models without sea ice.
!     I. Kirchner, MPI Hamburg, December 2000, time control
!     S. Legutke,  MPI,M&D, Jan  2002, time control for coupling revisited;
!                  exchange-field accumulation with timestep length weights
!     U. Schlese, M. Esch MPI, Jan 2003, 10m wind instead of ustar3
!     S. Legutke,  MPI,M&D, Apr  2003, coupling revisited;
!        - #ifdef cpl_hope   for coupling with HOPE
!        - #ifdef cpl_mpiom for coupling with MPI-OM
!                           
!     Method:
!
!     Fluxes calculated by ECHAM are accumulated and stored in arrays
!     which serve to transfer these data to *OASIS*.
!
USE mo_kind,         ONLY:  dp
USE mo_constants,    ONLY:  rhoh2o, alf
USE mo_time_control, ONLY:  l_putocean, delta_time
USE mo_couple,       ONLY:  couple_a2o_time
!
IMPLICIT NONE
!
!  scalar arguments
!
  INTEGER, INTENT (IN) :: kproma
!
! array arguments
!
  REAL(dp) ::  pahflw(kproma),  pahfsw(kproma)                         &
            , ptrflw(kproma),  psoflw(kproma)                          &
            , pahfice(kproma), pqres(kproma)                           &
            , pevapw(kproma),  pevapi(kproma)                          &
            , pustrw(kproma),  pvstrw(kproma)                          &
            , pustri(kproma),  pvstri(kproma)                          &
            , palake(kproma),  pslf(kproma),    pseaice(kproma)        &
            , pwind10w(kproma)                                         &
            , pawhea(kproma),  pawsol(kproma),  pawfre(kproma)         &
            , pawust(kproma),  pawvst(kproma),  pawsta(kproma)         &
            , paicon(kproma),  paiqre(kproma),  paifre(kproma)         &
            , paiust(kproma),  paivst(kproma)                          &
            , prsfc(kproma),   pssfc(kproma),   prsfl(kproma)          &
            , pssfl(kproma)
!
! local scalars
!
  INTEGER :: jl

  REAL(dp) ::  zzf1, zzf2, zrcouple

!     accumulate variables for coupling
!
   DO jl=1,kproma
      IF (pslf(jl).LT.1._dp) THEN
        zzf1=1._dp-pseaice(jl)
        zzf2=      pseaice(jl)

#if defined  __cpl_mpiom
          pawhea(jl) = pawhea(jl)+( pahflw(jl)+pahfsw(jl)              &
                                   +ptrflw(jl)+psoflw(jl)              &
                                  -(pssfl(jl)+pssfc(jl))*alf*zzf1)     &
                                  *delta_time
          pawust(jl) = pawust(jl)+pustrw(jl)*delta_time
          pawvst(jl) = pawvst(jl)+pvstrw(jl)*delta_time
#endif
#ifdef __cpl_hope
          pawhea(jl) = pawhea(jl)+( pahflw(jl)+pahfsw(jl)              &
                                   +ptrflw(jl)+psoflw(jl)              &
                                   -(pssfl(jl)+pssfc(jl))*alf)         &
                                   *zzf1*delta_time
          pawust(jl) = pawust(jl)+pustrw(jl)*zzf1*delta_time
          pawvst(jl) = pawvst(jl)+pvstrw(jl)*zzf1*delta_time
#endif
          pawsol(jl) = pawsol(jl)+psoflw(jl)*delta_time
          pawsta(jl) = pawsta(jl)+pwind10w(jl)*delta_time

#if defined __cpl_mpiom
          paicon(jl) = paicon(jl)+ pahfice(jl)*delta_time
          paiqre(jl) = paiqre(jl)+ pqres(jl)*delta_time
          paiust(jl) = paiust(jl)+ pustri(jl)*delta_time
          paivst(jl) = paivst(jl)+ pvstri(jl)*delta_time
#endif
#ifdef __cpl_hope
          paicon(jl) = paicon(jl)+ pahfice(jl)*zzf2*delta_time
          paiqre(jl) = paiqre(jl)+ pqres(jl)  *zzf2*delta_time
#ifdef __cpl_fluxes4
          paiust(jl) = paiust(jl)+pustri(jl)*zzf2*delta_time
          paivst(jl) = paivst(jl)+pvstri(jl)*zzf2*delta_time
#else
          paiust(jl) = paiust(jl)+ pustri(jl)*delta_time
          paivst(jl) = paivst(jl)+ pvstri(jl)*delta_time
#endif
#endif
        pawfre(jl) = pawfre(jl)+(prsfl(jl)+prsfc(jl))*delta_time       &
                               +(pssfl(jl)+pssfc(jl)+pevapw(jl))       &
                               *zzf1*delta_time
        paifre(jl) = paifre(jl)+(pssfc(jl)+pssfl(jl)+pevapi(jl))       &
                               *zzf2*delta_time
      END IF
    END DO
!
!    prepare coupling fields before transfer:
!            average and convert freshwater fluxes to [m/s]
!
   IF (l_putocean) THEN

     zrcouple = 1.0_dp/couple_a2o_time

     DO jl = 1,kproma
#ifdef __cpl_fluxes4
       pawhea(jl) = (pawhea(jl)+paicon(jl)+paiqre(jl))*zrcouple
       pawfre(jl) = (pawfre(jl)+paifre(jl))*zrcouple/rhoh2o
       pawust(jl) = (pawust(jl)+paiust(jl))*zrcouple
       pawvst(jl) = (pawvst(jl)+paivst(jl))*zrcouple
#else
       pawhea(jl) = pawhea(jl)*zrcouple
       pawfre(jl) = pawfre(jl)*zrcouple/rhoh2o
       pawust(jl) = pawust(jl)*zrcouple
       pawvst(jl) = pawvst(jl)*zrcouple
       paicon(jl) = paicon(jl)*zrcouple
       paiqre(jl) = paiqre(jl)*zrcouple
       paifre(jl) = paifre(jl)*zrcouple/rhoh2o
       paiust(jl) = paiust(jl)*zrcouple
       paivst(jl) = paivst(jl)*zrcouple
#endif
       pawsol(jl) = pawsol(jl)*zrcouple
       pawsta(jl) = pawsta(jl)*zrcouple
     END DO

   END IF

   RETURN
   END
