SUBROUTINE surf ( kproma, kbdim, klev, klevp1       &
! - 1D from mo_memory_g3b
          , ptsl,       ptslm,     ptslm1           &
          , pws,        pwl,       pwsmx            &
          , psn,        psnmel,    pgld             &
          , psnc,       pu10,      pv10             &
          , prunoff,    progl,     pdrain           &
          , papmegl,    psnacl,    porostd          &
          , prgcgn,     psodif,    pslm             &
          , pgrndcapc,  pgrndhflx, pgrndflux        &
! - 2D from mo_memory_g3b (soil variables)
          , ptsoil,     pgrndd,    pgrndc           &
          , ptm1,       pqm1,      ptte             &
          , paphm1                                  &
! - variables internal to physics
          , pcvs,       pcvw,      pwlmx            &
          , pevapl,     pevapot                     &
          , prsfl,      prsfc                       &
          , pssfl,      pssfc                       &
          , pros_hd,    pdrain_hd                   &
          , palac                                   &
          , lpland,     lpglac                      &
!---wiso-code
! - 1D from mo_memory_g3b - water isotopes
          , pwisows,    pwisowl,   pwisosn          &
          , pwisosnmel,pwisogld,   pwisosnc         &
          , pwisorunoff,pwisorogl, pwisodrain       &
          , pwisoapmegl,pwisosnacl                  &
! - variables internal to physics - water isotopes
          , pwisoevapl, pwisoevapot                 &
          , pwisorsfl,  pwisorsfc                   &
          , pwisossfl,  pwisossfc                   &
          , pwisoros_hd,pwisodrain_hd               &
          , pwisoalac)
!---wiso-code-end
!
!     ------------------------------------------------------------------
! ptm1          Air temp [K] at timestep t-dt (filtered)
! pu10          Wind (u-component) at 10m height [m/s] ... from 'vdiff'
! pv10          Wind (v-component) at 10m height [m/s] ... from 'vdiff'
! ptsl          Sfc temp [K] at timestep t+dt (unfiltered)
! ptslm         Sfc temp [K] at timestep t (unfiltered)
! ptslm1        Sfc temp [K] at timestep t-dt (filtered)
! pws           Soil water content [m]
! pwl           Water content [m] in skin reservoir
!               (vegetation and bare soil)
! pwsmx         Water holding capacity [m] of the soil
! pslm          Land sea mask
! psn           Snow depth [m water equivalent] at the ground
! psnc          Snow depth [m water equivalent] at the canopy
! psnmel        Snow melt [kg/m**2] (accumulated for diagnostics)
! pgld          Glacier depth (including snow) [m water equivalent]
! prunoff       Total runoff [kg/m**2] at non-glacier points (accumul.)
! progl         Glacier runoff (rain+snow/ice melt) [kg/m**2] (accumul.)
! pdrain        Drainage at non-glacier points [kg/m**2] (accumul.)
! papmegl       Precip-Evap at glacier points [kg/m**2] (accumul.)
! psnacl        Snow budget at non-glacier points [kg/m**2] (accumul.)
! porostd       Subgrid standard diviation [m] used in runoff scheme
! prgcgn        Volumetric heat capacity of the soil [j/m**3/K]
! psodif        Thermal diffusivity of the soil [m**2/s]
! pgrndcapc     Heat capacity of the uppermost soil layer [j/m**2/K]
! pgrndhflx     Soil heat flux at the surface [W/m**2]
! ptsoil        Temperature [K] in the five soil layers
! pgrndc,pgrndd Coefficients used in the ptsoil calculation (*soiltemp*)
! pcvs          Fractional snow cover (function of psn in *physc*)
! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
!               of vegetation index and leaf area index)
! pevapl        Total evaporation, including sublimation [kg/m**2/s]
! pevapot       Potential evaporation/sublimation [kg/m**2/s]
! prsfl         Large scale rainfall [kg/m**2/s]
! prsfc         Convective rainfall [kg/m**2/s]
! pssfl         Large scale snowfall [kg/m**2/s]
! pssfc         Convective snowfall [kg/m**2/s]
! pros_hd       Runoff for HD-Model (does NOT include drainage) [m]
! pdrain_hd     Drainage for HD-Model [m]
! palac         Precipitation minus sublimation at glacier points
! lpland        Logical land mask (including glaciers)
! lpglac        Logical glacier mask
!
!---wiso-code
!
! pwisows	tracer for Soil water content [m]
! pwisowl	tracer for Water content [m] in skin reservoir
!               (vegetation and bare soil)
! pwisosn	tracer for Snow depth [m water equivalent] at the ground
! pwisosnc	tracer for Snow depth [m water equivalent] at the canopy
! pwisosnmel	tracer for Snow melt [kg/m**2] (accumulated for diagnostics)
! pwisogld	tracer for Glacier depth (including snow) [m water equivalent]
! pwisorunoff	tracer for Total runoff [kg/m**2] at non-glacier points (accumul.)
! pwisorogl	tracer for Glacier runoff (rain+snow/ice melt) [kg/m**2] (accumul.)
! pwisodrain	tracer for Drainage at non-glacier points [kg/m**2] (accumul.)
! pwisoapmegl	tracer for Precip-Evap at glacier points [kg/m**2] (accumul.)
! pwisosnacl	tracer for Snow budget at non-glacier points [kg/m**2] (accumul.)
! pwisoevapl	tracer for Total evaporation, including sublimation [kg/m**2/s]
! pwisoevapot	tracer for Potential evaporation/sublimation [kg/m**2/s]
! pwisorsfl	tracer for Large scale rainfall [kg/m**2/s]
! pwisorsfc	tracer for Convective rainfall [kg/m**2/s]
! pwisossfl	tracer for Large scale snowfall [kg/m**2/s]
! pwisossfc	tracer for Convective snowfall [kg/m**2/s]
! pwisoros_hd	tracer for Runoff for HD-Model (does NOT include drainage) [m]
! pwisodrain_hd	tracer for Drainage for HD-Model [m]
! pwisoalac	tracer for Precipitation minus sublimation at glacier points
!
!---wiso-code-end
!
! The following local variables represent the respective fluxes
! integrated over one timestep (delta_time) and divided by the density of
! water (rhoh2o). Units are m water equivalent.
!
! zraind        Total rain
! zsnowd        Total snow
! zevttd        Total evaporation
! zevsnd        Sublimation from snow
! zevwld        Evaporation from the skin reservoir
! zevwsd        Evaporation from the soil and from the skin reservoir
! zros          Total runoff (including drainage) at non-glacier points
! zdrain        Drainage at non-glacier points
! zrogl         Runoff at glacier points (rain and melt, but no calving)
! zsnmel        Snow/ice melt at land and glacier points
! zsncmelt      Snow melt in the canopy
! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
! zmlres        Residual melt water available for infiltration into the
!               non-frozen soil after being intercepted by the
!               skin reservoir
!
!---wiso-code
!
! zwisoraind	tracer for Total rain
! zwisosnowd	tracer for Total snow
! zwisoevttd	tracer for Total evaporation
! zwisoevsnd	tracer for Sublimation from snow
! zwisoevwld	tracer for Evaporation from the skin reservoir
! zwisoevwsd	tracer for Evaporation from the soil and from the skin reservoir
! zwisoros      tracer for Total runoff (including drainage) at non-glacier points
! zwisodrain	tracer for Drainage at non-glacier points
! zwisorogl		tracer for Runoff at glacier points (rain and melt, but no calving)
! zwisosnmel	tracer for Snow/ice melt at land and glacier points
! zwisosncmelt	tracer for Snow melt in the canopy
! zwisosn		tracer for Snow budget at non-glacier points (snowfall-subl-melt)
! zwisomlres	tracer for Residual melt water available for infiltration into the
!				non-frozen soil after being intercepted by the skin reservoir
!
! zgld_tmp		Temporary file for Glacier depth (including snow) [m water equivalent]
! zsn_tmp		Temporary file for Snow depth [m water equivalent] at the ground
! zb1_tmp		Temporary file for zb1
!
!---wiso-code-end
!
!     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
!              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
!              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
!              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
!              DURING TIMESTEP DELTA_TIME.
!
!     J.F.GELEYN     E.C.M.W.F.     08/06/82.
!     MODIFIED BY
!     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
!     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
!     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
!                                 COUPLING BETWEEN LAND SURFACE
!                                 AND ATMOSPHERE.
!     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
!     MODIFIED BY M. ESCH        MPI - APR  1999
!     MODIFIED BY E. ROECKNER    MPI - JAN  2001
!     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
!     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
!                                                for snow changed
!     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
!
!     MODIFICATION
!
!     PURPOSE
!
!     INTERFACE.
!
!          *SURF* IS CALLED FROM *PHYSC*.
!
!     METHOD.
!
!     EXTERNALS.
!
!          NONE.
!
!     REFERENCE.
!
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
!
USE mo_kind,           ONLY: dp
USE mo_control,        ONLY: ngl
USE mo_param_switches, ONLY: lsurf
USE mo_physc2,         ONLY: jpgrnd, cwlmax
USE mo_constants,      ONLY: rhoh2o, g, cpd, vtmpc2, alf, tmelt
USE mo_vegetation,     ONLY: cvinter
USE mo_time_control,   ONLY: lstart, delta_time
USE mo_semi_impl,      ONLY: eps

!---wiso-code

USE mo_wiso,           ONLY: nwiso, twisorhoh2o, tnat, calc_delta

!---wiso-code-end
!
 IMPLICIT NONE
 
  INTEGER :: kproma, kbdim, jl, klev, klevp1
  REAL(dp) ::                                                             &
       ptsl(kbdim),          ptslm(kbdim),         ptslm1(kbdim)         &
     , pws(kbdim),           pwl(kbdim),           pwsmx(kbdim)          &
     , psn(kbdim),           psnmel(kbdim),        pgld(kbdim)           &
     , psnc(kbdim),          pu10(kbdim),          pv10(kbdim)           &
     , prunoff(kbdim),       progl(kbdim),         pdrain(kbdim)         &
     , papmegl(kbdim),       psnacl(kbdim),        porostd(kbdim)        &
     , prgcgn(kbdim),        psodif(kbdim),        pslm(kbdim)           &
     , pgrndcapc(kbdim),     pgrndhflx(kbdim),     pgrndflux(kbdim)
  REAL(dp) ::                                                             &
       ptsoil(kbdim,jpgrnd), pgrndd(kbdim,jpgrnd), pgrndc(kbdim,jpgrnd)  &
     , ptm1(kbdim,klev),     pqm1(kbdim,klev),     ptte(kbdim,klev)      &
     , paphm1(kbdim,klevp1)
  REAL(dp) ::                                                             &
       pcvs(kbdim),          pcvw(kbdim),          pwlmx(kbdim)          &
     , pevapl(kbdim),        pevapot(kbdim)                             &
     , prsfl(kbdim),         prsfc(kbdim)                               &
     , pssfl(kbdim),         pssfc(kbdim)                               &
     , pros_hd(kbdim),       pdrain_hd(kbdim),     palac(kbdim)
  LOGICAL ::                                                          &
       lpland(kbdim),        lpglac(kbdim)
!
!  local arrays
!
  REAL(dp) ::                                                             &
       zraind (kbdim),       zsnowd(kbdim),         zevttd(kbdim)        &
     , zevsnd(kbdim),        zevwsd(kbdim),         zevwld(kbdim)        &
     , zros(kbdim),          zdrain(kbdim),         zrogl(kbdim)         &
     , zsnmel(kbdim),        zsn(kbdim)                                 &
     , zmlres(kbdim),        zsncmelt(kbdim)                             &
     , zdp(kbdim),           zlfdcp(kbdim)
!
!  local scalars
!
  REAL(dp) ::                                                             &
       zorvari, zorvars, zdrmin, zdrmax, zdrexp, zsmelt, zsnmlt,      &
       zmprcp, zwlp, zwptr, zwdtr, zwslim, zconw2, zconw3, zconw4,    &
       zroeff, zbws, zb1, zbm, zconw1, zlyeps, zvol, zinfil, zprfl,   &
       zlysic, zwsup, zc2, zc3, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zdtime, zrcp, zsncfac
!

!---wiso-code

! declaration parameters water isotopes
!
  INTEGER :: jt
  REAL(dp) ::                                                                      &
       pwisows(kbdim,nwiso),    pwisowl(kbdim,nwiso),      pwisosn(kbdim,nwiso)    &
     , pwisosnmel(kbdim,nwiso), pwisogld(kbdim,nwiso),     pwisosnc(kbdim,nwiso)   &
     , pwisorunoff(kbdim,nwiso),pwisorogl(kbdim,nwiso),    pwisodrain(kbdim,nwiso) &
     , pwisoapmegl(kbdim,nwiso),pwisosnacl(kbdim,nwiso)
  REAL(dp) ::                                                                      &
       pwisoevapl(kbdim,nwiso), pwisoevapot(kbdim,nwiso)                           &
     , pwisorsfl(kbdim,nwiso),	pwisorsfc(kbdim,nwiso)                             &
     , pwisossfl(kbdim,nwiso),  pwisossfc(kbdim,nwiso)                             &
     , pwisoros_hd(kbdim,nwiso),pwisodrain_hd(kbdim,nwiso),pwisoalac(kbdim,nwiso)
       
! water isotopes: local arrays  
!
  REAL(dp) ::                                                                      &
       zwisoraind(kbdim,nwiso), zwisosnowd(kbdim,nwiso),   zwisoevttd(kbdim,nwiso) &
     , zwisoevsnd(kbdim,nwiso), zwisoevwsd(kbdim,nwiso),   zwisoevwld(kbdim,nwiso) &
     , zwisoros(kbdim,nwiso),   zwisodrain(kbdim,nwiso),   zwisorogl(kbdim,nwiso)  &
     , zwisosnmel(kbdim,nwiso), zwisosn(kbdim,nwiso)                               &
     , zwisomlres(kbdim,nwiso), zwisosncmelt(kbdim,nwiso)                          &
     , zsnowd_tmp(kbdim),       zmprcp_tmp(kbdim)                                  &
     , zsncp_tmp(kbdim),        zevsnd_tmp(kbdim),         zsubl_tmp(kbdim)        &
     , zsnc_tmp(kbdim)                                                             &
     , zsnc_tmp2(kbdim),        zsncwind_tmp(kbdim)                                &
     , zwl_tmp(kbdim),          zsnmlt_tmp(kbdim)                                  &
     , zraind_tmp(kbdim)                                                           &
     , zwlp_tmp(kbdim),         zevwld_tmp(kbdim),         zevap_tmp(kbdim)        &
     , zgld_tmp(kbdim),	        zsn_tmp(kbdim),            zb1_tmp(kbdim)          &
     , ztsl_tmp(kbdim),         zws_tmp(kbdim),            zprfl_tmp(kbdim)        &
     , zros_tmp(kbdim),         zwsup_tmp(kbdim)

  LOGICAL :: lo_sub(kbdim), lo_evap(kbdim), lo_zros(kbdim)

!  water isotopes: local scalars
!
  REAL(dp) ::                                                                      &
       zwisosnmlt, zwisomprcp, zwisowlp,                                           &
       zwisoinfil, zwisoprfl,  zwisowsup,    zwisosncp,                            &
       zwisosubl,  zwisoevap,  zwisosncwind, zdelta,                               &
       zwisomin, zwisosec

!---wiso-code-end
!
!      Parameters
!
  zdtime = delta_time
  zorvari=100._dp
  zorvars=1000._dp*64._dp/ngl
  zdrmin=0.001_dp/(3600._dp*1000._dp)
  zdrmax=0.1_dp/(3600._dp*1000._dp)
  zdrexp=1.5_dp
  zc2=1.87E5_dp
  zc3=1.56E5_dp
  zsncfac=rhoh2o*g/zdtime

!---wiso-code

  zwisomin= 1.e-12_dp
  zwisosec= 1.e-10_dp
  
!---wiso-code-end

!
!     ------------------------------------------------------------------
!
!*    1.     Convert water fluxes to [m water equivalent * timestep]
!
  DO 110 jl=1,kproma
     zrcp=1._dp/(cpd*(1._dp+vtmpc2*MAX(0.0_dp,pqm1(jl,klev))))
     zlfdcp(jl)=alf*zrcp
     zdp(jl)=paphm1(jl,klev+1)-paphm1(jl,klev)
     zsnmel(jl)=0._dp
     zsncmelt(jl)=0._dp
     zros(jl)=0._dp
     zdrain(jl)=0._dp
     zsn(jl)=0._dp
     zmlres(jl)=0._dp
     palac(jl)=0._dp
     zrogl(jl)=0._dp
     zraind(jl)=(prsfl(jl)+prsfc(jl))             *zdtime/rhoh2o
     zsnowd(jl)=(pssfl(jl)+pssfc(jl))             *zdtime/rhoh2o
     zevttd(jl)=pevapl(jl)                        *zdtime/rhoh2o
     zevsnd(jl)=pcvs(jl)*pevapot(jl)              *zdtime/rhoh2o
     zevwld(jl)=(1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
     zevwsd(jl)=zevttd(jl)-zevsnd(jl)
 110 END DO
!
!---wiso-code

! Convert water fluxes to [m water equivalent * timestep] - water isotopes
!
 DO jt=1,nwiso
  DO jl=1,kproma
     zwisosnmel(jl,jt)=0._dp
     zwisosncmelt(jl,jt)=0._dp
     zwisoros(jl,jt)=0._dp
     zwisodrain(jl,jt)=0._dp
     zwisosn(jl,jt)=0._dp
     zwisomlres(jl,jt)=0._dp
     pwisoalac(jl,jt)=0._dp
     zwisorogl(jl,jt)=0._dp
     zwisoraind(jl,jt)=(pwisorsfl(jl,jt)+pwisorsfc(jl,jt))*zdtime/twisorhoh2o(jt)
     zwisosnowd(jl,jt)=(pwisossfl(jl,jt)+pwisossfc(jl,jt))*zdtime/twisorhoh2o(jt)
!    calculate isotope ratio of snow reservoir
     zdelta=tnat(jt)
     IF (psn(jl).gt.zwisomin .and. pwisosn(jl,jt).gt.zwisomin) zdelta = pwisosn(jl,jt)/psn(jl)
     IF (abs(1._dp-zdelta).lt.zwisosec) zdelta=1._dp  ! cut off rounding errors

     zwisoevsnd(jl,jt)=pcvs(jl)*zdelta*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
     
!    calculate isotope ratio of skin reservoir
     zdelta=tnat(jt)
     IF (pwl(jl).gt.zwisomin .and. pwisowl(jl,jt).gt.zwisomin) zdelta = pwisowl(jl,jt)/pwl(jl)
     IF (abs(1._dp-zdelta).lt.zwisosec) zdelta=1._dp  ! cut off rounding errors
     zwisoevwld(jl,jt)=(1._dp-pcvs(jl))*pcvw(jl)*zdelta*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)

     zwisoevttd(jl,jt)=pwisoevapl(jl,jt)*zdtime/twisorhoh2o(jt)
     zwisoevwsd(jl,jt)=zwisoevttd(jl,jt)-zwisoevsnd(jl,jt)
  END DO
 END DO

!---wiso-code-end
!
  IF(lsurf) THEN

!
!     ------------------------------------------------------------------
!
!*    2.     Budgets of snow (canopy and ground) and glaciers
!
!*    2.1    Snow changes in the canopy (interception of snowfall,
!            sublimation, melting, unloading due to wind)
!
      DO 210 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            zsn(jl)=zsnowd(jl)+zevsnd(jl)
            zsncmax=MAX(0.0_dp,pwlmx(jl)-cwlmax)
            zsnowd_tmp(jl)=zsnowd(jl)                                          !---wiso-code: remember total snow
            zmprcp=MIN(zsnowd(jl)*cvinter,zsncmax-psnc(jl))
            zmprcp_tmp(jl)=zmprcp                                              !---wiso-code: remember interception
            zsncp=psnc(jl)+zmprcp
            zsnowd(jl)=zsnowd(jl)-zmprcp
            zsncp_tmp(jl)=zsncp                                                !---wiso-code: remember snow canopy
            zevsnd_tmp(jl)=zevsnd(jl)                                          !---wiso-code: remember atmospheric flux
            lo_sub(jl) = (zevsnd(jl).GT.0.0_dp)                                !---wiso-code: flag for sublimation vs. deposition
            psnc(jl)=MIN(MAX(0._dp,zsncp+zevsnd(jl)),zsncmax)
            zsubl_tmp(jl)=psnc(jl)-zsncp                                       !---wiso-code: remember sublimation
            zevsnd(jl)=zevsnd(jl)-(psnc(jl)-zsncp)
            zexpt=MAX(0._dp,ptm1(jl,klev)+3._dp-tmelt)*zdtime/zc2
            zexpw=SQRT(pu10(jl)**2+pv10(jl)**2)*zdtime/zc3
            zsnc_tmp(jl)=psnc(jl)                                              !---wiso-code: remember new snow canopy
            zsncmelt(jl)=psnc(jl)*(1._dp-EXP(-zexpt))
            psnc(jl)=psnc(jl)-zsncmelt(jl)
            zsnc_tmp2(jl)=psnc(jl)                                             !---wiso-code: remember new snow canopy
            zsncwind=psnc(jl)*(1._dp-EXP(-zexpw))
            zsncwind_tmp(jl)=zsncwind                                          !---wiso-code: remember unloading due to wind
            psnc(jl)=psnc(jl)-zsncwind
            zsnowd(jl)=zsnowd(jl)+zsncwind
            ptte(jl,klev)=ptte(jl,klev)-                               &
                            zsncmelt(jl)*zsncfac*zlfdcp(jl)/zdp(jl)
!
!   pwl(jl)=pwl(jl)+zsncmelt(jl) see section 2.5
!
         ELSE
            psnc(jl)=0._dp
         END IF
 210  END DO
!
!---wiso-code

   ! Snow changes in the canopy (interception of snowfall,
   ! sublimation, melting, unloading due to wind) - water isotopes (no additional fractionation)
   DO jt=1,nwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            zwisosn(jl,jt)=zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)

            zdelta=tnat(jt)
            IF (zsnowd_tmp(jl).gt.zwisomin .and. zwisosnowd(jl,jt).gt.zwisomin) zdelta=zwisosnowd(jl,jt)/zsnowd_tmp(jl)
            IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
            zwisomprcp=zmprcp_tmp(jl) * zdelta
	        zwisosncp=pwisosnc(jl,jt)+zwisomprcp
            zwisosnowd(jl,jt)=zwisosnowd(jl,jt)-zwisomprcp

            zdelta=tnat(jt)
            IF (lo_sub(jl)) THEN                           ! deposition, assume identical delta value as the atmospheric flux
              IF (zevsnd_tmp(jl).gt.zwisomin .and. zwisoevsnd(jl,jt).gt.zwisomin) zdelta=zwisoevsnd(jl,jt)/zevsnd_tmp(jl)
            ELSE                                           ! sublimation, assume identical delta value as the canopy layer
              IF (zsncp_tmp(jl).gt.zwisomin .and. zwisosncp.gt.zwisomin) zdelta=zwisosncp/zsncp_tmp(jl)
            ENDIF  
            IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors

            zwisosubl=zsubl_tmp(jl) * zdelta
            pwisosnc(jl,jt)=zwisosncp+zwisosubl
            zwisoevsnd(jl,jt)=zwisoevsnd(jl,jt)-zwisosubl
            
            zdelta=tnat(jt)
            IF (zsnc_tmp(jl).gt.zwisomin .and. pwisosnc(jl,jt).gt.zwisomin) zdelta=pwisosnc(jl,jt)/zsnc_tmp(jl)
            IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
            zwisosncmelt(jl,jt)=zsncmelt(jl) * zdelta
            pwisosnc(jl,jt)=pwisosnc(jl,jt)-zwisosncmelt(jl,jt)

            zdelta=tnat(jt)
            IF (zsnc_tmp2(jl).gt.zwisomin .and. pwisosnc(jl,jt).gt.zwisomin) zdelta=pwisosnc(jl,jt)/zsnc_tmp2(jl)
            IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
            zwisosncwind=zsncwind_tmp(jl) * zdelta
            pwisosnc(jl,jt)=pwisosnc(jl,jt)-zwisosncwind
            zwisosnowd(jl,jt)=zwisosnowd(jl,jt)+zwisosncwind
         ELSE
            pwisosnc(jl,jt)=0._dp
         END IF
      END DO
   END DO

!---wiso-code-end

!
!*    2.2    Snowfall and sublimation on land (excluding glaciers)
!
      DO 220 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            psn(jl)=psn(jl)+zsnowd(jl)+zevsnd(jl)
            zsn_tmp(jl)=psn(jl) ! wiso-code
            IF (psn(jl).LT.0._dp) THEN
               zevwsd(jl)=zevwsd(jl)+psn(jl)
               psn(jl)=0._dp
            END IF
         ELSE
            psn(jl)=0._dp
         END IF
 220  END DO
!
!---wiso-code

   ! Snowfall and sublimation on land (excluding glaciers) - water isotopes
   DO jt=1,nwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            pwisosn(jl,jt)=pwisosn(jl,jt)+zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)
            IF (zsn_tmp(jl).LT.0._dp) THEN
               zwisoevwsd(jl,jt)=zwisoevwsd(jl,jt)+pwisosn(jl,jt)
               pwisosn(jl,jt)=0._dp
            END IF
         ELSE
            pwisosn(jl,jt)=0._dp
         END IF
      END DO
   END DO

!---wiso-code-end
!
!*    2.3    Snowfall and sublimation on glaciers and diagnostics
!
      DO 230 jl=1,kproma
         IF (lpglac(jl)) THEN
            pgld(jl)=pgld(jl)+zsnowd(jl)+zevsnd(jl)
            palac(jl)=zraind(jl)+zsnowd(jl)+zevttd(jl)
            zrogl(jl)=zraind(jl)
         END IF
 230  END DO
!
!---wiso-code

   ! Snowfall and sublimation on glaciers and diagnostics - water isotopes
   DO jt=1,nwiso
      DO jl=1,kproma
         IF (lpglac(jl)) THEN
            pwisogld(jl,jt)=pwisogld(jl,jt)+zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)
            pwisoalac(jl,jt)=zwisoraind(jl,jt)+zwisosnowd(jl,jt)+zwisoevttd(jl,jt)
            zwisorogl(jl,jt)=zwisoraind(jl,jt)
         END IF
      END DO
   END DO

!---wiso-code-end
!
!*    2.4    Snow and glacier melt
!
   IF (.NOT. lstart) THEN
      DO 240 jl=1,kproma
         zgld_tmp(jl)=pgld(jl)                                                 !---wiso-code: remember glacier depth
         ztsl_tmp(jl)=ptsl(jl)                                                 !---wiso-code: remember land surface temperature
         zsn_tmp(jl)=psn(jl)                                                   !---wiso-code: remember snow depth
         IF (lpland(jl).AND.ptsl(jl).GT.tmelt) THEN
            IF (lpglac(jl)) THEN
               zsnmel(jl)=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
               pgld(jl)=pgld(jl)-zsnmel(jl)
               zrogl(jl)=zrogl(jl)+zsnmel(jl)
               ptsl(jl)=tmelt
            ELSE IF (psn(jl).GT.0._dp) THEN
               zsmelt=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
               zsnmel(jl)=MIN(psn(jl),zsmelt)
               ptsl(jl)=ptsl(jl)-zsnmel(jl)*alf*rhoh2o/pgrndcapc(jl)
               psn(jl)=MAX(psn(jl)-zsnmel(jl),0._dp)
            END IF
         END IF
 240  END DO
   END IF
!
!---wiso-code

  ! Snow and glacier melt - water isotopes (no additional fractionation)
   IF (.NOT. lstart) THEN
     DO jt=1,nwiso
       DO jl=1,kproma
         IF (lpland(jl).AND.ztsl_tmp(jl).GT.tmelt) THEN
            IF (lpglac(jl)) THEN
               zdelta=tnat(jt)
               IF (zgld_tmp(jl).gt.zwisomin .and. pwisogld(jl,jt).gt.zwisomin) zdelta=pwisogld(jl,jt)/zgld_tmp(jl)
               IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
	           zwisosnmel(jl,jt)=zsnmel(jl) * zdelta
               pwisogld(jl,jt)=pwisogld(jl,jt)-zwisosnmel(jl,jt)
               zwisorogl(jl,jt)=zwisorogl(jl,jt)+zwisosnmel(jl,jt)
            ELSE IF (zsn_tmp(jl).GT.0._dp) THEN
               zdelta=tnat(jt)
               IF (zsn_tmp(jl).gt.zwisomin .and. pwisosn(jl,jt).gt.zwisomin) zdelta=pwisosn(jl,jt)/zsn_tmp(jl)
               IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
	           zwisosnmel(jl,jt)=zsnmel(jl) * zdelta
               pwisosn(jl,jt)=MAX(pwisosn(jl,jt)-zwisosnmel(jl,jt),0._dp)
	        END IF   
	     END IF   
       END DO
     END DO
   END IF

!---wiso-code-end

!
!*    2.5    Snow budget and meltwater (glacier-free land only)
!
      DO 250 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            pwl(jl)=pwl(jl)+zsncmelt(jl)
            zwl_tmp(jl)=pwl(jl)                                                !---wiso-code: remember skin reservoir
            zsnmlt_tmp(jl)=MAX(0._dp,pwl(jl)-pwlmx(jl))                        !---wiso-code: remember additional meltwater
            zsnmlt=zsnmel(jl)+MAX(0._dp,pwl(jl)-pwlmx(jl))
            pwl(jl)=MIN(pwlmx(jl),pwl(jl))
            zmlres(jl)=zsnmlt
            zsnmel(jl)=zsnmel(jl)+zsncmelt(jl)
            zsn(jl)=zsn(jl)-zsnmel(jl)
         END IF
 250  END DO
!	         
!---wiso-code

   ! Snow budget and meltwater (glacier-free land only) - water isotopes (no additional fractionation)
   DO jt=1,nwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            pwisowl(jl,jt)=pwisowl(jl,jt)+zwisosncmelt(jl,jt)
            zdelta=tnat(jt)
            IF (zwl_tmp(jl).gt.zwisomin .and. pwisowl(jl,jt).gt.zwisomin) zdelta=pwisowl(jl,jt)/zwl_tmp(jl)
            IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors
            zwisosnmlt=zwisosnmel(jl,jt) + (zsnmlt_tmp(jl) * zdelta)
            pwisowl(jl,jt)=pwisowl(jl,jt) - (zsnmlt_tmp(jl) * zdelta)
            zwisomlres(jl,jt)=zwisosnmlt
            zwisosnmel(jl,jt)=zwisosnmel(jl,jt)+zwisosncmelt(jl,jt)
            zwisosn(jl,jt)=zwisosn(jl,jt)-zwisosnmel(jl,jt)
         END IF
      END DO
   END DO

!---wiso-code-end

!
!     ------------------------------------------------------------------
!*    3.     Soil temperatures
!
      CALL soiltemp (   kproma,     kbdim                  &
                       ,ptsl,       ptsoil,    psn         &
                       ,pgrndc,     pgrndd,    pgrndcapc   &
                       ,pgrndhflx,  psodif,    prgcgn      &
                       ,lpland,     lpglac   )
!     ------------------------------------------------------------------
!
!*    4.     Water budget
!
!*    4.1    Skin reservoir (vegetation and bare soil)
!
      DO 410 jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
!
!*    4.1.1  Interception of rain
!
            zraind_tmp(jl)=zraind(jl)                                          !---wiso-code: remember total rain
            zmprcp=MIN(zraind(jl)*cvinter,pwlmx(jl)-pwl(jl))
            zmprcp_tmp(jl)=zmprcp                                              !---wiso-code: remember interception
            zwlp=pwl(jl)+zmprcp
            zraind(jl)=zraind(jl)-zmprcp
!
!*    4.1.2  Evaporation or dew collection
!
            zwlp_tmp(jl)=zwlp                                                  !---wiso-code: remember skin reservoir
            zevwld_tmp(jl)=zevwld(jl)                                          !---wiso-code: remember atmospheric flux
            lo_evap(jl) = (zevwld(jl).GT.0.0_dp)                               !---wiso-code: flag for evaporation vs dew collection
            pwl(jl)=MIN(MAX(0._dp,zwlp+zevwld(jl)),pwlmx(jl))
            zevap_tmp(jl)=pwl(jl)-zwlp                                         !---wiso-code: remember evaporation/dew collection
            zevwsd(jl)=zevwsd(jl)-(pwl(jl)-zwlp)
          ELSE
            pwl(jl)=0._dp
          END IF
 410  END DO
!
!---wiso-code

   ! Water budget - water isotopes
   ! Skin reservoir (vegetation and bare soil) - water isotopes
   DO jt=1,nwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN

            ! Interception of rain - water isotopes

            zdelta=tnat(jt)
            IF (zraind_tmp(jl).gt.zwisomin .and. zwisoraind(jl,jt).gt.zwisomin) zdelta=zwisoraind(jl,jt)/zraind_tmp(jl)
            IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
            zwisomprcp=zmprcp_tmp(jl) * zdelta
            zwisowlp=pwisowl(jl,jt)+zwisomprcp
            zwisoraind(jl,jt)=zwisoraind(jl,jt)-zwisomprcp

            ! Evaporation or dew collection - water isotopes

            zdelta=tnat(jt)
            IF (lo_evap(jl)) THEN                          ! dew collection, assume identical delta value as the atmospheric flux
              IF (zevwld_tmp(jl).gt.zwisomin .and. zwisoevwld(jl,jt).gt.zwisomin) zdelta=zwisoevwld(jl,jt)/zevwld_tmp(jl)
            ELSE                                           ! evaporation, assume identical delta value as the skin reservoir
              IF (zwlp_tmp(jl).gt.zwisomin .and. zwisowlp.gt.zwisomin) zdelta=zwisowlp/zwlp_tmp(jl)
            ENDIF  
            IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp   ! cut off rounding errors
            zwisoevap=zevap_tmp(jl) * zdelta
            pwisowl(jl,jt)=zwisowlp+zwisoevap
	        zwisoevwsd(jl,jt)=zwisoevwsd(jl,jt)-zwisoevap
          ELSE
            pwisowl(jl,jt)=0._dp
          END IF
      END DO
   END DO

!---wiso-code-end
!
!*    4.2    Soil reservoir
!
      DO 420 jl=1,kproma
          IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
             zwptr=0.90_dp*pwsmx(jl)
             zwdtr=0.90_dp*pwsmx(jl)
             zwslim=0.05_dp*pwsmx(jl)
             zconw2=pwsmx(jl)-zwdtr
             zconw3=zdrmax-zdrmin
             zconw4=pwsmx(jl)-zwptr
             zroeff=MAX(0._dp, porostd(jl)-zorvari)   &
                          /(porostd(jl)+zorvars)
             zbws=MAX(MIN(zroeff,0.5_dp),0.01_dp)
             zb1=1._dp+zbws
	         zb1_tmp(jl)=zb1    ! wiso-code
             zbm=1._dp/zb1
             zconw1=pwsmx(jl)*zb1
             zlyeps=0._dp
             zvol=0._dp
             zinfil=0._dp
!
!*    4.2.1  Surface runoff, infiltration and evaporation from soil
!
             IF (zevwsd(jl) >= 0.0_dp) THEN
               zprfl=zmlres(jl)+zraind(jl)+zevwsd(jl)
             ELSE
                pws(jl)=pws(jl)+zevwsd(jl)
               zprfl=zmlres(jl)+zraind(jl)
             END IF
             zprfl_tmp(jl) = zprfl                                             !---wiso-code: remember total precipitation value
             IF (ptsoil(jl,1).LT.tmelt) THEN
                zros(jl)=zprfl
             ELSE
                lo_zros(jl) = (zprfl.GT.0._dp.AND.pws(jl).GT.zwslim)           !---wiso-code: remember if runoff occurs
                IF (zprfl.GT.0._dp.AND.pws(jl).GT.zwslim) THEN
                   IF (pws(jl).GT.pwsmx(jl)) THEN
                      zlyeps=pws(jl)-pwsmx(jl)
                   ELSE
                      zlyeps=0._dp
                   END IF
                   zlysic=(pws(jl)-zlyeps)/pwsmx(jl)
                   zlysic=MIN(zlysic,1._dp)
                   zvol=(1._dp-zlysic)**zbm-zprfl/zconw1
                   zros(jl)=zprfl-(pwsmx(jl)-pws(jl))
                   IF (zvol.GT.0._dp) THEN
                      zros(jl)=zros(jl)+pwsmx(jl)*zvol**zb1
                   END IF
                   zros(jl)=MAX(zros(jl),0._dp)
                   zros_tmp(jl) = zros(jl)                                     !---wiso-code: remember runoff value
                   zinfil=zprfl-zros(jl)
                ELSE
                   zros(jl)=0._dp
                   zinfil=zprfl
                END IF
                pws(jl)=pws(jl)+zinfil
             END IF
!
!*    4.2.2  Drainage and total runoff
!
             zws_tmp(jl)=pws(jl)                                               !---wiso-code: remember value of pws before drainage
             IF (pws(jl).LE.zwslim) THEN
                zdrain(jl)=0._dp
             ELSE
                IF (ptsoil(jl,1).GT.tmelt) THEN
                   zdrain(jl)=zdrmin*pws(jl)/pwsmx(jl)
                   IF (pws(jl).GT.zwdtr) THEN
                      zdrain(jl)=zdrain(jl)+zconw3*                  &
                                ((pws(jl)-zwdtr)/zconw2)**zdrexp
                   END IF
                   zdrain(jl)=zdrain(jl)*zdtime
                   zdrain(jl)=MIN(zdrain(jl),pws(jl)-zwslim)
                   pws(jl)=pws(jl)-zdrain(jl)
                ELSE
                   zdrain(jl)=0._dp
                END IF
             END IF
             zwsup=MAX(pws(jl)-pwsmx(jl),0._dp)
             zwsup_tmp(jl) = zwsup                                             !---wiso-code: remember soil water surplus
             pws(jl)=pws(jl)-zwsup
             zros(jl)=zros(jl)+zdrain(jl)+zwsup
          ELSE
             pws(jl)=0._dp
          END IF
 420   END DO
!
!---wiso-code

   ! Soil reservoir - water isotopes
   DO jt=1,nwiso
      DO jl=1,kproma
          IF (lpland(jl).AND..NOT.lpglac(jl)) THEN

             zwslim=0.05_dp*pwsmx(jl)
             zwisoinfil=0._dp

             ! Surface runoff, infiltration and evaporation from soil - water isotopes
             IF (zevwsd(jl) >= 0.0_dp) THEN
               zwisoprfl=zwisomlres(jl,jt)+zwisoraind(jl,jt)+zwisoevwsd(jl,jt)
             ELSE
               pwisows(jl,jt)=pwisows(jl,jt)+zwisoevwsd(jl,jt)
               zwisoprfl=zwisomlres(jl,jt)+zwisoraind(jl,jt)
             END IF

             ! assume that surface runoff has same delta value as precipitation (no mixing with soil water pool ws)
             zdelta=tnat(jt)
             IF (zprfl_tmp(jl).gt.zwisomin .and. zwisoprfl.gt.zwisomin) zdelta=zwisoprfl/zprfl_tmp(jl)
             IF ((1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors

             IF (ptsoil(jl,1).LT.tmelt) THEN
                zwisoros(jl,jt)=zwisoprfl
             ELSE
                IF (lo_zros(jl)) THEN
                   zwisoros(jl,jt)=zros_tmp(jl)*zdelta
                   zwisoinfil=zwisoprfl-zwisoros(jl,jt)
                ELSE
                   zwisoros(jl,jt)=0._dp
                   zwisoinfil=zwisoprfl
                END IF
                pwisows(jl,jt)=pwisows(jl,jt)+zwisoinfil
             END IF

             ! Drainage and total runoff - water isotopes

             ! drainage water and soil water surplus have same delta value as soil water ws
             zdelta=tnat(jt)
             IF (zws_tmp(jl).gt.zwisomin .and. pwisows(jl,jt).gt.zwisomin) zdelta=pwisows(jl,jt)/zws_tmp(jl)
             IF (ABS(1._dp-zdelta).lt.zwisosec) zdelta=1._dp ! cut off rounding errors

             IF (zws_tmp(jl).LE.zwslim) THEN
                zwisodrain(jl,jt)=0._dp
             ELSE
                IF (ptsoil(jl,1).GT.tmelt) THEN
                   zwisodrain(jl,jt)=zdrain(jl)*zdelta
                   pwisows(jl,jt)=pwisows(jl,jt)-zwisodrain(jl,jt)
                ELSE
                   zwisodrain(jl,jt)=0._dp
                END IF
             END IF
             zwisowsup=zwsup_tmp(jl)*zdelta
             pwisows(jl,jt)=pwisows(jl,jt)-zwisowsup
             zwisoros(jl,jt)=zwisoros(jl,jt)+zwisodrain(jl,jt)+zwisowsup
          ELSE
             pwisows(jl,jt)=0._dp
          END IF
      END DO
   END DO

!---wiso-code-end
!
!*     4.2.3  Runoff and drainage for the HD-Model
!
      DO 423 jl=1,kproma
        pros_hd(jl)=zros(jl)-zdrain(jl)
        pdrain_hd(jl)=zdrain(jl)
 423  END DO
!
!---wiso-code

   ! Runoff and drainage for the HD-Model - water isotopes
   DO jt=1,nwiso
      DO jl=1,kproma
        pwisoros_hd(jl,jt)=zwisoros(jl,jt)-zwisodrain(jl,jt)
        pwisodrain_hd(jl,jt)=zwisodrain(jl,jt)
      END DO
   END DO

!---wiso-code-end
!
!     ------------------------------------------------------------------
!
!*    5.     Time filter for surface temperature
!
   IF (.NOT.lstart) THEN
      DO 510 jl=1,kproma
         ptslm1(jl)=ptslm(jl)+eps*(ptslm1(jl)-2._dp*ptslm(jl)+ptsl(jl))
         ptslm(jl)=ptsl(jl)
 510  END DO
   ELSE
      DO 511 jl=1,kproma
         ptslm1(jl)=ptslm(jl)
 511  END DO
   END IF

!
!     ------------------------------------------------------------------
!
!*    6.     Accumulate fluxes for diagnostics
!
!     6.1    Water fluxes
!
      DO 601 jl=1,kproma
         prunoff(jl)= prunoff(jl) +zros(jl)   *rhoh2o*pslm(jl)
         psnmel(jl) = psnmel(jl)  +zsnmel(jl) *rhoh2o*pslm(jl)
         papmegl(jl)= papmegl(jl) +palac(jl)  *rhoh2o*pslm(jl)
         pdrain(jl) = pdrain(jl)  +zdrain(jl) *rhoh2o*pslm(jl)
         psnacl(jl) = psnacl(jl)  +zsn(jl)    *rhoh2o*pslm(jl)
         progl(jl)  = progl(jl)   +zrogl(jl)  *rhoh2o*pslm(jl)
 601  END DO
!
!---wiso-code

   ! Water fluxes - water isotopes
   DO jt=1,nwiso
      DO jl=1,kproma
         pwisorunoff(jl,jt)= pwisorunoff(jl,jt) +zwisoros(jl,jt)   *twisorhoh2o(jt)*pslm(jl)
         pwisosnmel(jl,jt) = pwisosnmel(jl,jt)  +zwisosnmel(jl,jt) *twisorhoh2o(jt)*pslm(jl)
         pwisoapmegl(jl,jt)= pwisoapmegl(jl,jt) +pwisoalac(jl,jt)  *twisorhoh2o(jt)*pslm(jl)
         pwisodrain(jl,jt) = pwisodrain(jl,jt)  +zwisodrain(jl,jt) *twisorhoh2o(jt)*pslm(jl)
         pwisosnacl(jl,jt) = pwisosnacl(jl,jt)  +zwisosn(jl,jt)    *twisorhoh2o(jt)*pslm(jl)
         pwisorogl(jl,jt)  = pwisorogl(jl,jt)   +zwisorogl(jl,jt)  *twisorhoh2o(jt)*pslm(jl)
      END DO
   END DO

!Corrections of minor water pools
   DO jt = 1,nwiso
     DO jl = 1,kproma
       IF (pwl(jl).LT.zwisomin) pwisowl(jl,jt)=pwl(jl)*tnat(jt)
       IF (pws(jl).LT.zwisomin) pwisows(jl,jt)=pws(jl)*tnat(jt)
       IF (psn(jl).LT.zwisomin) pwisosn(jl,jt)=psn(jl)*tnat(jt)
       IF (psnc(jl).LT.zwisomin) pwisosnc(jl,jt)=psnc(jl)*tnat(jt)
       IF (pgld(jl).LT.zwisomin) pwisogld(jl,jt)=pgld(jl)*tnat(jt)
     END DO
   END DO

!---wiso-code-end

!
!     6.2     Ground heat flux
!
      DO 602 jl=1,kproma
         pgrndflux(jl)=pgrndflux(jl)+pgrndhflx(jl)*pslm(jl)*zdtime
 602  END DO

!
!     ------------------------------------------------------------------
!
!*    7.     Set variables for lsurf = .FALSE. (*surf* by-passed)
!
  ELSE
!
      DO 700 jl=1,kproma
         ptslm1(jl)=ptslm(jl)
         ptsl(jl)=ptslm(jl)
 700  END DO
!
  END IF
!
  RETURN
END SUBROUTINE surf
