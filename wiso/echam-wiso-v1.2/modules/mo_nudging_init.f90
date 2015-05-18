MODULE mo_nudging_init
!BOP
  ! !MODULE: mo_nudging_init (layer 1)

  ! !DESCRIPTION: 
  ! initialisation, function control and clean up of nudging procedure 

  ! !REVISION HISTORY: 
  ! I. Kirchner, MPI Hamburg, April-2001
  ! R. Johanni, IPP Garching, May-2002, parallel version
  ! I. Kirchner, MPI Hamburg, Aug-2002, revision

!BOX
  USE mo_nudging_constants, ONLY: &
       lnudgdbx, lnudgini, lnudgimp, lnudgfrd, lnudgwobs, ldamplin, &
       ltintlin, nudgdamp, nudgdsize, lnudgpat, lnudgcli, &
       nudgd, nudgv, nudgt, nudgp, &
       nudgtrun, nudgsmin, nudgsmax, nudglmin, nudglmax, &
       lnudgstop, lsite, nudgmin, nudgmax

  USE mo_nudging_sst,  ONLY: nsstinc, nsstoff, ndg_freez, ndg_file_sst
  USE mo_nudging_io,   ONLY: ndg_file_vor, ndg_file_div, ndg_file_stp
  USE mo_time_control, ONLY: dt_nudg_start, dt_nudg_stop

  IMPLICIT NONE

  PRIVATE
!EOX

  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC  :: NudgingInit

  ! !PUBLIC DATA MEMBERS:
  INTEGER, PARAMETER, PUBLIC ::   &!--- function switch for initialization
       NDG_INI_IO     = 0, &! initialize namelist parameters
       NDG_INI_STREAM = 1, &! initialize stream
       NDG_INI_MEM    = 2, &! allocate memory
       NDG_CLEAN_MEM  = 3, &! deallocate memory
       NDG_CLOSE      = 4   ! close nudging
!BOX
  CHARACTER(len=256) :: mess
!EOX
!EOP
!BOP
    ! !NAMELIST: NDGCTL
    ! !DESCRIPTION:
    ! external user interface of nudging

    ! !INPUT PARAMETERS: 
    NAMELIST /ndgctl/&
         lnudgini,          &! data synchronisation
         lnudgcli,          &! use nudging data cyclic
         dt_nudg_start,     &! nudging start date
         dt_nudg_stop,      &! nudging stop date
         ndg_file_stp,      &! nudging data, temperature and log surface pressure
         ndg_file_div,      &! nudging data, divergence
         ndg_file_vor,      &! nudging data, vorticity
         nsstinc,           &! sst handling, read new sst after NSSTINC hours
         nsstoff,           &! read sst at hour NSSTOFF first
         ndg_freez,         &! sea ice mask detection value
         ndg_file_sst,      &! sst data file
         lnudgimp,          &! nudging method, implicit/explicit
         lnudgpat,          &! pattern nudging
         lnudgfrd,          &! define place for NMI filter
         nudglmin,          &! vertical level separation, index of uppermost layer
         nudglmax,          &! index of deepest layer
         nudgtrun,          &! spectral domain selection
         nudgsmin,          &! lowest nudged wavenumber
         nudgsmax,          &! highest nudged wavenumber
         nudgp,             &! relaxation time definition for surface parameter
         nudgt,             &! relaxation time definition for temperature
         nudgd,             &! relaxation time definition for divergence
         nudgv,             &! relaxation time definition for vorticity
         ltintlin,          &! time interpolation linear/cubic spline
         ldamplin,          &! change relaxation between synoptic times
         nudgdamp,          &! 1= no damping, 0= full damping
         nudgdsize,         &! environment radius for damping
         lnudgdbx,          &! diagnostics
         lnudgwobs,         &! store reference fields
         lsite               ! tendency diagnostics
!EOP
CONTAINS

  !======================================================================
!BOP
  !
  ! !IROUTINE:  NudgingInit
  ! !INTERFACE:

  SUBROUTINE NudgingInit(ifunction)

    ! !DESCRIPTION: 
    ! setup nudging parameter and fields

    ! !USES:

    !* general modules
    USE mo_kind,          ONLY: dp
    USE mo_exception,     ONLY: message, finish
    USE mo_control,       ONLY: &
         lnudge, ltdiag, n2sp, nlev, nm, nlevp1, lnmi, nsp
    USE mo_filename,      ONLY: find_next_free_unit, out_expname, str_filter
    USE mo_physc2,        ONLY: ctfreez
    USE mo_namelist,      ONLY: position_nml, nnml, POSITIONED
    USE mo_time_conversion, ONLY: OPERATOR(<)
    USE mo_time_control, ONLY: &
         dt_nudg_start, dt_nudg_stop, lstart, lresume, &
         init_nudgingtime_a, &
         init_nudgingtime_b, &
         init_nudgingtime_c, &
         lfirst_cycle, get_date_components, inp_convert_date, &
         nudg_start, next_date, nudg_stop, ndg_inp_date, &
         delta_time, ndg_date0, ndg_date1

    !* MPI functions
    USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_io, p_bcast
    USE mo_truncation,    ONLY: nmp, nnp
    USE mo_decomposition, ONLY: &
         ldc => local_decomposition, gdc => global_decomposition
    USE mo_transpose,     ONLY: scatter_sp

    !* nudging modules
    USE mo_nudging_constants, ONLY: &
         ndg_version, &
         lnudgdbx, lnudgini, lnudgimp, lnudgfrd, lnudgwobs, ldamplin, &
         ltintlin, nudgdamp, nudgdsize, lnudgpat, lnudgcli, &
         nudgd, nudgv, nudgt, nudgp, &
         nudgtrun, nudgsmin, nudgsmax, nudglmin, nudglmax, &
         lnudgstop, nmaxc, nfact, lsite, ndunit, &
         linp_vor, linp_div, linp_tem, linp_lnp, &
         lnudg_run, nudgmin, nudgmax, &
         ilev_v_max, ilev_d_max, ilev_t_max, &
         ilev_v_min, ilev_d_min, ilev_t_min, &
         ino_v_lev,  ino_d_lev,  ino_t_lev, &
         NDG_WINDOW_ALL, NDG_WINDOW_CUT, NDG_WINDOW_CUT0

    USE mo_nudging_sst, ONLY: &
         nsstinc, nsstoff, ndg_freez, ndg_file_sst, NudgingSSTClose
  
    USE mo_nudging_io, ONLY: &
         ndg_file_vor, ndg_file_div, ndg_file_stp, ndgblock, &
         lndgstep0, lndgstep1, lndgstep2, lndgstep3, CloseBlock

    USE mo_nudging_buffer, ONLY: &
         flagn, nudgdo, nudgvo, nudgto, nudgpo, &
         NudgingAllocMem, NudgingDeallocMem, &
         NudgingAllocRefMem, NudgingDeallocRefMem, &
         NudgingStreamInit, nudg, NdgInitCounter

    USE mo_nudging_utils, ONLY: cpbread, HEAD_LEN, WORD_LEN

    USE mo_memory_base,   ONLY: set_stream

!EOP
!BOC
!BOX
    INTEGER, INTENT(in)  :: ifunction    ! function selector

    INTEGER  :: &
         i, kret, ibytes, snsp, &
         ihead(8), is1, is2, jm, ii, sum1, sum2, &
         yr, mo, dy, hr, mi, se, &
         ipos ! error return value from position_nml
    REAL(dp)      :: twodt
    LOGICAL            :: found
    CHARACTER(len=8)   :: yhead(HEAD_LEN)
    CHARACTER(len=256) :: cfile
    INTEGER            :: kfile
    REAL(dp), POINTER :: flagn_global(:,:,:), flagn_local(:,:,:)

    INTRINSIC MAX, MIN

    snsp = ldc%snsp

    IF (.NOT. lnudge) RETURN

    SELECT CASE(ifunction)

! ==============================================================================
!
!EOX
! ***** initialise constant parameters **** PART 1 nudging init

    CASE(NDG_INI_IO) !==========================================================
!BOX

      !*** read namelist *******************************************************

      IF (p_parallel_io) THEN
        CALL position_nml ('NDGCTL', status=ipos)
        SELECT CASE (ipos)
        CASE (POSITIONED);  READ (nnml,ndgctl)
        END SELECT
      ENDIF

      IF (p_parallel) THEN
        CALL bc_nudging
        IF (lnudgpat) CALL finish('NudgingInit',&
             'Pattern nudging not allowed in parallel mode.')
      END IF

      CALL message('',ndg_version)



!EOX
      !*** define data channel names *******************************************
!BOX
      WRITE(mess,*) ' File name template VOR : ',TRIM(ndg_file_vor)
      CALL message('',mess)
      WRITE(mess,*) ' File name template DIV : ',TRIM(ndg_file_div)
      CALL message('',mess)
      WRITE(mess,*) ' File name template STP : ',TRIM(ndg_file_stp)
      CALL message('',mess)

      IF (ltdiag .AND. p_parallel_io) THEN
         ndunit = find_next_free_unit(90,99)     ! search next free unit in range 90-99
         cfile = TRIM(out_expname) // '_ndgcorr'
         WRITE(mess,*) &
              ' Correlation results on unit  ', ndunit,' file: ',TRIM(cfile)
         CALL message('',mess)
         INQUIRE(file=cfile,exist=found)
         IF (found) THEN
            OPEN(unit=ndunit,file=cfile,position='append',form='formatted')
         ELSE
            OPEN(unit=ndunit,file=cfile,status='replace',form='formatted')
         END IF
         WRITE(ndunit,'(a)') '### start/append file with tendency correlations'
      END IF
      CALL message('','')



!EOX
      !*** define SST dependent parameters *************************************

      !--- nsstinc equals 0 means no external sst is used
!BOX
      nsstinc = MAX(MIN(nsstinc,24),0)
      nsstoff = MAX(MIN(nsstoff,23),0)
      IF (ndg_freez < 0) THEN
         CALL message('Warning',&
              'Predefined NDG_FREEZ not correct ... reset to CTFREEZ')
         ndg_freez = ctfreez
      END IF

      !--- initialize sst memory

      IF (nsstinc==0) THEN
         CALL message('',' Use standard ECHAM SST')
      ELSE IF (lfirst_cycle) THEN
         WRITE(mess,*) ' File name template SST : ',TRIM(ndg_file_sst)
         CALL message('',mess)
         WRITE(mess,*) ' Use sea ice detection limit NDG_FREEZ = ',ndg_freez
         CALL message('',mess)
      ENDIF



!EOX
      !*** general nudging options *********************************************
!BOX
      CALL message('','Nudging options:')

      IF (lnudgpat) THEN
         lnudgimp = .FALSE.
         CALL message('','Force explicit method.')
      END IF
      IF (lstart) THEN
         lnudgini = .TRUE.
         CALL message('','Force time synchronization with nudging data.')
      END IF
      WRITE(mess,*) '  pattern assimilation            LNUDGPAT  = ',lnudgpat
      CALL message('',mess)
      WRITE(mess,*) '  Adjust date/time                LNUDGINI  = ',lnudgini
      CALL message('',mess)
      WRITE(mess,*) '  Implicit method                 LNUDGIMP  = ',lnudgimp
      CALL message('',mess)
      WRITE(mess,*) '  Debbugging                      LNUDGDBX  = ',lnudgdbx
      CALL message('',mess)
      WRITE(mess,*) '  Store reference data in output  LNUDGWOBS = ',lnudgwobs
      CALL message('',mess)
      WRITE(mess,*) '  Calculate SITE measure          LSITE     = ',lsite
      CALL message('',mess)
      WRITE(mess,*) '  Use nudging data cyclic         LNUDGCLI  = ',lnudgcli
      CALL message('',mess)
      WRITE(mess,*) '  Linear time interpolation       LTINTLIN  = ',ltintlin
      CALL message('',mess)

      CALL message('','   Time window options ....')
      WRITE(mess,*) '    Linear damping in time        LDAMPLIN  = ',ldamplin
      CALL message('',mess)
      nudgdamp = MIN(MAX(nudgdamp,0.0_dp),1.0_dp)
      WRITE(mess,*) '    Reduce damping between reference times ',&
           nudgdamp,' (fraction)'
      CALL message('',mess)
      nudgdsize = MIN(MAX(nudgdsize,0.0_dp),0.5_dp)
      WRITE(mess,*) '    Damping radius                         ',&
           nudgdsize,' (fraction)'
      CALL message('',mess)



!EOX
      !*** define nudging start and stop date **********************************
!BOX
      lnudgstop = init_nudgingtime_a()

!EOX
      !*** synchronize model date with nudging files ***************************
!BOX
      IF (lfirst_cycle) THEN

         DO
            !--- compose the file name used for time adjustment

            IF (lnudgini) THEN
               !--- use nudging start date for file name evaluation
               CALL message('','Nudging initial time adjustment ...')
               CALL get_date_components(nudg_start,yr,mo,dy,hr,mi,se)

            ELSE IF (lresume) THEN
               !--- use next model date for name evaluation
               CALL message('','Nudging rerun time adjustment ...')
               CALL get_date_components(next_date,yr,mo,dy,hr,mi,se)

            ELSE

               EXIT                                    ! no evaluation needed
            END IF

            !--- which data set must be available ?
            DO i = 1,nmaxc
               IF (nudgd(i) >= 0._dp) THEN; linp_div = .TRUE.; EXIT; END IF
               IF (nudgv(i) >= 0._dp) THEN; linp_vor = .TRUE.; EXIT; END IF
               IF (nudgt(i) >= 0._dp) THEN; linp_tem = .TRUE.; EXIT; END IF
            END DO
                     
            IF (linp_div) THEN
               cfile = str_filter(ndg_file_div,yr,mo,dy,hr,mi,se,ndgblock)
            ELSE IF (linp_vor) THEN
               cfile = str_filter(ndg_file_vor,yr,mo,dy,hr,mi,se,ndgblock)
            ELSE IF (linp_tem) THEN
               cfile = str_filter(ndg_file_stp,yr,mo,dy,hr,mi,se,ndgblock)
            ELSE
               CALL finish('NudgingInit','No external data defined.')
            END IF

            WRITE(mess,*) '   Adjust date using file: ',TRIM(cfile)
            CALL message('',mess)

            IF (p_parallel_io) THEN
              INQUIRE(file=cfile,exist=found)
              IF (.NOT.found) &
                 CALL finish('NudgingInit','Nudging data file not found.')

              CALL pbopen  (kfile, cfile, 'r', kret)
              CALL cpbread (kfile, yhead, HEAD_LEN*WORD_LEN, kret)
              CALL util_i8toi4 (yhead(1), ihead(1),HEAD_LEN)
              ! check dimension
              IF (ihead(5)*ihead(6) /= n2sp*nlev) CALL finish('NudgingInit',&
                   'Dimension of nudging data mismatch with model.')

            ENDIF

            IF (p_parallel) CALL p_bcast(ihead, p_io)
            CALL inp_convert_date(ihead(3),ihead(4)*10000, ndg_date0)

            IF (p_parallel_io) THEN
              ! skip first record and read second header
              ibytes = (ihead(5)*ihead(6) + HEAD_LEN) * WORD_LEN
              CALL pbseek  (kfile, ibytes, 0, kret)
              CALL cpbread (kfile, yhead, HEAD_LEN*WORD_LEN, kret)
              CALL util_i8toi4 (yhead(1), ihead(1), HEAD_LEN)
              CALL pbclose(kfile,kret)
            ENDIF

            IF (p_parallel) CALL p_bcast(ihead, p_io)
            CALL inp_convert_date(ihead(3),ihead(4)*10000, ndg_date1)

            IF (lnudgini) THEN
               ! ndg_inp_date defines the nudging data file name
               ! which will be used later for file name generation
               ndg_inp_date = ndg_date0

               IF (ltintlin) THEN
                  ! linear interpolation needs only two synoptic times
                  ! therefore nudging will start at the first time
                  IF (nudg_start < ndg_date0) nudg_start = ndg_date0
               ELSE
                  ! nonlinear time interpolation needs four time levels
                  ! therefore nudging will start at the second time
                  IF (nudg_start < ndg_date1) nudg_start = ndg_date1
               END IF

               ! calculate nudging start date and adjust echam time manager
               CALL init_nudgingtime_b

             ELSE
               ! for the rerun case the first nudging date needed must be
               ! recalculated otherwise the automatic evaluation of file 
               ! names will fail
               CALL init_nudgingtime_c(ltintlin)

            END IF

            EXIT
         END DO
      END IF


      CALL message('','Nudging basic initialization done.')
      CALL message('','')
!EOX
   CASE(NDG_INI_STREAM) ! =========================================================

      !*** define the output stream ********************************************
!BOX
      CALL NudgingStreamInit

      CALL message('','Nudging output stream initialization done.')
      CALL message('','')


! ==============================================================================
!
!EOX
! *****  Setup nudging coefficients and work space **** PART2 nudging init

   CASE(NDG_INI_MEM) ! =========================================================
!BOX
      IF (lnudg_run) THEN
         IF (lnudgstop) THEN
            IF (nudg_stop < next_date) CALL StopNudging
         END IF
         RETURN               ! initialisation was performed
      END IF

      IF (next_date < nudg_start) RETURN  ! no nudging for next time step

      CALL message('','***********************************************')
      CALL message('','')
      CALL message('','Start nudging mode.')
      CALL message('','')
      CALL message('','***********************************************')

      lnudg_run = .TRUE.

!EOX
      !*** setup spectral window for nudging ***********************************
!BOX
      WRITE(mess,*) 'Nudging window options:'
      CALL message('',mess)

      nudgsmin = MIN(MAX(nudgsmin,-1),nm)
      nudgsmax = MIN(MAX(nudgsmax, 0),nm)
      nudgsmax = MAX(nudgsmax,nudgsmin)

      WRITE(mess,*) '  Modified wavenumber(s) ',nudgsmin,' ... ',nudgsmax
      CALL message('',mess)

      IF (nudgsmin<0) THEN
         CALL message('','     global average included')
         nudgmin  = 1
         nudgsmin = 0

      ELSE IF (nudgsmin == 0) THEN
         CALL message('','     global average NOT included')
         nudgmin = 2

      ELSE
         nudgmin = nmp(nudgsmin+1)+1
      ENDIF

      nudgmax = nmp(nudgsmax+1)+nnp(nudgsmax+1)
      nudgmax = MAX(nudgmax,nudgmin)
      WRITE(mess,*) '  Modified spectral index SPECC = ',nudgmin, ' ... ',nudgmax
      CALL message('',mess)


!EOX
      !*** setup level window for nudging **************************************
!BOX
      nudglmin = MIN(MAX(nudglmin,1),nlev)
      nudglmax = MIN(MAX(nudglmax,1),nlev)
      nudglmax = MAX(nudglmax,nudglmin)


!EOX
      !*** compose nudging flag field for level-spectral domain ****************
!BOX
      !     For parallel version: Allocate a global flag field, 
      !     set the global field and scatter it to the local field

      nudgtrun = MAX(MIN(nudgtrun,2),0)
      IF (.NOT.ALLOCATED(flagn)) ALLOCATE(flagn(nlevp1,2,snsp))
      flagn_local => flagn; flagn = 0

      ALLOCATE(flagn_global(nlevp1,2,nsp)); flagn_global(:,:,:) = 0

      SELECT CASE(nudgtrun)
      CASE (NDG_WINDOW_ALL)
         WRITE(mess,*) '  Truncation (',nudgtrun,') ... USE all meridional parts'
         flagn_global(nudglmin:nudglmax,:,nudgmin:nudgmax) = 1
         flagn_global(nlevp1           ,:,nudgmin:nudgmax) = 1

      CASE (NDG_WINDOW_CUT)
         WRITE(mess,*) '  Truncation (',nudgtrun,')... cut off meridional parts'
         DO jm=nudgsmin,nudgsmax
            is1 = nmp(jm+1) + 1
            is2 = nmp(jm+1) + 1 + nudgsmax - jm
            flagn_global(nudglmin:nudglmax,:,is1:is2) = 1
            flagn_global(nlevp1           ,:,is1:is2) = 1
         END DO

      CASE (NDG_WINDOW_CUT0)
         WRITE(mess,*) '  Truncation (',nudgtrun,')... cut off except wave 0'
         DO jm=nudgsmin,nudgsmax
            IF (jm==0) THEN
               is1 = nmp(1) + 1
               is2 = nmp(1) + nnp(1)
            ELSE
               is1 = nmp(jm+1) + 1
               is2 = nmp(jm+1) + 1 + nudgsmax - jm
            END IF
            flagn_global(nudglmin:nudglmax,:,is1:is2) = 1
            flagn_global(nlevp1           ,:,is1:is2) = 1
         END DO
      END SELECT

      CALL message('',mess)
      WRITE(mess,*) &
           '      [',NDG_WINDOW_ALL,' ... all meridional parts used]'
      CALL message('',mess)
      WRITE(mess,*) &
           '      [',NDG_WINDOW_CUT,' ... triangular cut-off]'
      CALL message('',mess)
      WRITE(mess,*) &
           '      [',NDG_WINDOW_CUT0,' ... triangular cut-off except wave 0]'
      CALL message('',mess)

      !--- check flag field
      sum1 = 0; sum2 = 0
      DO i=1,nlevp1
         DO ii=1,nsp
            sum1 = sum1 + flagn_global(i,1,ii)+flagn_global(i,2,ii)
            sum2 = sum2 + (1-flagn_global(i,1,ii)) + (1-flagn_global(i,2,ii))
         END DO
      END DO
      IF ((sum1+sum2)/=(n2sp*nlevp1)) THEN
         WRITE(mess,*) ' FLAGN= ',sum1,' FLAGO= ',sum2,' NO= ',n2sp*nlevp1
         CALL message('',mess)
         CALL finish('NudgingInit','ERROR nudging flag field mismatch')
      END IF
      CALL message('','')

      CALL scatter_sp(flagn_global,flagn_local,gdc)

      DEALLOCATE(flagn_global) ! not needed any more


!EOX
      !*** which data channels are needed ? ************************************
!BOX  
      !--- find maximum level again, but using now correct number of levels
      linp_tem = .FALSE.
      linp_div = .FALSE.
      linp_vor = .FALSE.
      DO i = 1,nlev
         IF (nudgt(i) >= 0._dp) THEN; linp_tem = .TRUE.; ilev_t_max = i; END IF
         IF (nudgd(i) >= 0._dp) THEN; linp_div = .TRUE.; ilev_d_max = i; END IF
         IF (nudgv(i) >= 0._dp) THEN; linp_vor = .TRUE.; ilev_v_max = i; END IF
      END DO
           
      !--- clear upper part in coefficient arrays
      IF (linp_tem .AND. ilev_t_max < nlev) nudgt(ilev_t_max+1:nlev) = 0.0_dp
      IF (linp_div .AND. ilev_d_max < nlev) nudgd(ilev_d_max+1:nlev) = 0.0_dp
      IF (linp_vor .AND. ilev_v_max < nlev) nudgv(ilev_v_max+1:nlev) = 0.0_dp
 
      !--- initialize lower bound
      ilev_t_min = ilev_t_max
      ilev_d_min = ilev_d_max
      ilev_v_min = ilev_v_max
 
      !-- find minimum level
      DO i = ilev_t_max,1,-1; IF (nudgt(i) < 0._dp) EXIT; ilev_t_min = i
      END DO
      DO i = ilev_d_max,1,-1; IF (nudgd(i) < 0._dp) EXIT; ilev_d_min = i
      END DO
      DO i = ilev_v_max,1,-1; IF (nudgv(i) < 0._dp) EXIT; ilev_v_min = i
      END DO


!EOX
      !*** rescale nudging coefficients and merge with level window ************
!BOX
      !--- allocate weights of observations
      IF (nlev > nmaxc) CALL finish('NudgingInit','number of levels too large')
      IF (.NOT. ALLOCATED(nudgdo)) &
           ALLOCATE(nudgdo(nlev), nudgvo(nlev), nudgto(nlev))
      nudgto(:) = 0.0_dp ; nudgdo(:) = 0.0_dp ; nudgvo(:) = 0.0_dp

      IF (nudglmin > 1) THEN
         nudgt(1:nudglmin-1) = 0.0_dp
         nudgd(1:nudglmin-1) = 0.0_dp
         nudgv(1:nudglmin-1) = 0.0_dp
      ENDIF
      IF (nudglmax < nlev) THEN
         nudgt(nudglmax+1:nlev) = 0.0_dp
         nudgd(nudglmax+1:nlev) = 0.0_dp
         nudgv(nudglmax+1:nlev) = 0.0_dp
      ENDIF

      !--- clear lower part in coefficient arrays
      IF (ilev_t_min > 1) nudgt(1:ilev_t_min-1) = 0.0_dp
      IF (ilev_d_min > 1) nudgd(1:ilev_d_min-1) = 0.0_dp
      IF (ilev_v_min > 1) nudgv(1:ilev_v_min-1) = 0.0_dp
 
      IF (linp_tem) ino_t_lev = ilev_t_max - ilev_t_min + 1
      IF (linp_div) ino_d_lev = ilev_d_max - ilev_d_min + 1
      IF (linp_vor) ino_v_lev = ilev_v_max - ilev_v_min + 1
      IF (nudgp >= 0._dp) THEN
         linp_lnp = .TRUE.
         ino_t_lev = ino_t_lev + 1
         IF (.NOT. linp_tem) ilev_t_min = 1
      ELSE
         nudgp = 0.0_dp
      END IF

      CALL message('','Define external data amount:')
      IF (linp_tem) THEN
         WRITE(mess,*) &
              '   Use temperature LEVELS = (',ilev_t_min,':',ilev_t_max,')'
         CALL message('',mess)
      END IF
      IF (linp_lnp) THEN
         WRITE(mess,*) '   Use surface fields (log surface pressue)'
         CALL message('',mess)
      END IF
      IF (linp_div) THEN
         WRITE(mess,*) &
              '   Use divergence  LEVELS = (',ilev_d_min,':',ilev_d_max,')'
         CALL message('',mess)
      END IF
      IF (linp_vor) THEN
         WRITE(mess,*) &
              '   Use vorticity   LEVELS = (',ilev_v_min,':',ilev_v_max,')'
         CALL message('',mess)
      END IF

      !--- correct nudging coefficients
      IF (.NOT. linp_tem) nudgt(:) = 0.0_dp
      IF (.NOT. linp_lnp) nudgp    = 0.0_dp
      IF (.NOT. linp_div) nudgd(:) = 0.0_dp
      IF (.NOT. linp_vor) nudgv(:) = 0.0_dp

      !--- merge level intervalls in nudging coefficient field
      IF (linp_vor .AND. linp_div .AND. linp_tem) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_v_min, ilev_d_min, ilev_t_min))
         nudglmax = MIN(nudglmax, MAX(ilev_v_max, ilev_d_max, ilev_t_max))

      ELSE IF (linp_vor .AND. linp_div ) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_v_min, ilev_d_min))
         nudglmax = MIN(nudglmax, MAX(ilev_v_max, ilev_d_max))

      ELSE IF (linp_vor .AND. linp_tem) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_v_min, ilev_t_min))
         nudglmax = MIN(nudglmax, MAX(ilev_v_max, ilev_t_max))

      ELSE IF (linp_div .AND. linp_tem) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_d_min, ilev_t_min))
         nudglmax = MIN(nudglmax, MAX(ilev_d_max, ilev_t_max))

      ELSE IF (linp_vor) THEN
         nudglmin = MAX(nudglmin, ilev_v_min)
         nudglmax = MIN(nudglmax, ilev_v_max)

      ELSE IF (linp_div) THEN
         nudglmin = MAX(nudglmin, ilev_d_min)
         nudglmax = MIN(nudglmax, ilev_d_max)

      ELSE IF (linp_tem) THEN
         nudglmin = MAX(nudglmin, ilev_t_min)
         nudglmax = MIN(nudglmax, ilev_t_max)

      END IF

      WRITE(mess,*) '  Assimilate external data into ',&
           ' LEVEL(s) = ',nudglmin,' ... ',nudglmax
      CALL message('',mess)



!EOX
      !*** prepare nudging weights level dependend *****************************

!BOX
      ! set weights for explicit nudging: new = N(read)*2*dt*nfact

      !--- rescale coefficients
      twodt = 2._dp*delta_time
      nudgpo    = nudgp         * nfact * twodt
      nudgto(:) = nudgt(1:nlev) * nfact * twodt
      nudgdo(:) = nudgd(1:nlev) * nfact * twodt
      nudgvo(:) = nudgv(1:nlev) * nfact * twodt

      !--- set weights for model data
      IF (lnudgpat) THEN
         nudgp          = 1._dp
         nudgt(1:nlev)  = 1._dp
         nudgd(1:nlev)  = 1._dp
         nudgv(1:nlev)  = 1._dp
      ELSE
         nudgp          = 1._dp-nudgpo
         nudgt(1:nlev)  = 1._dp-nudgto(:)
         nudgd(1:nlev)  = 1._dp-nudgdo(:)
         nudgv(1:nlev)  = 1._dp-nudgvo(:)
      END IF
             
      !--- convert coefficients for implicit nudging
      IF (lnudgimp) THEN
         !--- using implicit nudging reset weights
         nudgp         = 1._dp/(1._dp+nudgpo)
         nudgt(1:nlev) = 1._dp/(1._dp+nudgto(:))
         nudgd(1:nlev) = 1._dp/(1._dp+nudgdo(:))
         nudgv(1:nlev) = 1._dp/(1._dp+nudgvo(:))
             
         !--- reset weights for observations
         nudgpo    = nudgpo    * nudgp
         nudgto(:) = nudgto(:) * nudgt(1:nlev) 
         nudgdo(:) = nudgdo(:) * nudgd(1:nlev) 
         nudgvo(:) = nudgvo(:) * nudgv(1:nlev) 
      ENDIF

      !--- print out the time constant part of the nudging weights
      CALL message('',' ')
      WRITE(mess,*) &
           'Nudging coefficients (1/sec) and weights for model and reference (%)'
      CALL message('',mess)
      WRITE(mess,'(a5,3a25,7x)') 'lev','TEM','DIV','VOR'
      CALL message('',mess)

      IF (lnudgimp) THEN
         CALL message('','    IMPLICIT nudging used')
         DO i=1,nlev
            WRITE(mess,'(i5,3(e12.3,2f10.3))') i,&
                 (1.0_dp-nudgt(i))/(nudgt(i)*twodt),nudgt(i)*100.0_dp,nudgto(i)*100.0_dp,&
                 (1.0_dp-nudgd(i))/(nudgd(i)*twodt),nudgd(i)*100.0_dp,nudgdo(i)*100.0_dp,&
                 (1.0_dp-nudgv(i))/(nudgv(i)*twodt),nudgv(i)*100.0_dp,nudgvo(i)*100.0_dp
            CALL message('',mess)
         ENDDO
         WRITE(mess,'(a,e12.3,2f10.3)') ' surface fields = ',&
              (1.0_dp-nudgp)/(nudgp*twodt),nudgp*100.0_dp,nudgpo*100._dp
         CALL message('',mess)

      ELSE
         CALL message('','    EXPLICIT nudging used')
         DO i=1,nlev
            WRITE(mess,'(i5,3(e12.3,2f10.3))') i,&
                 nudgto(i)/twodt,nudgt(i)*100.0_dp,nudgto(i)*100.0_dp,&
                 nudgdo(i)/twodt,nudgd(i)*100.0_dp,nudgdo(i)*100.0_dp,&
                 nudgvo(i)/twodt,nudgv(i)*100.0_dp,nudgvo(i)*100.0_dp
            CALL message('',mess)
         ENDDO
         WRITE(mess,'(a,e12.3,2f10.3)') ' surface fields = ',&
              nudgpo/twodt,nudgp*100.0_dp,nudgpo*100.0_dp
         CALL message('',mess)

      ENDIF

      CALL message('','')


!EOX
      !*** allocate memory *****************************************************
!BOX
      IF (.NOT. lndgstep3) CALL NudgingAllocRefMem

      CALL NudgingAllocMem

      CALL set_stream(nudg, lpout=.TRUE.)  ! switch output stream on

      CALL NdgInitCounter

      lnudgini = .FALSE.



! ==============================================================================
!
!EOX
!****** clean up work space memory

   CASE(NDG_CLEAN_MEM)
!BOX
      IF (lnudg_run) THEN
         CALL message('','Clean up nudging memory:')
         CALL NudgingDeallocMem
         lnudg_run = .FALSE.
      ELSE
         CALL message('','No nudging memory allocated.')
      END IF


! ==============================================================================
!
!EOX
!****** clean up i/o memory and close nudging

   CASE (NDG_CLOSE)
!BOX
      IF (lnmi) lnmi = .FALSE.
      IF (ltdiag .AND. p_parallel_io) CLOSE(ndunit)

      CALL  NudgingDeallocRefMem

      lndgstep0 = .FALSE.
      lndgstep1 = .FALSE.
      lndgstep2 = .FALSE.
      lndgstep3 = .FALSE.

      CALL CloseBlock

! Changes to NudgingiInit routine to allow "sea ice bug" correction via the nsstinc parameter
! For details see E-Mail by Sebastian Rast, MPI Met, HH (05.01.2011)
! M. Werner, AWI, 01/2011
!      CALL NudgingSSTClose(.TRUE.)
      IF (nsstinc > 0) THEN 
        CALL NudgingSSTClose(.TRUE.)
      END IF

      CALL message('','Close nudging data files.')

      CALL set_stream(nudg, lpout=.FALSE., lrerun=.FALSE.)  ! switch output stream off

      lnudge = .FALSE.

      CALL message('','***********************************************')
      CALL message('','')
      CALL message('','Stop nudging mode.')
      CALL message('','')
      CALL message('','***********************************************')


! ==============================================================================

   CASE default

      CALL finish('NudgingInit','Wrong initialization mode.')
      
   END SELECT



 CONTAINS

   SUBROUTINE bc_nudging
     
     USE mo_mpi, ONLY: p_bcast, p_io

     CALL p_bcast (lnudgini,      p_io)
     CALL p_bcast (lnudgcli,      p_io)
     CALL p_bcast (dt_nudg_start, p_io)
     CALL p_bcast (dt_nudg_stop,  p_io)

     CALL p_bcast (ndg_file_vor,  p_io)
     CALL p_bcast (ndg_file_div,  p_io)
     CALL p_bcast (ndg_file_stp,  p_io)

     CALL p_bcast (nsstinc,       p_io)
     CALL p_bcast (nsstoff,       p_io)
     CALL p_bcast (ndg_freez,     p_io)
     CALL p_bcast (ndg_file_sst,  p_io)

     CALL p_bcast (lnudgimp,      p_io)
     CALL p_bcast (lnudgpat,      p_io)
     CALL p_bcast (lnudgfrd,      p_io)

     CALL p_bcast (nudglmin,      p_io)
     CALL p_bcast (nudglmax,      p_io)

     CALL p_bcast (nudgtrun,      p_io)
     CALL p_bcast (nudgsmin,      p_io)
     CALL p_bcast (nudgsmax,      p_io)

     CALL p_bcast (nudgd,         p_io)
     CALL p_bcast (nudgv,         p_io)
     CALL p_bcast (nudgt,         p_io)
     CALL p_bcast (nudgp,         p_io)

     CALL p_bcast (ltintlin,      p_io)
     CALL p_bcast (ldamplin,      p_io)
     CALL p_bcast (nudgdamp,      p_io)
     CALL p_bcast (nudgdsize,     p_io)

     CALL p_bcast (lnudgdbx,      p_io)
     CALL p_bcast (lnudgwobs,     p_io)
     CALL p_bcast (lsite,         p_io)

   END SUBROUTINE bc_nudging

  END SUBROUTINE NudgingInit
!EOX
!EOC
! ==============================================================================
!BOP
  ! !IROUTINE:  StopNudging
  ! !INTERFACE:

  SUBROUTINE StopNudging

    ! !DESCRIPTION: 
    ! close nudging, clean up memory

!EOP
    CALL NudgingInit(NDG_CLEAN_MEM)
    CALL NudgingInit(NDG_CLOSE)

  END SUBROUTINE StopNudging
END MODULE mo_nudging_init
