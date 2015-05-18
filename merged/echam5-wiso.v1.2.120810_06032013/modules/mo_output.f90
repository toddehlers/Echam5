MODULE mo_output

  ! Description:
  !
  ! This module contains the output routines based on CDI version 1.1.0.
  !
  ! For a detailed description of the contained subroutines look on their
  ! header.
  ! 
  ! Most important new features:
  ! - Use of a code table in section 1 which allows the use of 128 tables
  !   each containg 128 different variables.
  ! - Uses the sub center entry (232 for ECHAM). Center is still 98 for
  !   ECMWF
  ! - The century parameter is correctly used. Now, in 360 day mode, 
  !   25599 years can be simulated before an overflow in the date occures.
  ! - Uses correct dates in 365/366 day mode. Allows for clean forecast and
  !   nudging runs.
  ! - The output is now real standard !!!
  !
  ! Authors:
  !
  ! L. Kornblueh,   MPI, December 1998
  ! L. Kornblueh,   MPI, April    1999, added NWP forecast mode
  ! U. Schulzweida, MPI, April    2000, EMOS compatible
  ! I. Kirchner,    MPI, December 2000, time control
  ! A. Rhodin,      MPI, Mai      2001, condensed to one output routine
  ! A. Rhodin,  DWD/MPI, October  2001, NetCDF calls, parallel GRIB encoding
  ! L. Kornblueh,   MPI, October  2001, changed subcenter to 232 
  !                                     (ECMWF assigned)
  ! A. Rhodin,  DWD/MPI, February 2001, bug fixes for parallel mode
  ! U. Schulzweida, MPI, May      2002, blocking (nproma)
  ! R. Johanni, IPP Garching, May-2002, parallel nudging
  ! I. Kirchner,    MPI, August   2002, lpout flag, add backup_output_streams
  ! U. Schulzweida, MPI, February 2003, change codegb5 interface to gribex
  ! A. Rhodin,      DWD, March    2003, bug fix for no_cycles > 1
  ! L. Kornblueh    MPI, April    2003, global communicator specified
  ! U. Schulzweida  MPI, November 2007, changed output to CDI library

  USE mo_kind,         ONLY: dp
  USE mo_time_control, ONLY: ev_putdata, l_putdata, & 
                             get_interval_seconds, next_date, start_date, &
                             get_date_components, get_forecast_hours,     &
                             write_date
  USE mo_netcdf,       ONLY: IO_dim_ids
  USE mo_decomposition,ONLY: ld => local_decomposition,  &
                             gd => global_decomposition

  IMPLICIT NONE

  INTEGER :: gaussianID, spectralID, taxisIDa, taxisIDr
  INTEGER :: belowsurID, surfaceID, hybridID, hybrid_hID
  INTEGER :: instID, modelID
  INTEGER :: local_tableID, nudging_tableID, tracer_tableID, chem_tableID

  INTEGER ::  ksec1(43)

  INTEGER, PARAMETER :: local_table     = 128 !  local code table
  INTEGER, PARAMETER :: nudging_table   = 129 !  nudging code table
  INTEGER, PARAMETER :: tracer_table    = 131 !  tracer code table
  INTEGER, PARAMETER :: chem_table      = 199 !  chemie code table
  INTEGER, PARAMETER :: center_id       =  98 !  identification of centre
  INTEGER            :: model_id              !  model identification
  INTEGER, PARAMETER :: grid_type       = 255 !  grid definition
  INTEGER, PARAMETER :: nflag           = 128 !  flag(see code table 1)
  INTEGER            :: code_parameter        !  data field code 
                                              !  (see code table 2 )

  INCLUDE 'cdi.inc'


  ! reference time of data

  INTEGER            :: year                  !  year of century 
  INTEGER            :: month                 !  month 
  INTEGER            :: day                   !  day 
  INTEGER            :: hour                  !  hour
  INTEGER            :: minute                !  minute
  INTEGER            :: second                !  second

  INTEGER            :: time_unit       =   0 ! unit of time range 
                                              ! (see code table 4)
  INTEGER            :: time_p1         =   0 ! time range 1
  INTEGER            :: time_p2         =   0 ! time range 2
  INTEGER            :: range_flag      =  10 ! time range flag
                                              ! (see code table 5)
  INTEGER            :: century               ! century
  INTEGER, PARAMETER :: subcenter_id    = 232 ! subcenter

  INTEGER            :: forecast_hours  =   0 ! number of forecast hours in
                                              ! NWP mode

  INTEGER, PARAMETER :: closeID =  -1

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE set_output_time

    USE mo_control,      ONLY: lnwp
    USE mo_time_control, ONLY: get_forecast_hours

    IF (lnwp) THEN
       CALL get_date_components(start_date,year,month,day,hour,minute,second)
       forecast_hours = get_forecast_hours()
    ELSE
       CALL get_date_components(next_date,year,month,day,hour,minute,second)
       forecast_hours = 0
    END IF
    century = year/100+1
    year    = MOD(year,100)

  END SUBROUTINE set_output_time

!------------------------------------------------------------------------------
  SUBROUTINE close_output_streams

    USE mo_memory_base, ONLY: ostreams, nstreams, &! output streams
                              nofiles, cfilenames  ! subjob interface
    USE mo_mpi,         ONLY: p_pe, p_io           ! processor id's
    !---------------------------------------------------------------------
    ! Loop over all the output streams and close the associated files, set
    ! unit numbers to closeID.
    !---------------------------------------------------------------------
    INTEGER                     :: i

    nofiles = 0
    !-----------------------------
    ! loop over all output streams
    !-----------------------------
    DO i=1,nstreams
      IF (ostreams(i)% fileID /= closeID) THEN
        IF (p_pe == p_io) CALL streamClose(ostreams(i)% fileID)
        IF (p_pe == p_io) CALL vlistDestroy(ostreams(i)% vlistID)

        WHERE (ostreams% filetype == ostreams(i)% filetype  &
             .AND. ostreams% fileID   == ostreams(i)% fileID)   &
             ostreams% fileID = closeID

        WHERE (ostreams% filetype == ostreams(i)% filetype  &
             .AND. ostreams% fileID   == ostreams(i)% fileID)   &
             ostreams% vlistID = closeID

        ! store subjob command
        nofiles = nofiles + 1
        cfilenames(nofiles) = 'OSTREAM' // &
             TRIM(ostreams(i)%post_suf) // '=' &
             // TRIM(ostreams(i)%filename)
      END IF
    END DO

    !--------------------------------------------
    ! zero fileIDs on all processor elements
    !--------------------------------------------
    ostreams% fileID = closeID

  END SUBROUTINE close_output_streams

!------------------------------------------------------------------------------
  SUBROUTINE backup_output_streams

    USE mo_exception,   ONLY: message
    USE mo_memory_base, ONLY: ostreams, nstreams   ! output streams
    USE mo_mpi,         ONLY: p_pe, p_io           ! processor id's
    !---------------------------------------------------------------------
    ! make a copy of all output streams
    !---------------------------------------------------------------------
    INTEGER         :: i, ilenc
    CHARACTER (512) :: ycopy, my_mess

    ! the command is machine dependent, now only for linux testet
    CHARACTER (512), PARAMETER :: ycopy_cmd = 'cp'     ! backup command

    INTEGER, EXTERNAL :: util_system

    IF (p_pe == p_io) THEN

      stream_loop: DO i=1,nstreams

        IF (ostreams(i)% fileID /= closeID) THEN
          ycopy = TRIM(ycopy_cmd) // ' ' &
               // TRIM(ostreams(i)%filename) // ' ' &
               // TRIM(ostreams(i)%filename) // '.bak'
          ilenc = MIN(LEN_TRIM(ycopy),512)
          WRITE(my_mess,*) 'backup: ',TRIM(ostreams(i)%filename),' <',TRIM(ycopy),'>'
          CALL message('backup_output_streams',my_mess)
          IF (util_system(ycopy(:ilenc)) /= 0) &
               CALL message('backup_output_streams','copy failed')
!!!               CALL finish('backup_output_streams','copy failed')
        END IF

      END DO stream_loop

    END IF

  END SUBROUTINE backup_output_streams

!------------------------------------------------------------------------------
  SUBROUTINE open_output_streams

  ! Loop over all the output streams and open the associated files. Set
  ! unit numbers (file IDs) for all streams associated with a file.

    USE mo_filename,    ONLY: out_expname, NONE    ! experiment name
    USE mo_netcdf, ONLY: global_att ! global attributes
    USE mo_memory_base, ONLY: ostreams, nstreams, &! output streams
                              GRIB, NETCDF, NETCDF64,  &! allowed file types
                              GAUSSIAN, SPECTRAL   ! grid representations
    USE mo_linked_list, ONLY: print_stream,       &! print output stream
                              list_element         ! output stream entry
    USE mo_doctor,      ONLY: nout, nerr           ! standard output,error unit
    USE mo_control,     ONLY: lnwp, nn,           &! forecast mode flag
                              lcolumn 
    USE mo_util_string, ONLY: separator            ! format string (----)
    USE mo_filename,    ONLY: compose_filenames,  &!
                              standard_grib_file, &!
                              find_next_free_unit
    USE mo_exception,   ONLY: finish
    USE mo_mpi,         ONLY: p_pe, p_io, p_bcast

    INTEGER                                :: i ,j, ia  ! loop indices
    CHARACTER (LEN(standard_grib_file)+4)  :: base  !
    CHARACTER (LEN(standard_grib_file)+16) :: file  !
    LOGICAL                                :: first = .TRUE.
    INTEGER                                :: used(255) ! used codes
    INTEGER                                :: iunit ! codefile unit
    INTEGER                                :: status

    ! Derive base of filename
    ! filename is: standard_grib_file[+forecast_hours][suffix][.nc]
    ! output streams with the same 'suffix' use the same file

    CALL compose_filenames
    IF (lnwp) THEN
      CALL set_output_time
      IF (forecast_hours > 744) THEN
        WRITE (nerr,*) 'NWP mode makes no sense for this time range.'
        WRITE (nerr,*) 'Please change to the climate mode.'
        CALL finish('open_output_streams','Run terminated.')
      END IF
    ENDIF
    base = standard_grib_file
    !------------------------------
    ! print output streams (header)
    !------------------------------
    IF (p_pe == p_io) THEN
      WRITE(nout,separator)
      WRITE(nout,*)
      WRITE(nout,'("  open_output_streams:")')
      WRITE(nout,*)
      WRITE(nout,*)' basename=',TRIM(base)
      WRITE(nout,*)
      WRITE(nout,'(4x,a32,1x,a8,1x,a12,1x,a5,1x,a5)')               &
        'file                            ','stream  ','fileID',&
        'lpost','lrest'
      WRITE(nout,*)
    ENDIF
    !--------------------------------------------------------------------
    ! Loop over all the output streams and open the associated files. Set
    ! file IDs for all streams associated with a file.
    !--------------------------------------------------------------------
    ostreams% first = .FALSE.
    DO i = 1, nstreams
      !-------------------------------
      ! Skip if file is already opened
      !-------------------------------
      IF (ostreams(i)% fileID /= closeID) CYCLE
      !------------------------------------------------
      ! Derive filename
      ! Skip if column model is running and GRIB output
      !------------------------------------------------
      SELECT CASE (ostreams(i)% filetype)
      CASE (GRIB)
        IF (lcolumn) CYCLE
        file = TRIM(base)//ostreams(i)% post_suf
        !---------------------
        ! open code table file
        !---------------------
        IF (p_pe == p_io) THEN
          iunit = find_next_free_unit (80, 100)
          OPEN (iunit,file=TRIM(file)//'.codes')
        ENDIF
      CASE (NETCDF)
        file = TRIM(base)//TRIM(ostreams(i)% post_suf)//'.nc'
      CASE (NETCDF64)
        file = TRIM(base)//TRIM(ostreams(i)% post_suf)//'.nc'
      CASE default
        CALL finish('open_output_streams','unknown filetype')
      END SELECT
      !------------------------------------------------------------
      ! loop over all output streams corresponding to the same file
      !------------------------------------------------------------
      used = 0
      DO j = i, nstreams
        IF (ostreams(j)% post_suf /= ostreams(i)% post_suf)   CYCLE
        IF (ostreams(j)% filetype /= ostreams(i)% filetype)   &
          CALL finish('open_output_streams',              &
                      'different file types for same file')
        !----------------------------
        ! check for valid grib codes
        ! print code table
        !----------------------------
        CALL test_codes 
        !----------
        ! Open file
        !----------
        IF (i==j) THEN
          IF (ostreams(i)% lpost) THEN
            ostreams(i)% first = .TRUE.

            IF (p_pe==p_io) THEN
              SELECT CASE (ostreams(i)%filetype)
              CASE (GRIB)
                ostreams(i)%fileID = streamOpenWrite(file, FILETYPE_GRB)
                IF ( ostreams(i)% ztype /= NONE ) THEN
                  CALL streamDefZtype(ostreams(i)%fileID, ostreams(i)% ztype);
                  CALL streamDefZlevel(ostreams(i)%fileID, 0);
                END IF
              CASE (NETCDF)
                ostreams(i)%fileID = streamOpenWrite(file, FILETYPE_NC)
              CASE (NETCDF64)
                ostreams(i)%fileID = streamOpenWrite(file, FILETYPE_NC2)
              CASE default
                CALL finish('open_output_streams','unknown filetype')
              END SELECT

              IF ( ostreams(i)%fileID .LT. 0 ) THEN
                WRITE(0,*) cdiStringError(ostreams(i)%fileID)
                CALL finish ('open_output_streams', 'Open failed on '//TRIM(file))
              END IF

              !--------------------------------------------------------------
              ! due to performance reasons switch off fill mode (netCDF only)
              !--------------------------------------------------------------
              ! CALL nf(nf_set_fill(stream%fileID, nf_nofill, old_mode))

              ostreams(i)%vlistID = vlistCreate()
              ostreams(i)%timestep = 0

              status = vlistDefAttTxt(ostreams(i)%vlistID, CDI_GLOBAL, 'title', len(TRIM(out_expname)), TRIM(out_expname))
              DO ia = 1, SIZE(global_att)
                IF (global_att(ia)% name /= '') &
                  status = vlistDefAttTxt(ostreams(i)%vlistID, CDI_GLOBAL, &
                  global_att(ia)%name, len(TRIM(global_att(ia)%text)), TRIM(global_att(ia)%text) )
              END DO
              status = vlistDefAttInt(ostreams(i)%vlistID, CDI_GLOBAL, 'truncation', 1, nn)

              IF ( ostreams(i)%filetype == GRIB ) THEN
                CALL vlistDefTaxis(ostreams(i)%vlistID, taxisIDa)
              ELSE
                CALL vlistDefTaxis(ostreams(i)%vlistID, taxisIDr)
              END IF
            ENDIF

            CALL p_bcast (ostreams(i)%fileID, p_io)
            CALL p_bcast (ostreams(i)%vlistID, p_io)
            CALL p_bcast (ostreams(i)%timestep, p_io)
            ostreams(i)% filename = file
          ELSE
            EXIT
          END IF
        END IF
        !---------------------------------------------------
        ! set file IDs of all associated output streams
        !---------------------------------------------------
        IF(.NOT. ostreams(i)% lpost) CYCLE

        ostreams(j)%fileID   = ostreams(i)%fileID
        ostreams(j)%vlistID  = ostreams(i)%vlistID
        ostreams(j)%timestep = ostreams(i)%timestep
 
        IF (p_pe==p_io) THEN

          CALL addStreamToVlist(ostreams(j), ostreams(j)%vlistID)

          WRITE(nout,'(4x,a32,1x,a8,1x,i12,3x,l1,5x,l1)')      &
            file, ostreams(j)% name, ostreams(j)% fileID, &
            ostreams(j)% lpost, ostreams(j)% lrerun
        ENDIF
      END DO

      IF (p_pe==p_io.AND.ostreams(i)%lpost) &
           CALL streamDefVlist(ostreams(i)%fileID, ostreams(i)%vlistID)

      IF (p_pe==p_io) THEN
        !----------------------
        ! close code-table file
        !----------------------
        IF (ostreams(i)% filetype == GRIB) THEN
          IF (SUM(used) > 0) THEN
            CLOSE (iunit)
          ELSE
            CLOSE (iunit, status='DELETE')
          ENDIF
        ENDIF
      ENDIF
    END DO

    !---------------------
    ! print output streams
    !---------------------
    IF (p_pe==p_io) THEN
      WRITE(nout,*)
      IF (first) THEN
        first = .FALSE.
        DO i = 1, nstreams
          CALL print_stream (ostreams (i))
        END DO
      ENDIF
    ENDIF
  CONTAINS
!..............................................................................
    !
    ! check grib codes, write code table
    !
    SUBROUTINE test_codes
      TYPE (list_element) ,POINTER :: next
      LOGICAL :: lpost, lrerun
      INTEGER :: newcode, minl(1)
      !
      ! flags for postprocessing, restart, initialisation
      !
      lpost  = .FALSE.
!     linit  = .FALSE.
      lrerun = .FALSE.
      !
      ! loop over stream entries
      !
      next => ostreams(j)% first_list_element
      DO
        IF (.NOT.ASSOCIATED(next)) EXIT
        !
        ! update flags
        !
        lrerun = lrerun .OR. next% field% info% lrerun
        lpost  = lpost .OR. next% field% info% lpost
        !
        ! check if fields can be written
        !
        IF (next% field% info% lpost) THEN
          !
          ! check for NETCDF output
          !
          IF (ostreams(j)% filetype == NETCDF) THEN
            SELECT CASE (next% field% info% repr)
            CASE (GAUSSIAN)
              SELECT CASE (next% field% info% ndim)
              CASE (2,3,4)
              CASE default
                next% field% info% lpost = .FALSE.
              END SELECT
            CASE (SPECTRAL)
              IF (lcolumn) next% field% info% lpost = .FALSE.
            CASE default
              next% field% info% lpost = .FALSE.
            END SELECT
          ENDIF
          !
          ! check for GRIB output
          !
          IF (ostreams(j)% filetype == GRIB) THEN
            !
            ! check dimensions
            !
            SELECT CASE (next% field% info% repr)
            CASE (GAUSSIAN)
              SELECT CASE (next% field% info% ndim)
              CASE (2,3)
              CASE default
                next% field% info% lpost = .FALSE.
              END SELECT
            END SELECT   
            !
            ! check codes for GRIB output
            !
            IF (next% field% info% gribcode == 0 &
            .OR.next% field% info% gribcode >255 ) THEN
              !
              ! no valid code
              !
              IF (p_pe==p_io) &
                WRITE(nout,*) 'no gribcode for ', next% field% info% name
            ELSE IF (next% field% info% gribcode < 0) THEN
              !
              ! automatic code determination
              !
              minl    = MINLOC(used)
              newcode = minl(1)
              IF (used(newcode)==0) THEN
                next% field% info% gribcode = newcode
              ELSE
                IF (p_pe==p_io) &
                  WRITE(nout,*) 'no gribcode for ', next% field% info% name
                  next% field% info% gribcode = 0
              ENDIF
            ENDIF
            !
            ! check for codes used twice
            !
            newcode = next% field% info% gribcode
            IF (newcode > 0) THEN
              IF (used(newcode)/=0) THEN
                IF (p_pe==p_io) &
                  WRITE(nout,*)'gribcode used twice for ',next% field% info% name
                next% field% info% gribcode = 0
              ELSE
                used (newcode) = used (newcode) + 1
              ENDIF
            ENDIF
            !
            ! write code table
            !
            IF (p_pe==p_io) THEN
              IF (next% field% info% gribcode /= 0) THEN
                WRITE(iunit,'(i4,i4,1x,a,1x,f7.2,1x,f7.2,1x,a," [",a,"]")') &
                  next% field% info% gribcode,      &
                  next% field% info% klev,          &
                  next% field% info% name,          &
                  0., 1.,                     &
                  TRIM(next% field% info% longname),&
                  TRIM(next% field% info% units)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        next => next% next_list_element
      END DO
      !
      ! update flags
      !
!      ostreams(j)% lpost  = ostreams(j)% lpost  .AND. lpost
!      ostreams(j)% lrerun = ostreams(j)% lrerun .AND. lrerun
    END SUBROUTINE test_codes

  END SUBROUTINE open_output_streams

!------------------------------------------------------------------------------
  SUBROUTINE addStreamToVlist(stream, vlistID)

    USE mo_exception,    ONLY: finish
    USE mo_linked_list,  ONLY: t_stream,         & ! output stream data type
                               memory_info,      & ! stream entry data type
                               list_element,     & ! linked list element type
                               GAUSSIAN,         & ! Gaussian grid  indicator
                               SPECTRAL            ! spectral grid  indicator

    USE mo_netCDF,       ONLY: BELOWSUR, SURFACE, HYBRID, HYBRID_H
    USE mo_memory_base,  ONLY: GRIB

    TYPE (t_stream) ,INTENT(in) :: stream ! output stream
    INTEGER,INTENT(inout) :: vlistID
      
    TYPE (list_element) ,POINTER    :: le
    TYPE (list_element) ,TARGET     :: first
    TYPE (memory_info)  ,POINTER    :: info

    INTEGER                :: varID, gridID, zaxisID
    INTEGER                :: n
    INTEGER                :: prec          ! float prec. (default)
    CHARACTER(len=5)       :: axis
    CHARACTER(len=32)      :: grid_type

    first%next_list_element => stream%first_list_element
    !---------------------------------------
    ! 1st loop, define additional dimensions
    !---------------------------------------
    le => first
    DO ! loop over elements in linked list
      le => le%next_list_element
      IF (.NOT.ASSOCIATED(le)) EXIT
      info => le%field%info
      IF (.NOT. info%lpost)      CYCLE ! skip if lpost flag not set
      IF (.NOT.(info%repr == GAUSSIAN .OR. info%repr == SPECTRAL)) CYCLE
      !-----------------------------------------------------------
      ! Only a standard 2D field (nlon,ngl) and a
      ! standard 3D field (nlon,...,ngl) 
      ! and spectral representations are allowed for netcdf/GRIB
      ! output. For all other fields, info%lpost is set to .FALSE.
      ! surface pressure is written in any case. 
      !-----------------------------------------------------------
      n = info%ndim
      IF (info%repr == SPECTRAL .AND. info% gdim(n) == 1) n=n-1
      !----------------------------------------
      ! write message for grid types not written
      !-----------------------------------------
      IF (.NOT.info%lpost) THEN
        PRINT *,'   ',TRIM(info%name), ' is non-standard: info%gdim_* = ' &
          ,info%gdim(1), info%gdim(2), info%gdim(3), info%gdim(4)
        CYCLE
      ENDIF
    END DO
    !------------------------------------------
    ! 2nd loop, define variables and attributes
    !------------------------------------------
    le => first
    DO ! loop over elements in linked list
      le => le%next_list_element
      IF (.NOT.ASSOCIATED(le)) EXIT
      info => le%field%info
      IF (.NOT. info%lpost)      CYCLE ! skip if lpost flag not set

      n = info% ndim
      SELECT CASE (info%repr)
      CASE (GAUSSIAN)
        info%gridID = gaussianID
        grid_type = 'gaussian'
        IF ((n==4)) THEN
          CALL finish('addStreamToVlist', '5-dimensional arrays unsupported!');
        ELSE IF ((n==3)) THEN
          !------------------------------
          ! 3d data, transpose dimensions
          !------------------------------
          axis       = 'tzyx'
        ELSE
          !-------------
          ! regular data
          !-------------
          axis      = 'tyx'
        ENDIF
      CASE (SPECTRAL)
        info%gridID = spectralID
        IF (info% gdim(1) == 1) THEN
          n=2
          axis      = 't--'
        ELSE
          axis      = 'tz--'
        ENDIF
        grid_type = 'spectral, triangular truncation'
      CASE default
        CYCLE
      END SELECT

      IF ( info%levelindx == SURFACE ) THEN
        info%zaxisID = surfaceID
      ELSE IF ( info%levelindx == HYBRID ) THEN
        info%zaxisID = hybridID
      ELSE IF ( info%levelindx == HYBRID_H ) THEN
        info%zaxisID = hybrid_hID
      ELSE IF ( info%levelindx == BELOWSUR ) THEN
        info%zaxisID = belowsurID
      END IF

      gridID  = info%gridID
      zaxisID = info%zaxisID

      IF ( gridID  == -1 ) THEN
        WRITE(0,*) 'GRID undefined for ', info%name
        CALL finish('addStreamToVlist', 'GRID definition missing')
      END IF
      IF ( zaxisID == -1 ) THEN
        WRITE(0,*) 'ZAXIS undefined for ', info%name
        CALL finish('addStreamToVlist', 'ZAXIS definition missing')
      END IF

      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE)
      info%IO_var_stid = varID ! store for later use

      prec = DATATYPE_FLT32
      IF ( stream% filetype == GRIB ) THEN
        prec =  info%gribbits
      ELSE
        IF (info%gribbits > 32) prec = DATATYPE_FLT64
      END IF
        
      CALL vlistDefVarDatatype(vlistID, varID, prec)
      CALL vlistDefVarName(vlistID, varID, info%name)

      IF (info% longname/='')  CALL vlistDefVarLongname(vlistID, varID, info%longname)
      IF (info% units/='')     CALL vlistDefVarUnits(vlistID, varID, info%units)
      IF (info% gribcode > 0)  CALL vlistDefVarCode(vlistID, varID, info%gribcode)

      IF (info% gribtable > 0) THEN
        SELECT CASE (info%gribtable)
        CASE(local_table)
          CALL vlistDefVarTable(vlistID, varID, local_tableID)
        CASE(nudging_table)
          CALL vlistDefVarTable(vlistID, varID, nudging_tableID)
        CASE(tracer_table)
          CALL vlistDefVarTable(vlistID, varID, tracer_tableID)
        CASE(chem_table)
          CALL vlistDefVarTable(vlistID, varID, chem_tableID)
        END SELECT
      END IF

      CALL vlistDefVarInstitut(vlistID, varID, instID) 
      CALL vlistDefVarModel(vlistID, varID, modelID) 
      
      !------------------------------------------
      ! print tracer attributes into netcdf file:
      !------------------------------------------
      IF (info% tracidx > 0) CALL write_tracer_header(info% tracidx, vlistID, varID)

    END DO

  END SUBROUTINE addStreamToVlist

  SUBROUTINE write_tracer_header(tracidx, vlistID, varID)
    USE mo_tracdef, ONLY: trlist               ! tracer info variable
    INTEGER ,INTENT(in) :: tracidx             ! tracer index
    INTEGER ,INTENT(in) :: vlistID, varID
    INTEGER             :: status

    status = vlistDefAttInt(vlistID, varID, 'index',      1, tracidx)
    status = vlistDefAttFlt(vlistID, varID, 'molar_mass', 1, trlist%ti(tracidx)%moleweight)
    status = vlistDefAttFlt(vlistID, varID, 'Henry',      1, trlist%ti(tracidx)%henry)
    status = vlistDefAttFlt(vlistID, varID, 'dryreac',    1, trlist%ti(tracidx)%dryreac)
    status = vlistDefAttInt(vlistID, varID, 'ndrydep',    1, trlist%ti(tracidx)%ndrydep)
    status = vlistDefAttInt(vlistID, varID, 'ntran',      1, trlist%ti(tracidx)%ntran)
    status = vlistDefAttInt(vlistID, varID, 'nvdiff',     1, trlist%ti(tracidx)%nvdiff)
    status = vlistDefAttInt(vlistID, varID, 'nconv',      1, trlist%ti(tracidx)%nconv)
    status = vlistDefAttInt(vlistID, varID, 'nwetdep',    1, trlist%ti(tracidx)%nwetdep)
    status = vlistDefAttInt(vlistID, varID, 'nsoluble',   1, trlist%ti(tracidx)%nsoluble)
    status = vlistDefAttInt(vlistID, varID, 'nphase',     1, trlist%ti(tracidx)%nphase)
    status = vlistDefAttInt(vlistID, varID, 'mode',       1, trlist%ti(tracidx)%mode)

  END SUBROUTINE write_tracer_header

!------------------------------------------------------------------------------
  SUBROUTINE init_output

    USE mo_constants, ONLY: api
    USE mo_control,   ONLY: nvclev, nlon, ngl, nhgl, nn, lnmi, lnudge
    USE mo_control,   ONLY: vct
    USE mo_gaussgrid, ONLY: gl_gmu
    USE mo_doctor,    ONLY: nvers
    USE mo_exception, ONLY: message
    USE mo_netCDF,    ONLY: BELOWSUR
    USE mo_netCDF,    ONLY: io_dim_ids,       & ! dimension table
                            io_dim              ! table entry data type
    USE mo_time_control,   ONLY: start_date,       & ! reference for time axis
                                 get_date_components ! split date into components

    TYPE (io_dim) ,POINTER :: pdim
    REAL(dp) ,POINTER :: levels(:), xvals(:), yvals(:)
    INTEGER :: i
    INTEGER :: year, month, day, hour, minute, second ! date/time variables
    INTEGER :: idate, itime
    CHARACTER (32)  :: ymodel

    ! set model version

    model_id = nvers + 10
    WRITE(ymodel,'(a,f3.1)') 'ECHAM', REAL(nvers,dp)/10_dp

    ! set fixed values in GRIB block 1, common to all GRIB sets

    ! to set standard output to analysis mode lnwp must be .false. and
    ! lanalysis .true.. Later with working adiabatic NMI this must be 
    ! changed to range_flag = 1. range_flag = 0 is default for output
    ! of the OI or 3DVAR analysis, this are given external. In case a 
    ! 4DVAR is running inside ECHAM this value must be adjusted to 
    ! range_flag = 0, as it is now.
    ! For usual climate mode runs range_flag is set to 10, which means
    ! The time given in the GRIB header in section 1 is the valid time
    ! and time_p1 and time_p2 do not mean anything. 

    IF (lnmi .AND. .NOT.lnudge) range_flag = 0    

    ksec1(18) = range_flag

    IF ( nvclev .GT. 127 ) CALL message('GRIB', 'VCT not defined for more than 126 levels!')

    gaussianID = gridCreate(GRID_GAUSSIAN, nlon*ngl)
    CALL gridDefXsize(gaussianID, nlon)
    CALL gridDefYsize(gaussianID, ngl)
    ALLOCATE(xvals(nlon))
    ALLOCATE(yvals(ngl))
    DO i = 1, nlon
      xvals(i) = (i-1)*360.0_dp/nlon
    END DO
    DO i = 1, nhgl
      yvals(i) = ASIN(gl_gmu(i))*180.0_dp/api
      yvals(ngl-i+1) = -yvals(i)
    END DO
    CALL gridDefXvals(gaussianID, xvals)
    CALL gridDefYvals(gaussianID, yvals)
    DEALLOCATE(xvals)
    DEALLOCATE(yvals)

    spectralID = gridCreate(GRID_SPECTRAL, (nn+1)*(nn+2))
    CALL gridDefTrunc(spectralID, nn)

    surfaceID  = zaxisCreate(ZAXIS_SURFACE, 1)
    belowsurID = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, 5)
    hybridID   = zaxisCreate(ZAXIS_HYBRID, nvclev-1)
    hybrid_hID = zaxisCreate(ZAXIS_HYBRID_HALF, nvclev)

    ALLOCATE(levels(1))
    levels(1) = 0.0_dp
    CALL zaxisDefLevels(surfaceID, levels)
    DEALLOCATE(levels)

    pdim => IO_dim_ids(BELOWSUR)
    IF (ASSOCIATED (pdim% value)) THEN
      CALL zaxisDefName(belowsurID, pdim%dim_name)
      CALL zaxisDefLevels(belowsurID, pdim%value)
    END IF

    ALLOCATE(levels(nvclev))
    DO i = 1, nvclev
      levels(i) = i
    END DO
    CALL zaxisDefLevels(hybridID, levels)
    CALL zaxisDefLevels(hybrid_hID, levels)
    DEALLOCATE(levels)

    CALL zaxisDefVct(hybridID,   2*nvclev, vct(1:2*nvclev))
    CALL zaxisDefVct(hybrid_hID, 2*nvclev, vct(1:2*nvclev))

    taxisIDa = taxisCreate(TAXIS_ABSOLUTE);
    taxisIDr = taxisCreate(TAXIS_RELATIVE);

    CALL get_date_components(start_date,year,month,day,hour,minute,second)
    idate = ABS(year)*10000+month*100+day
    if ( year < 0 ) idate = -idate
    itime = hour*100+minute

    CALL taxisDefRdate(taxisIDr, idate)
    CALL taxisDefRtime(taxisIDr, itime);
    CALL taxisDefTunit(taxisIDr, TUNIT_DAY);

    instID  = institutDef(center_id, subcenter_id, "MPIMET", "Max-Planck-Institute for Meteorology");
    modelID = modelDef(instID,  model_id, ymodel);
    local_tableID   = tableDef(modelID, local_table,   "echam5");
    nudging_tableID = tableDef(modelID, nudging_table, "echam5");
    tracer_tableID  = tableDef(modelID, tracer_table,  "echam5");
    chem_tableID    = tableDef(modelID, chem_table,    "echam5");

  END SUBROUTINE init_output

!------------------------------------------------------------------------------
  SUBROUTINE out_streams
    !-------------------------------------------------
    ! loop over all output streams
    ! write if flag lpost is set and
    !       if output time (flag L_PUTDATA) is reached
    !------------------------------------------------- 
    USE mo_memory_base,  ONLY: ostreams, &! output streams
                               nstreams, &! number of active streams
                               NETCDF     ! file type flag value
    USE mo_mpi,          ONLY: p_pe, p_io ! processor element indices
    USE mo_filename,     ONLY: out_filetype, GRIB, NETCDF 
    USE mo_filename,     ONLY: out_ztype, SZIP

    INTEGER       :: i,j
    LOGICAL       :: time_written, write_info

    write_info = .TRUE.

    !-----------------------------------------------
    ! loop over all output streams
    ! pick up first stream associated with each file
    !-----------------------------------------------
    DO i=1,nstreams
      IF (.NOT. ostreams(i)% first) CYCLE
      time_written = .FALSE.
      !-----------------------------------------------
      ! loop over all streams associated with the file
      !-----------------------------------------------
      DO j=i,nstreams
        IF (ostreams(j)%filetype == ostreams(i)%filetype .AND. &
            ostreams(j)%fileID   == ostreams(i)%fileID   ) THEN 
          !---------------------------
          ! check condition for output
          !---------------------------
          IF(l_putdata (ostreams(j)% post_idx) &
             .AND.      ostreams(j)% lpost     &
             .AND.      ostreams(j)% lpout     ) THEN
            !
            IF (Write_info) THEN
              SELECT CASE (out_filetype)
              CASE (NETCDF)
                CALL write_date(next_date,'Write netCDF output for : ')
              CASE (GRIB)
                IF ( out_ztype == SZIP ) THEN
                  CALL write_date(next_date,'Write GRIB1/SZIP output for : ')
                ELSE
                  CALL write_date(next_date,'Write GRIB1 output for : ')
                END IF
              END SELECT
              write_info = .FALSE.
            ENDIF

            !------------------------------------
            ! increase time slice in NETCDF files
            !------------------------------------
            IF ( (p_pe==p_io) .AND.  .NOT. time_written) THEN
              CALL write_time_to_stream(ostreams(i))
              time_written = .TRUE.
            ENDIF

            !----------------
            ! write variables
            !----------------
            CALL out_stream (ostreams(j))
          ENDIF
        END IF
      END DO
    END DO

  END SUBROUTINE out_streams

!------------------------------------------------------------------------------
  SUBROUTINE out_stream (stream)
  !
  ! Description:
  !
  ! Control postprocessing of output fields.
  ! Generic routine for all streams
  !
    !
    ! Modules used
    !   
    USE mo_exception,     ONLY: finish
    USE mo_control,       ONLY: lnwp, lcolumn
    USE mo_linked_list,   ONLY: list_element, memory_info, t_stream, &
                                GAUSSIAN, SPECTRAL
    USE mo_memory_base,   ONLY: GRIB, NETCDF         ! allowed file types
    USE mo_transpose,     ONLY: gather_gp, gather_sp
    USE mo_mpi,           ONLY: p_pe, p_io
    !
    ! Argument
    !
    TYPE (t_stream) ,INTENT(in) :: stream
    !
    !  Local scalars: 
    !
    REAL(dp)          :: quot
    INTEGER           :: gridtype, no
    CHARACTER(len=16) :: yname
    INTEGER           :: nglat     ! number of latitudes on this pe
    INTEGER           :: nglon     ! number of longitudes on this pe
    INTEGER           :: level_idx ! index to level definition table
    !
    ! variables of derived type used in linked list
    !
    TYPE (memory_info)  ,POINTER :: info
    TYPE (list_element) ,POINTER :: element
    TYPE (list_element) ,TARGET  :: start
    !
    !  Local arrays: 
    !
    REAL(dp),POINTER     :: ptr4d (:,:,:,:) ! field distributed over processors
    REAL(dp),POINTER     :: z4d   (:,:,:,:) ! field gathered on I/O processor
    !
    !  External subroutines 
    !
    EXTERNAL pbwrite
    !
    nglon = ld%nglon
    nglat = ld%nglat
    !
    ! Definition of grib blocks
    !
    CALL set_output_time

    IF (lnwp) THEN
      ksec1(15) = time_unit
      ksec1(16) = forecast_hours
      ksec1(18) = range_flag
    END IF
    !
    ! Loop over all fields in linked list
    !
    element => start
    element% next_list_element => stream% first_list_element    
    DO
      element => element% next_list_element
      IF(.NOT.ASSOCIATED(element)) EXIT
      !-----------------------------------------------------
      ! retrieve information from actual linked list element
      !-----------------------------------------------------
      info           => element% field% info 
      ptr4d          => element% field% ptr(:,:,:,:)
      yname          = info% name
      code_parameter = info% gribcode
      gridtype       = info% repr
      level_idx      = info% levelindx
      !------------------
      ! skip this field ?
      !------------------
      IF ( .NOT. info% lpost           ) CYCLE
      IF ( yname                == ' ' ) CYCLE
      IF ( stream% filetype     == GRIB &
           .AND. code_parameter <= 0   ) CYCLE
      IF ( level_idx < 1 .OR. level_idx > SIZE(IO_dim_ids) ) CYCLE
      !------------------------------------------
      ! rescale field if accumulation flag is set
      !------------------------------------------
      IF (info% laccu) THEN
        no = get_interval_seconds(ev_putdata(stream%post_idx))
        IF (no > 0) THEN
          quot = 1.0_dp/no
        ELSE
          quot = 1.0_dp
        END IF
        ptr4d(:,:,:,:) = ptr4d(:,:,:,:) * quot
      ENDIF

      !----------------------------------------------------
      ! Allocate temporary global array on output processor
      ! Gather field from other processors
      !----------------------------------------------------
      NULLIFY  (z4d)
      IF (.NOT.lcolumn) THEN
        IF (p_pe==p_io) &
             ALLOCATE (z4d (info%gdim(1),info%gdim(2),info%gdim(3),info%gdim(4)) )
        SELECT CASE (gridtype)
        CASE (GAUSSIAN)
          CALL gather_gp (z4d, ptr4d, gd)
        CASE (SPECTRAL)
          CALL gather_sp (z4d, ptr4d, gd)
        CASE default
          CALL finish('out_stream','unknown grid type')
        END SELECT
      ENDIF
      !----------------------
      ! Write data
      !----------------------
      IF (p_pe==p_io) THEN
        IF (lcolumn) THEN
          CALL write_var (info, stream, ptr4d(:,:,:,:))
        ELSE
          CALL write_var (info, stream, z4d(:,:,:,:))
        ENDIF
      END IF

      !-------------------------------------------------
      ! reset field if accumulation or reset flag is set
      !-------------------------------------------------
      IF (info% laccu .OR. info% reset /= 0._dp) THEN
        ptr4d(:,:,:,:) = info% reset
      ENDIF
      !-----------------------------------
      ! Deallocate temporary global arrays
      !-----------------------------------
      IF (ASSOCIATED (z4d))    DEALLOCATE (z4d)
    END DO

!    IF ( (p_pe==p_io) .AND. stream% filetype==NETCDF) &
!      CALL write_end_netcdfstream(stream)

  END SUBROUTINE out_stream

!!$  SUBROUTINE write_end_netcdfstream (stream)
!!$    TYPE (t_stream) ,INTENT(in) :: stream    ! output stream description
!!$
!!$    INTEGER :: ncid                     ! NetCDF file IDs
!!$
!!$    ncid = stream%fileID
!!$
!!$    CALL nf(nf_sync(ncid)) ! write buffer to file
!!$
!!$  END SUBROUTINE write_end_netcdfstream

!------------------------------------------------------------------------------
  SUBROUTINE cleanup_output
    !
    ! Deallocate module variables
    !

    CALL gridDestroy(gaussianID)
    CALL gridDestroy(spectralID)

    CALL zaxisDestroy(surfaceID)
    CALL zaxisDestroy(belowsurID)
    CALL zaxisDestroy(hybridID)
    CALL zaxisDestroy(hybrid_hID)

  END SUBROUTINE cleanup_output

!------------------------------------------------------------------------------
  SUBROUTINE write_var (info, stream, xzy)

    USE mo_exception,    ONLY: finish
    USE mo_linked_list,  ONLY: t_stream,         & ! output stream data type
                               memory_info,      & ! stream entry data type
                               GAUSSIAN,         & ! Gaussian grid  indicator
                               SPECTRAL            ! spectral grid  indicator

    TYPE (memory_info) ,INTENT(in) :: info         ! field description
    TYPE (t_stream)    ,INTENT(in) :: stream       ! output stream description
    REAL(dp)           ,INTENT(in) :: xzy(:,:,:,:) ! output field

    REAL(dp)    ,ALLOCATABLE :: xyz(:,:,:,:) ! transposed field
    INTEGER :: jy, jz                    ! indices used for transposition
    INTEGER :: fileID                    ! File ID
    INTEGER :: varID                     ! Variable ID
    INTEGER :: n                         ! rank of field to write


    fileID  = stream%fileID
    n     = info% ndim
    !-----------------------------------------------------------
    ! variable ID may be invalid for no_cycle > 0 (rerun cycle).
    !  request ID from NetCDF file in this case. 
    !-----------------------------------------------------------
    varID   = info%IO_var_stid

    !----------------------------------------
    ! write 3D,2D Gaussian or spectral fields
    !----------------------------------------
    SELECT CASE (info%repr)
    CASE (GAUSSIAN)

      SELECT CASE (n)
      CASE (3)
        !-------------------------------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,lat) - ONLY FOR single netcdf
        ! output. Parallel netcdf data hasn't been 'gather'd, so is still
        ! scrambled and needs to be 'reorder'd and unpacked -  but the 
        ! COARDS convention for netcdf requires (lon,lat,lev). Applies
        ! to /all/ fields here. 
        !-------------------------------------------------------------------

        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,3),SIZE(xzy,2),1))
        FORALL (jy=1:SIZE(xzy,3),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,1)=xzy(:,jz,jy,1) ! switch lat and lev
        END FORALL
        CALL streamWriteVar(fileID, varID, xyz(:,:,:,1), 0)
        DEALLOCATE (xyz)

      CASE (4)
        !-----------------------------------------------
        ! The array xzy is sorted xzy(lon,lev,?,lat) but
        ! we require (lon,lat,lev,?).
        !-----------------------------------------------

        ALLOCATE (xyz (SIZE(xzy,1),SIZE(xzy,4),SIZE(xzy,2),SIZE(xzy,3)))
        FORALL (jy=1:SIZE(xzy,4),jz=1:SIZE(xzy,2))
          xyz(:,jy,jz,:)=xzy(:,jz,:,jy) ! switch lat and lev
        END FORALL

        CALL finish('4-dimensional arrays unsupported for output!')

        DEALLOCATE (xyz)

      CASE default

        CALL streamWriteVar(fileID, varID, xzy(:,:,1,1), 0)

      END SELECT


    CASE (SPECTRAL)
      IF (info% gdim(1) == 1) n=n-1

      SELECT CASE (n)
      CASE (3)

        ALLOCATE (xyz (SIZE(xzy,2),SIZE(xzy,3),SIZE(xzy,1),1))
        FORALL (jz=1:SIZE(xzy,1))
          xyz(:,:,jz,1)=xzy(jz,:,:,1) ! switch lat and lev
        END FORALL

        CALL streamWriteVar(fileID, varID, xyz(:,:,:,1), 0)

        DEALLOCATE (xyz)

      CASE (2)
 
        CALL streamWriteVar(fileID, varID, xzy(1,:,:,1), 0)

      END SELECT

    END SELECT

  END SUBROUTINE write_var

!------------------------------------------------------------------------------
  SUBROUTINE write_time_to_stream (stream)

    USE mo_time_control,    ONLY: next_date, get_date_components
    USE mo_time_conversion, ONLY: TC_get
    USE mo_linked_list,     ONLY: t_stream ! output stream data type


    TYPE (t_stream)      ,INTENT(inout) :: stream    ! output stream description

    INTEGER :: start_day, start_sec, present_day, present_sec
    INTEGER :: fileID                                 ! File ID
    INTEGER :: year, month, day, hour, minute, second ! date/time variables
    REAL(dp):: yyyymmdd
    INTEGER :: idate, itime, iret

    fileID = stream%fileID
    !-----------------------------------------
    ! convert time_days format into 2 integers
    !-----------------------------------------
    CALL TC_get(start_date,    start_day,   start_sec)
    CALL TC_get(next_date,   present_day, present_sec)

    !--------------------------------
    ! write alternative time yyyymmdd
    !--------------------------------
    CALL get_date_components(next_date,year,month,day,hour,minute,second)
    yyyymmdd = ABS(year)*10000+month*100+day&
             + (hour*3600+minute*60+second)/86400._dp
    IF (year<0) yyyymmdd = -yyyymmdd

    idate = ABS(year)*10000+month*100+day
    if ( year < 0 ) idate = -idate
    itime = hour*100+minute

    CALL taxisDefVdate(taxisIDa, idate)
    CALL taxisDefVtime(taxisIDa, itime)

    CALL taxisDefVdate(taxisIDr, idate)
    CALL taxisDefVtime(taxisIDr, itime)

    iret = streamDefTimestep(fileID, stream%timestep)
    stream%timestep = stream%timestep + 1

  END SUBROUTINE write_time_to_stream

END MODULE mo_output
