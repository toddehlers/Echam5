#ifdef __SX__
#define VECTOR 1
#endif

SUBROUTINE fftd

#ifdef FFT991
  USE mo_kind,          ONLY: dp
  USE mo_fft991,        ONLY: fft991cy
#elif MKL_DFT
  USE mo_mkl_dft,       ONLY: mkl_dft
#elif ESSL_DFT
  USE mo_essl_dft,      ONLY: essl_dft
#else
  USE mo_fft992,        ONLY: fft992
#endif
  USE mo_buffer_fft,    ONLY: fftz, nvar
  USE mo_decomposition, ONLY: dc => local_decomposition
  USE mo_exception,     ONLY: finish
#if defined (VECTOR) && defined (_OPENMP)
  USE omp_lib,          ONLY: omp_get_thread_num, &
                              omp_get_max_threads
#endif

  IMPLICIT NONE

  INTEGER :: inc, isign
  INTEGER :: nlon, nlp2, nlev, nlat
#ifdef FFT991
  REAL(dp) :: zwork((dc%nlon+2) * dc%nflevp1 * dc%nflat * nvar)
#elif defined (VECTOR) && defined (_OPENMP)
  INTEGER :: ivar, tid, nthreads, chunk, rest
  INTEGER, ALLOCATABLE, SAVE :: istart(:), icount(:)
#elif defined (MKL_DFT) || defined (ESSL_DFT)
#else
  INTEGER :: ilat, ivar
#endif
  LOGICAL :: col_1d

!-- 2. Direct Fourier transforms

!-- 2.1 Set constants

  inc    = 1
  isign  = -1
  nlon   = dc% nlon
  nlp2   = nlon + 2
  nlev   = dc% nflevp1
  nlat   = dc% nflat
  col_1d = dc% col_1d

!-- 2.2 fft(vo, d, t, alps, u, v, dtl, dtm, dalpsl, dalpsm)

  IF (.NOT.col_1d) THEN
#ifdef FFT991

    CALL fft991cy(fftz,zwork,inc,nlp2,nlon,nvar*nlev*nlat,isign)

#else

#ifndef _OPENMP

#ifdef VECTOR

    CALL fft992(fftz,inc,nlp2,nlon,nvar*nlev*nlat,isign)

#elif MKL_DFT

    CALL mkl_dft(fftz,inc,nlp2,nlon,nvar*nlev*nlat,isign)

#elif ESSL_DFT

    CALL essl_dft(fftz,nlon,nlp2,nvar*nlev*nlat,isign)

#else
    
    DO ivar = 1, nvar
      DO ilat = 1, nlat
        CALL fft992(fftz(1,1,ilat,ivar),inc,nlp2,nlon,nlev,isign)
      ENDDO
    ENDDO

#endif

#else 

#ifdef VECTOR

    IF (.NOT. ALLOCATED(istart)) THEN
      nthreads = omp_get_max_threads()
      ALLOCATE(istart(0:nthreads), icount(0:nthreads))
      istart(0) = 1
      DO tid = 0, nthreads-1
        chunk = nlat/nthreads 
        rest  = MOD(nlat, nthreads)
        if (tid < rest) chunk = chunk+1
        icount(tid) = chunk
        istart(tid+1) = istart(tid)+chunk
      ENDDO
    ENDIF

!$OMP PARALLEL PRIVATE(tid)
    tid = omp_get_thread_num()
    DO ivar = 1, nvar
        CALL fft992(fftz(1,1,istart(tid),ivar),inc,nlp2,nlon,nlev*icount(tid),isign)
    END DO
!$OMP END PARALLEL

#elif MKL_DFT

!$OMP PARALLEL PRIVATE(ilat)
    DO ivar = 1, nvar
!$OMP DO
      DO ilat = 1, nlat
        CALL mkl_dft(fftz(1,1,ilat,ivar),inc,nlp2,nlon,nlev,isign)
      ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL

#elif ESSL_DFT

    CALL finish('fftd','OpenMP parallel ESSL DFT not supported yet.')

#else

!$OMP PARALLEL PRIVATE(ilat)
    DO ivar = 1, nvar
!$OMP DO
      DO ilat = 1, nlat
        CALL fft992(fftz(1,1,ilat,ivar),inc,nlp2,nlon,nlev,isign)
      ENDDO
!$OMP END DO
    ENDDO
!$OMP END PARALLEL

#endif

#endif

#endif

  ENDIF

END SUBROUTINE fftd
