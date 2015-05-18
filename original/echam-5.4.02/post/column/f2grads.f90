program f2grads
  !-------------
  ! modules used
  !-------------
  use mo_system                   ! getarg, iarg
  use mo_1d_io, only: wp,&        ! working precision
                      ln,        &! length of identifier string
                      readfile    ! read forcing file
  use mo_grads, only: t_ctl,     &! grads data type definition
                      init_ctl,  &! init grads data type
                      add_var,   &! enter variable to ctl info
                      write_var, &! write var. to grads
                      write_ctl   ! write grads .ctl file
  implicit none
  !----------
  ! variables
  !----------
  character(len=128)    :: ifile, ofile        ! file names
  type (t_ctl)          :: ctl                 ! grads data type
  logical               :: exist_i             ! input files exist
  integer               :: ios                 ! iostat return value
  character(len=ln)     :: name                ! name of variable read
  integer               :: rank, shap(3), index(2), lonlat(2)
  real(wp) ,pointer     :: xi(:)               ! fieldt input
  real(wp) ,pointer     :: xo(:,:,:)           ! fieldt output
  integer               :: narg                ! number of arguments provided
  character             :: copt                ! option (-p or -z)
  integer               :: nstep               ! time step counter
  integer               :: i                   ! do loop index
  integer               :: nlev                ! number of levels
  integer               :: nsum                ! # of summands for level calc. 
  real(wp)              :: ps                  ! surface pressure
  real(wp) ,allocatable :: ak(:)               ! Vertical Coordinate Table
  real(wp) ,allocatable :: bk(:)               ! Vertical Coordinate Table
  character (len=16)    :: iname
  integer               :: ierr                ! error return argument
  !--------------
  ! get arguments
  !--------------
  ifile = ''
  copt  = ''
  iname = ''
  narg = iargc()
  do i = 1,narg
    call getarg (i,ifile)
    if(iname == '################') then
      iname = ifile
    else
      select case (ifile)
      case ('-n')
        iname = '################'
      case ('-p')
        copt  = 'p'
        ifile = ''
      case ('-z')
        copt  = 'z'
        ifile = ''
      end select
    endif
  end do
  !---------------------------
  ! mark arguments not present
  !---------------------------
  if(ifile=='') ifile = '?'
  ofile=trim(ifile)//'.grads'
  !-----------------------
  ! do input files exist ?
  !-----------------------
  inquire (file=ifile,exist=exist_i)
  !----------------
  ! print arguments
  !----------------
  write(0,*)
  write(0,*)'convert forcing file to grib file'
  write(0,*)'  from : ',trim(ifile)
  write(0,*)'  to   : ',trim(ofile)
  write(0,*)
  !-----------------------------
  ! print error message and exit
  !-----------------------------
  if (.not.exist_i) write(0,*)'         ',trim(ifile),' DOES NOT EXIST'
  if (.not.exist_i .or. ifile=='?' ) then
    write(0,*)
    write(0,*)'USAGE     : f2grads [-p|-z] ifile'
    write(0,*)
    write(0,*)'       -p : use approximate pressure coordinates [hPa]'
    write(0,*)'       -z : use approximate geopotential height coordinates [m]'
    write(0,*)'  default : use model levels'
    write(0,*)
    stop
  endif
  !-----------
  ! open files
  !-----------
  open(1,file=ifile,status='old',form='unformatted',iostat=ios)
                                              if (ios/=0) goto 93
  open(2,file=ofile             ,form='unformatted',iostat=ios)
                                              if (ios/=0) goto 92
  !-------------------
  ! loop: read - write
  !-------------------
  call init_ctl(ctl,ofile)
  nullify (xi)
  index  = -1
  lonlat = (/0,0/)
  nsum   = 0
  do i=1,2
    nstep  = 0
    nlev   = 0
    do
      !------------
      ! read column
      !------------
      call readfile (xi, name, 1, index, lonlat, rank=rank, shap=shap)
      if(.not.associated(xi)) goto 91
      select case (name)
      case ('AK')
        if(nlev==0) call init_levels (shap(2)-1)
        if(.not.allocated(ak)) then
          allocate (ak(size(xi)))
          ak = xi
        endif
      case ('BK')
        if(nlev==0) call init_levels (shap(2)-1)
        if(.not.allocated(bk)) then
          allocate (bk(size(xi)))
          bk = xi
        endif
      case ('DTIME')
        write(ctl% tdefi(2),'(i5,"mn")') nint(xi/60._wp)
      case ('NSTEP')
        nstep = nstep + 1
        if(nstep==1) then
        endif
      case default 
      !------
      ! write
      !------
      if(nstep>0.and.shap(1)==1.and.shap(3)<=1) then
        allocate (xo(1,1,shap(2)))
        xo(1,1,:) = xi
        ierr=0
        if(iname==''.or.iname==name) then
          if(i==1)call add_var  (ctl,   name,size(xi),comment=name)
          if(i==2)call write_var(ctl,xo,name,t=nstep ,comment=name,iostat=ierr)
        endif
        if(ierr/=0) then
          write (6,*)'WARNING: nstep=',nstep,', cannot write: ',trim(name)
        endif
        deallocate (xo)
        !
        ! process levels
        !
        if(nlev==0) call init_levels (shap(2)-1)
        if (i==1) then
          if(nlev>1 .and. copt=='p' .and. name=='PSREF') then
            ps   = ps + xi(1)
            nsum = nsum + 1
          endif
          if(nlev>1 .and. copt=='z'.and.(name=='geom1'.or.name=='GEOM1')) then
            ctl% zlevels = ctl% zlevels + xi
            nsum = nsum + 1
          endif
        endif
      endif
      end select
    end do
91  continue
    !-------------
    ! EOF on input
    !-------------
    print *
    write(0,*)'EOF reading ',trim(ifile)
    rewind (1)
  end do
  close (1)
  !----------------------------
  ! set levels, write .ctl-file
  !----------------------------
  ctl% zdefi(1) = 1
  if(nlev > 1) then
    if (nsum>0) then
      select case (copt)
      case ('p')
        if(allocated(ak) .and. allocated(bk)) then
          ps = exp(ps/nsum)
          ctl% zlevels          = (ak(1:nlev) + ps*bk(1:nlev)) / 200._wp
          ctl% zlevels(1:nlev-1)= ctl% zlevels(1:nlev-1) + ctl% zlevels(2:nlev)
          ctl% zlevels(nlev)    = ps / 100._wp
        endif
      case ('z')
        ctl% zlevels = ctl% zlevels / nsum / 9.81
      case default
        ctl% zlevels = (/(i,i=1,nlev)/)
      end select
    else
      ctl% zlevels = (/(i,i=1,nlev)/)
    endif
    ctl% zdefi(1) = ctl% zlevels(1)  
  endif
  call write_ctl (ctl)
  write(0,*) 'f2grads      :'
  write(0,*) '  levels     :',nlev
  write(0,*) '  variables  :',ctl% vars
  write(0,*) '  time slices:',ctl% tdefn
  stop
92 continue
  !----------------
  ! error condition
  !----------------
  write(0,*)
  write(0,*)'ERROR writing ',trim(ofile)
  stop
93 continue
  write(0,*)
  write(0,*)'ERROR opening ',trim(ifile)
  stop
contains
  subroutine init_levels (n)
  !
  ! set levels
  !
  integer :: n ! number of levels
    if(n > 1 .and. nlev == 0) then
      nlev = n
      ctl% zdef = 'levels'
      allocate (ctl% zlevels(n))
      ctl% zlevels = 0.
      ctl% zdefn   = nlev
    endif
  end subroutine init_levels
end program f2grads


