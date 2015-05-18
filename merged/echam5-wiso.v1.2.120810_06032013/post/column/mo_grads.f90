module mo_grads
!==============================================================================
!+ $Id: mo_grads.f90,v 1.5 2000/07/31 08:27:03 m214030 Exp $
!
! read and write grads data sets
!==============================================================================
  !-------------
  ! used modules
  !-------------
  use mo_1d_io,           only: sp, wp             ! single,working precision
  use unit_number_module, only: get_unit_number, & ! obtain Fortran unit number
                                return_unit_number ! release Fortran unit number
  use mo_system               ! abort              ! abort on errror condition
  implicit none
!==============================================================================
  !----------------
  ! public entities
  !----------------
  private
  public :: t_ctl, &   ! data type holding .ctl file information
            t_var      ! component of t_ctl holding information on datasets
  public :: read_ctl   ! read .ctl file information
  public :: write_ctl  ! write .ctl file
  public :: init_ctl   ! preset .ctl file
  public :: add_var    ! add variable in .ctl file
  public :: write_var  ! write variable to grads file
  public :: read_var   ! read variable from grads file
  public :: c_month    ! character string encoding for months
  public :: mmm2mm     ! month character string encoding to integer
  public :: destruct   ! destruct t_ctl data type
!==============================================================================
  !-----------
  ! data types
  !-----------
  type t_var
    character(len=32)         :: name          = ''
    integer                   :: levels        = 1
    character(len=16)         :: units         = '99'
    character(len=32)         :: comment       = 'no comment'
    integer                   :: record        = 0
  end type t_var

  type t_ctl
    character(len=128)        :: path          = ''
    character(len=128)        :: dset          = ''
    logical                   :: yrev          = .false.
    logical                   :: sequential    = .false.
    logical                   :: big_endian    = .false.
    logical                   :: little_endian = .false.
    logical                   :: byteswapped   = .false.
    integer                   :: fileheader    = 0
    integer                   :: theader       = 0
    integer                   :: xyheader      = 0
    real(sp)                  :: undef         = -huge(1._sp)
    integer                   :: xdefn         = 1
    integer                   :: ydefn         = 1
    integer                   :: zdefn         = 1
    integer                   :: tdefn         = 1
    character(len=8)          :: xdef          = 'linear'
    character(len=8)          :: ydef          = 'linear'
    character(len=8)          :: tdef          = 'linear'
    character(len=8)          :: zdef          = 'linear'
    real(sp)                  :: xdefi(2)      = (/0.,1./)
    real(sp)                  :: ydefi(2)      = (/0.,1./)
    real(sp)                  :: zdefi(2)      = (/0.,1./)
    character(len=16)         :: tdefi(2)      = (/'0z1jan0000','1hr       '/)
    integer                   :: vars          = 0
    integer                   :: records       = 0
    type (t_var) ,pointer     :: var(:)        => NULL()
    real(wp)     ,pointer     :: ylevels(:)    => NULL()
    real(wp)     ,pointer     :: zlevels(:)    => NULL()
  end type t_ctl
  !----------
  ! constants
  !----------
  character(len=3) ,parameter :: c_month (12) = &
                                 (/'jan','feb','mar','apr','may','jun',&
                                   'jul','aug','sep','oct','nov','dez'/)
  character(len=3) ,parameter :: u_month (12) = &
                                 (/'JAN','FEB','MAR','APR','MAY','JUN',&
                                   'JUL','AUG','SEP','OCT','NOV','DEZ'/)

  integer*4        ,parameter :: i4 = ((49*256+50)*256+51)*256+52
  character(len=*) ,parameter :: cb = '1234'
  character(len=*) ,parameter :: cl = '4321'
  type (t_ctl)     ,save      :: default_ctl
  !-----------
  ! interfaces
  !-----------
  interface write_var
    module procedure write_var2
    module procedure write_var3
  end interface write_var

  interface flip
    module procedure flip2
    module procedure flip3
  end interface flip

  interface init_ctl
    module procedure init_ctl
  end interface init_ctl

  interface destruct
    module procedure destruct_ctl
  end interface destruct
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  function mmm2mm (c) result (i)
  character (len=3) ,intent(in) :: c
  integer                       :: i
  !-------------------------------------------
  ! month character string encoding to integer
  !-------------------------------------------
    integer :: k
    i = 0
    do k=1,12
      if(c==c_month(k) .or. c==u_month(k)) then
        i = k
        exit
      endif
    end do
  end function mmm2mm    
!------------------------------------------------------------------------------
  subroutine init_ctl(ctl, file, ni, nj, ngl, ke, di, dj, lo1, la1, &
                      tdefn, tdefi, tdefd, undef, yrev)
  type (t_ctl)         ,intent(inout)        :: ctl
  character(len=*)     ,intent(in) ,optional :: file
  integer              ,intent(in) ,optional :: ni, nj, ngl, ke
  real(wp)             ,intent(in) ,optional :: di, dj, lo1, la1
  integer              ,intent(in) ,optional :: tdefn
  character(len=*)     ,intent(in) ,optional :: tdefi
  character(len=*)     ,intent(in) ,optional :: tdefd
  real(wp)             ,intent(in) ,optional :: undef
  logical              ,intent(in) ,optional :: yrev
    call destruct (ctl)
    ctl% big_endian = .true.
    if (present(file )) ctl% dset     = file
    if (present(tdefn)) ctl% tdefn    = tdefn
    if (present(tdefi)) ctl% tdefi(1) = tdefi
    if (present(tdefd)) ctl% tdefi(2) = tdefd
    if (present(undef)) ctl% undef    = undef
    if (present(ni   )) ctl% xdefn    = ni
    if (present(lo1  )) ctl% xdefi(1) = lo1
    if (present(di   )) ctl% xdefi(2) = di
    if (present(nj   )) ctl% ydefn    = nj
    if (present(la1  )) ctl% ydefi(1) = la1
    if (present(dj   )) ctl% ydefi(2) = dj
    if (present(ke   )) ctl% zdefn    = ke
    if (present(yrev )) ctl% yrev     = yrev
    if (present(ngl  )) then
      if(ngl>0) then
        ctl% ydef = 'levels'
        allocate (ctl% ylevels(ngl))
        ctl% ydefn = ngl
!        call gauaw (ga=ctl% ylevels)
        ctl% ylevels = 999.
        ctl% ylevels  = - 90._wp/asin(1._wp) * asin (ctl% ylevels)
        ctl% ydefi(1) = ctl% ylevels(1)
      endif
    endif
!    ctl% little_endian =       little()
!    ctl% big_endian    = .not. little()
  end subroutine init_ctl
!------------------------------------------------------------------------------
  subroutine destruct_ctl (ctl)
  type (t_ctl) ,intent(inout) :: ctl
    if (associated(ctl% var))     deallocate (ctl% var)
    if (associated(ctl% ylevels)) deallocate (ctl% ylevels)
    if (associated(ctl% zlevels)) deallocate (ctl% zlevels)
    ctl = default_ctl
  end subroutine destruct_ctl
!==============================================================================
  subroutine add_var (ctl, name, levels, comment)
  !---------------------------
  ! add a variable to the list 
  !---------------------------
  type (t_ctl)     ,intent(inout)        :: ctl
  character(len=*) ,intent(in)           :: name
  integer          ,intent(in) ,optional :: levels
  character(len=*) ,intent(in) ,optional :: comment
    type (t_var) ,allocatable  :: var(:)
    integer                    :: n, i
    !---------------------------------------------
    ! check if variable name is already registered
    !---------------------------------------------
    do i=1,ctl%vars
      if (ctl% var(i)% name == name) return
    end do
    !--------------
    ! increase list
    !--------------
    n = ctl%vars + 1
    if (associated(ctl% var)) then
      allocate (var(n-1));    var = ctl% var;       deallocate (ctl% var)
      allocate (ctl% var(n)); ctl% var(:n-1) = var; deallocate (var)
    else
      allocate (ctl% var(n))
    endif
    ctl%vars = n
    !----------------------
    ! insert new list entry
    !----------------------
    ctl% var(n)% name    = name
    ctl% var(n)% levels  = ctl% zdefn
    if (present(levels))  ctl% var(n)% levels  = levels
    if (present(comment)) ctl% var(n)% comment = comment
    ctl% var(n)% record = ctl% records + 1
    ctl% records = ctl% records + ctl% var(n)% levels
  end subroutine add_var
!------------------------------------------------------------------------------
  subroutine write_var2 (ctl, x, name, t, comment, yrev, iostat)
  !---------------------------------
  ! write a data set to a grads file
  !---------------------------------
  type (t_ctl)     ,intent(inout)        :: ctl     ! grads ctl structure
  real(wp)         ,intent(in)           :: x (:,:) ! field to write
  character(len=*) ,intent(in)           :: name    ! name of data set
  integer          ,intent(in) ,optional :: t       ! time slice
  character(len=*) ,intent(in) ,optional :: comment ! data set description
  logical          ,intent(in) ,optional :: yrev    ! flip n-s direction
  integer          ,intent(out),optional :: iostat
    integer  :: iunit, irecl, it, irec, i, ny
    logical  :: yf
    real(sp) :: buf (size(x,1),size(x,2))
    integer  :: ibuf(size(x,1),size(x,2))
    it = 1         ;if(present(t))     it = t
    yf = ctl% yrev ;if(present(yrev))  yf = yrev .neqv. ctl% yrev
    ny = size(x,2)
    if(present(iostat)) iostat = 0
    !---------------------
    ! check shape of array
    !---------------------
    if (any (shape(x) /= (/ctl% xdefn, ctl% ydefn/))) then
      if(present(iostat)) then
        iostat = 1
        return
      endif
      print *,'write_var: abort, shape(x) /= (/xdefn,ydefn/)'
      print *,'write_var: abort,',shape(x),'/=',ctl% xdefn, ctl% ydefn
      call abort('write_var: shape(x) /= (/xdefn,ydefn/)')
    endif
    !---------------------------
    ! add entry to ctl structure
    !---------------------------
    if (it==1) call add_var (ctl, name, levels=1, comment=comment)
    ctl% tdefn = max (ctl% tdefn, it)
    !------------------------
    ! determine record number
    !------------------------
    irec = 0
    do i=1,ctl% vars
      if(ctl%var(i)% name == name) then
        irec = ctl%var(i)% record
        exit
      endif
    end do
    if (irec == 0) then
      !---------------------
      ! variable not present
      !---------------------
      if(present(iostat)) then
        iostat = 1
        return
      else
        write(6,*) 'write_var2: variable is not in ctl list.'
        write(6,*) '            variable =',name
        do i=1,ctl% vars
          write(6,*) '           listentry =',ctl%var(i)% name,&
                                              ctl%var(i)% name == name
        end do
      endif
      call abort('write_var2: variable is not in ctl list.')
    endif
    irec = irec + (it-1) * ctl% records
    !---------------
    ! write data set
    !---------------
    if(yf) then
      do i=1,ny
        buf(:,i) = x(:,ny-i+1)
      end do
    else
      buf = x
    endif
    ibuf = reshape (transfer (buf, ibuf), (/size(x,1),size(x,2)/))
    if (byteswapped(ctl)) call flip(ibuf)
    inquire (iolength = irecl) ibuf
    iunit = get_unit_number()
    open(iunit,file=ctl%dset,access='direct',recl=irecl)
    write(iunit,rec=irec) ibuf
    close(iunit)
    call return_unit_number(iunit)
  end subroutine write_var2
!------------------------------------------------------------------------------
  subroutine write_var3 (ctl, x, name, t, comment, yrev, iostat)
  !---------------------------------
  ! write a data set to a grads file
  !---------------------------------
  type (t_ctl)     ,intent(inout)        :: ctl       ! grads ctl structure
  real(wp)         ,intent(in)           :: x (:,:,:) ! field to write
  character(len=*) ,intent(in)           :: name      ! name of data set
  integer          ,intent(in) ,optional :: t         ! time slice
  character(len=*) ,intent(in) ,optional :: comment   ! data set description
  logical          ,intent(in) ,optional :: yrev      ! flip n-s direction
  integer          ,intent(out),optional :: iostat
    integer  :: iunit, irecl, it, irec, i, ny
    logical  :: yf
    real(sp) :: buf (size(x,1),size(x,2),size(x,3))
    integer  :: ibuf(size(x,1),size(x,2),size(x,3))
    it = 1         ;if(present(t))     it = t
    yf = ctl% yrev ;if(present(yrev))  yf = yrev .neqv. ctl% yrev
    ny = size(x,2)
    if(present(iostat)) iostat = 0
    !---------------------
    ! check shape of array
    !---------------------
    if (any (shape(x(:,:,1)) /= (/ctl% xdefn, ctl% ydefn/))) then
      if(present(iostat)) then
        iostat = 1
        return
      endif
      call abort('write_var: shape(x) /= (/xdefn,ydefn/)')
    endif
    !---------------------------
    ! add entry to ctl structure
    !---------------------------
    if (it==1) call add_var (ctl, name, levels=size(x,3), comment=comment)
    ctl% tdefn = max (ctl% tdefn, it)
    !------------------------
    ! determine record number
    !------------------------
    irec = 0
    do i=1,ctl% vars
      if(ctl%var(i)% name == name) then
        irec = ctl%var(i)% record
        exit
      endif
    end do
    if (irec == 0) then
      !---------------------
      ! variable not present
      !---------------------
      if(present(iostat)) then
        iostat = 1
        return
      else
        write(6,*) 'write_var3: variable is not in ctl list.'
        write(6,*) '            variable =',name
        do i=1,ctl% vars
          write(6,*) '           listentry =',ctl%var(i)% name,&
                                              ctl%var(i)% name == name
        end do
      endif
      call abort('write_var3: variable is not in ctl list.')
    endif
    irec = irec + (it-1) * ctl% records
    !---------------
    ! write data set
    !---------------
    if(yf) then
      do i=1,ny
        buf(:,i,:) = x(:,ny-i+1,:)
      end do
    else
      buf = x
    endif
    ibuf = reshape (transfer (buf, ibuf), (/size(x,1),size(x,2),size(x,3)/))
    if (byteswapped(ctl)) call flip(ibuf)
    inquire (iolength = irecl) ibuf(:,:,1)
    iunit = get_unit_number()
    open(iunit,file=ctl%dset,access='direct',recl=irecl)
    do i=1,size(x,3)
      write(iunit,rec=irec) ibuf(:,:,i)
      irec=irec+1
    end do
    close(iunit)
    call return_unit_number(iunit)
  end subroutine write_var3
!------------------------------------------------------------------------------
  subroutine read_var (ctl, x, name, t, yrev, iostat)
  !---------------------------------
  ! write a data set to a grads file
  !---------------------------------
  type (t_ctl)     ,intent(in)           :: ctl     ! grads ctl structure
  real(wp)         ,intent(out)          :: x (:,:) ! field to write
  character(len=*) ,intent(in)           :: name    ! name of data set
  integer          ,intent(in) ,optional :: t       ! time slice
  logical          ,intent(in) ,optional :: yrev    ! flip n-s direction
  integer          ,intent(out),optional :: iostat  ! io stat return variable

    integer                  :: iunit, irecl, it, irec, i, ny
    real(sp)                 :: buf (size(x,1),size(x,2))
    integer                  :: ibuf(size(x,1),size(x,2))
    logical                  :: yf
    it = 1; if(present(t)) it = t
    yf = ctl% yrev ;if(present(yrev))  yf = yrev .neqv. ctl% yrev
    if(present(iostat)) iostat = 0
    !---------------------
    ! check shape of array
    !---------------------
    if (any (shape(x) /= (/ctl% xdefn, ctl% ydefn/))) then
      write(6,*) 'read_var: shape(x) /= (/xdefn,ydefn/)'
      write(6,*) '          shape(x)        =',shape(x)
      write(6,*) '          (/xdefn,ydefn/) =',ctl% xdefn, ctl% ydefn
      call abort('read_var: shape(x) /= (/xdefn,ydefn/)')
    endif
    !------------------------
    ! determine record number
    !------------------------
    irec = 0
    do i=1,ctl% vars
      if(ctl%var(i)% name == name) then
        irec = ctl%var(i)% record
        exit
      endif
    end do
    if (irec == 0) then
      !---------------------
      ! variable not present
      !---------------------
      if(present(iostat)) then
        iostat = 1
        return
      else
        write(6,*) 'read_var: variable is not in ctl list.'
        write(6,*) '          variable =',name
        do i=1,ctl% vars
          write(6,*) '         listentry =',ctl%var(i)% name,&
                                            ctl%var(i)% name == name
        end do
      endif
      call abort('read_var: variable is not in ctl list.')
    endif
    irec = irec + (it-1) * ctl% records
    !--------------
    ! read data set
    !--------------
    inquire (iolength = irecl) ibuf
    iunit = get_unit_number()
    if(ctl%dset(1:1) /= '^') then
      open(iunit,file=ctl%dset,access='direct',recl=irecl,status='old')
    else
      open(iunit,file=trim(ctl%path)//ctl%dset(2:),&
                               access='direct',recl=irecl,status='old')
    endif
    read(iunit,rec=irec) ibuf
    close(iunit)
    call return_unit_number(iunit)
    if (byteswapped(ctl)) call flip(ibuf)
    buf = reshape (transfer (ibuf, buf), (/size(x,1),size(x,2)/))
    if(yf) then
      ny = size(x,2)
      do i=1,ny
        x(:,ny-i+1) = buf(:,i)
      end do
    else
      x = buf
    endif
  end subroutine read_var
!------------------------------------------------------------------------------
  subroutine write_ctl (ctl, unit)
  type (t_ctl) ,intent(in)           :: ctl
  integer      ,intent(in) ,optional :: unit
    integer            :: i, ios, iu
    character(len=128) :: file
    !----------------------
    ! open file if required
    !----------------------
    if (present(unit)) then
      iu = unit
    else
      iu = get_unit_number()
      file = trim(ctl%dset)//'.ctl'
      open (iu,file=file,action='write',iostat=ios)
      if(ios/=0) call abort('write_ctl: cannot open '//trim(file))
    endif
    !----------
    ! write ctl
    !----------
    if(iu==6) then
      write(iu,'(a,1x,a)')            'path',trim(ctl% path)
    endif
    write(iu,'(a,1x,a)')              'dset',trim(ctl% dset)

    if(ctl% sequential .or. ctl% yrev .or. ctl% little_endian .or. &
       ctl% big_endian .or. ctl% byteswapped) then
      write                       (iu,'(a)',advance='no') 'options'
      if(ctl% sequential)    write(iu,'(a)',advance='no') ' sequential'
      if(ctl% yrev)          write(iu,'(a)',advance='no') ' yrev'
      if(ctl% little_endian) write(iu,'(a)',advance='no') ' little_endian'
      if(ctl% big_endian)    write(iu,'(a)',advance='no') ' big_endian'
      if(ctl% byteswapped)   write(iu,'(a)',advance='no') ' byteswapped'
      write                    (iu,'()')
    endif

    if(ctl% fileheader/=0) write(iu,'(a,1x,i8)') 'fileheader', ctl% fileheader
    if(ctl% theader   /=0) write(iu,'(a,1x,i8)') 'theader   ', ctl% theader
    if(ctl% xyheader  /=0) write(iu,'(a,1x,i8)') 'xyheader  ', ctl% xyheader

    write(iu,'(a,1x,g12.6)')          'undef',ctl% undef
    write(iu,'(a,1x,i6,1x,a,2g12.6)') 'xdef' ,ctl%xdefn,ctl%xdef,&
                                             ctl%xdefi(1),ctl%xdefi(2)
    select case (ctl%ydef)
    case ('linear','LINEAR')
      write(iu,'(a,1x,i6,1x,a,2g12.6)') 'ydef' ,ctl%ydefn,ctl%ydef,&
                                               ctl%ydefi(1),ctl%ydefi(2)
    case ('levels','LEVELS')
      write(iu,'(a,1x,i6,1x,a,2g12.6)') 'ydef' ,ctl%ydefn,ctl%ydef
      write(iu,'(15x,6f8.3)') real(ctl%ylevels,sp)
    end select
    select case (ctl%zdef)
    case ('linear','LINEAR')
      write(iu,'(a,1x,i6,1x,a,2g12.6)') 'zdef' ,ctl%zdefn,ctl%zdef,&
                                               ctl%zdefi(1),ctl%zdefi(2)
    case ('levels','LEVELS')
      write(iu,'(a,1x,i6,1x,a,2g12.6)') 'zdef' ,ctl%zdefn,ctl%zdef
      write(iu,'(15x,6f9.2)') real(ctl%zlevels,sp)
    end select
    write(iu,'(a,1x,i6,3(1x,a))')     'tdef' ,ctl%tdefn,ctl%tdef,&
                                             ctl%tdefi(1),ctl%tdefi(2)
    write(iu,'(a,1x,i4)')             'vars' ,ctl%vars
    do i=1,ctl%vars
      write(iu,"(a,1x,i4,1x,a,1x,a)") &
        ctl%var(i)%name, ctl%var(i)%levels,         &
        trim(ctl%var(i)%units),trim(ctl%var(i)%comment)
    end do
    write(iu,'(a)')                   'endvars'
    !----------------------
    ! close file if required
    !----------------------
    if (.not.present(unit)) then
      close (iu)
      call return_unit_number(iu)
    endif
  end subroutine write_ctl
!------------------------------------------------------------------------------
  subroutine read_ctl (ctl, file, iostat)
  type (t_ctl)     ,intent(inout)         :: ctl
  character(len=*) ,intent(in)            :: file
  integer          ,intent(out) ,optional :: iostat
    type(t_ctl)        :: empty
    integer            :: iu, ios, idx, idx2, i
    character(len=128) :: line
    character(len=16)  :: code
    iu = get_unit_number()
    open (iu,file=file,action='read',status='old',iostat=ios)
    if(present(iostat)) iostat = ios
    if(ios/=0 .and. .not.present(iostat)) &
      call abort('read_ctl: cannot open '//trim(file))
    if(associated(ctl% var    )) deallocate (ctl% var)
    if(associated(ctl% ylevels)) deallocate (ctl% ylevels)
    ctl = empty
    do
      read(iu,'(a)',iostat=ios) line
      if(ios/=0) exit
      idx=index(line,' ')
      if(idx>0) then
        code=line(:idx-1)
        select case (code)
        case ('dset','DSET')
          ctl% dset = line(idx+1:)
        case ('undef','UNDEF')
          read(line(idx+1:),*) ctl% undef
        case ('options','OPTIONS')
          ctl% sequential    = index (line,'sequential')    /= 0
          ctl% yrev          = index (line,'yrev')          /= 0
          ctl% little_endian = index (line,'little_endian') /= 0
          ctl% big_endian    = index (line,'big_endian')    /= 0
          ctl% byteswapped   = index (line,'byteswapped')   /= 0
        case ('xdef','XDEF')
          read(line(idx+1:),*) ctl%xdefn,ctl%xdef,ctl%xdefi(1),ctl%xdefi(2)
        case ('ydef','YDEF')
          read(line(idx+1:),*) ctl%ydefn,ctl%ydef
          select case (ctl%ydef)
          case ('linear','LINEAR')
            read(line(idx+1:),*) ctl%ydefn,ctl%ydef,ctl%ydefi(1),ctl%ydefi(2)
          case ('levels','LEVELS')
            allocate (ctl% ylevels (ctl%ydefn))
            read(iu,*) ctl% ylevels
          end select
        case ('zdef','ZDEF')
          read(line(idx+1:),*) ctl%zdefn,ctl%zdef
          select case (ctl%zdef)
          case ('linear','LINEAR')
            read(line(idx+1:),*) ctl%zdefn,ctl%zdef,ctl%zdefi(1),ctl%zdefi(2)
          case ('levels','LEVELS')
            allocate (ctl% zlevels (ctl%zdefn))
            read(iu,*) ctl% zlevels
          end select
        case ('tdef','TDEF')
          read(line(idx+1:),*) ctl%tdefn,ctl%tdef,ctl%tdefi(1),ctl%tdefi(2)
        case ('vars','VARS')
          read(line(idx+1:),*) ctl%vars
          allocate (ctl%var(ctl%vars))
          ctl% records = 0
          do i=1,ctl%vars
            read(iu,'(a)',iostat=ios) line
            if(ios/=0) exit
            ctl%var(i)%units     = '99'
            ctl% var(i)% comment = ''
            read(line,*,iostat=ios) ctl%var(i)%name, ctl%var(i)%levels
            ctl%var(i)%levels = max(1,ctl%var(i)%levels)
            idx =       verify(line         ,' ')  ! start of name
            idx = idx + index (line(idx +1:),' ')  ! end of name
            idx = idx + verify(line(idx +1:),' ')  ! start of levels
            idx = idx + index (line(idx +1:),' ')  ! end of levels
            idx = idx + verify(line(idx +1:),' ')  ! start of units 
            idx2= idx + index (line(idx +1:),' ')  ! end of units
            if(idx2>idx) ctl%var(i)%units = line(idx:idx2-1)
            idx = idx2+ verify(line(idx2+1:),' ')  ! start of comment
            ctl% var(i)% comment = line(idx:)
            ctl% var(i)% record = ctl% records + 1
            ctl% records = ctl% records + ctl% var(i)% levels
          end do
        case ('endvars','ENDVARS')
          exit
        case default
          print *,'read_ctl: UNKNOWN LINE:',trim(line)
        end select
      endif
    end do
    close (iu)
    call return_unit_number(iu)
  end subroutine read_ctl
!------------------------------------------------------------------------------
  function little()
    logical :: little
    little = (cl==transfer(i4,cl))
  end function little
!------------------------------------------------------------------------------
  function byteswapped (ctl)
  type (t_ctl) ,intent(in) :: ctl
  logical                  :: byteswapped
    byteswapped = ctl% byteswapped
    if (ctl% big_endian)    byteswapped = little()
    if (ctl% little_endian) byteswapped = .not. little()
  end function byteswapped
!------------------------------------------------------------------------------
  subroutine flip2 (x)
  integer ,intent(inout) :: x(:,:)
    character(len=4) :: c
    character(len=4) :: t
    integer          :: i,j
    do j=1,size(x,2)
      do i=1,size(x,1)
        c = transfer (x(i,j),c)
        t(1:1) = c(4:4)
        t(2:2) = c(3:3)
        t(3:3) = c(2:2)
        t(4:4) = c(1:1)
        x(i,j) = transfer (t,x(i,j))
      end do
    end do
  end subroutine flip2
!------------------------------------------------------------------------------
  subroutine flip3 (x)
  integer ,intent(inout) :: x(:,:,:)
    character(len=4) :: c
    character(len=4) :: t
    integer          :: i,j,k
    do k=1,size(x,3)
      do j=1,size(x,2)
        do i=1,size(x,1)
          c = transfer (x(i,j,k),c)
          t(1:1) = c(4:4)
          t(2:2) = c(3:3)
          t(3:3) = c(2:2)
          t(4:4) = c(1:1)
          x(i,j,k) = transfer (t,x(i,j,k))
        end do
      end do
    end do
  end subroutine flip3
!------------------------------------------------------------------------------
end module mo_grads
