module decompInitMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use spmdMod         , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils      , only : endrun
  use elm_varctl      , only : iulog, use_fates
  use FatesInterfaceTypesMod, only : fates_maxElementsPerSite
  use decompMod
  use topounit_varcon   , only : max_topounits, has_topounit
  use domainMod         , only: ldomain
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public decompInit_lnd          ! initializes lnd grid decomposition into clumps and processors
  public decompInit_clumps       ! initializes atm grid decomposition into clumps
  public decompInit_gtlcp         ! initializes g,l,c,p decomp info
  public decompInit_lnd_using_gp ! initialize lnd grid decomposition into clumps and processors using graph partitioning approach
  public decompInit_ghosts       ! initialize ghost/halo for land grid
  !
  ! !PRIVATE TYPES:
  private
  integer, pointer :: lcid(:)       ! temporary for setting ldecomp
  integer            :: nglob_x, nglob_y ! global sizes
  integer, parameter :: dbug=0           ! 0 = min, 1=normal, 2=much, 3=max
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine decompInit_lnd(lni,lnj,amask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use elm_varctl, only : nsegspc
    use decompMod  , only : gsMap_lnd_gdc2glo, nclumps, clumps
    use decompMod  , only : bounds_type, get_proc_bounds, procinfo
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: amask(:)
    integer , intent(in) :: lni,lnj   ! domain global size
    !
    ! !LOCAL VARIABLES:
    integer :: lns                    ! global domain size
    integer :: ln,lj                  ! indices
    integer :: ag,an,ai,aj            ! indices
    integer :: numg                   ! number of land gridcells
    logical :: seglen1                ! is segment length one
    real(r8):: seglen                 ! average segment length
    real(r8):: rcid                   ! real value of cid
    integer :: cid,pid                ! indices
    integer :: n,m,ng                 ! indices
    integer :: ier                    ! error code
    integer :: beg,end,lsize,gsize    ! used for gsmap init
    integer, pointer :: gindex(:)     ! global index for gsmap init
    integer, pointer :: clumpcnt(:)   ! clump index counter
    type(bounds_type) :: bounds       ! contains subgrid bounds data
    integer, allocatable :: proc_ncell(:) ! number of cells assigned to a process
    integer, allocatable :: proc_begg(:)  ! beginning cell index assigned to a process
    !------------------------------------------------------------------------------

    lns = lni * lnj

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write(iulog,*) 'decompInit_lnd(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    else
       write(iulog,*)'clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! allocate and initialize procinfo (from decompMod.F90) and clumps 
    ! beg and end indices initialized for simple addition of cells later 

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for procinfo%cid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    procinfo%nclumps = clump_pproc
    procinfo%cid(:)  = -1
    procinfo%ncells  = 0
    procinfo%ntunits  = 0
    procinfo%nlunits = 0
    procinfo%ncols   = 0
    procinfo%npfts   = 0
    procinfo%nCohorts = 0
    procinfo%begg    = 1
    procinfo%begt    = 1
    procinfo%begl    = 1
    procinfo%begc    = 1
    procinfo%begp    = 1
    procinfo%begCohort    = 1
    procinfo%endg    = 0
    procinfo%endt    = 0
    procinfo%endl    = 0
    procinfo%endc    = 0
    procinfo%endp    = 0
    procinfo%endCohort    = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for clumps'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    clumps(:)%owner   = -1
    clumps(:)%ncells  = 0
    clumps(:)%ntunits = 0
    clumps(:)%nlunits = 0
    clumps(:)%ncols   = 0
    clumps(:)%npfts   = 0
    clumps(:)%nCohorts = 0
    clumps(:)%begg    = 1
    clumps(:)%begt    = 1
    clumps(:)%begl    = 1
    clumps(:)%begc    = 1
    clumps(:)%begp    = 1
    clumps(:)%begCohort    = 1
    clumps(:)%endg    = 0
    clumps(:)%endt    = 0
    clumps(:)%endl    = 0
    clumps(:)%endc    = 0
    clumps(:)%endp    = 0
    clumps(:)%endCohort    = 0

    ! assign clumps to proc round robin 
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write(iulog,*) 'decompInit_lnd(): round robin pid error ',n,pid,npes
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    ! count total land gridcells
    numg = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          numg = numg + 1
       endif
    enddo

    if (npes > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of processes exceeds number ', &
            'of land grid cells',npes,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    if (nclumps > numg) then
       write(iulog,*) 'decompInit_lnd(): Number of clumps exceeds number ', &
            'of land grid cells',nclumps,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (float(numg)/float(nclumps) < float(nsegspc)) then
       seglen1 = .true.
       seglen = 1.0_r8
    else
       seglen1 = .false.
       seglen = dble(numg)/(dble(nsegspc)*dble(nclumps))
    endif

    if (masterproc) then
       write(iulog,*) ' decomp precompute numg,nclumps,seglen1,avg_seglen,nsegspc=', &
            numg,nclumps,seglen1,&
            sngl(seglen),sngl(dble(numg)/(seglen*dble(nclumps)))
    end if

    ! Assign gridcells to clumps (and thus pes) ---

    allocate(lcid(lns), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for lcid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    lcid(:) = 0
    ng = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          ng = ng  + 1

          !--- give to clumps in order based on nsegspc
          if (seglen1) then
             cid = mod(ng-1,nclumps) + 1
          else
             rcid = (dble(ng-1)/dble(numg))*dble(nsegspc)*dble(nclumps)
             cid = mod(int(rcid),nclumps) + 1
          endif
          lcid(ln) = cid

          !--- give gridcell cell to pe that owns cid ---
          !--- this needs to be done to subsequently use function
          !--- get_proc_bounds(begg,endg)
          if (iam == clumps(cid)%owner) then
             procinfo%ncells  = procinfo%ncells  + 1
          endif
          if (iam >  clumps(cid)%owner) then
             procinfo%begg = procinfo%begg + 1
          endif
          if (iam >= clumps(cid)%owner) then
             procinfo%endg = procinfo%endg + 1
          endif

          !--- give gridcell to cid ---
          !--- increment the beg and end indices ---
          clumps(cid)%ncells  = clumps(cid)%ncells  + 1
          do m = 1,nclumps
             if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                 (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
                clumps(m)%begg = clumps(m)%begg + 1
             endif

             if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
                 (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
                clumps(m)%endg = clumps(m)%endg + 1
             endif
          enddo

       end if
    enddo

    ! Set gsMap_lnd_gdc2glo 
    allocate(ldecomp%gdc2glo(numg), stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error for proc_ncell'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    ldecomp%gdc2glo(:) = 0
    allocate(clumpcnt(nclumps),stat=ier)
    if (ier /= 0) then
       write(iulog,*) 'decompInit_lnd(): allocation error1 for clumpcnt'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! clumpcnt is the start gdc index of each clump

    ag = 0
    clumpcnt = 0
    ag = 1
    do pid = 0,npes-1
    do cid = 1,nclumps
       if (clumps(cid)%owner == pid) then
         clumpcnt(cid) = ag
         ag = ag + clumps(cid)%ncells
       endif
    enddo
    enddo

    ! now go through gridcells one at a time and increment clumpcnt
    ! in order to set gdc2glo

    do aj = 1,lnj
    do ai = 1,lni
       an = (aj-1)*lni + ai
       cid = lcid(an)
       if (cid > 0) then
          ag = clumpcnt(cid)
          ldecomp%gdc2glo(ag) = an
          clumpcnt(cid) = clumpcnt(cid) + 1
       end if
    end do
    end do

    ! Initialize global gindex (non-compressed, includes ocean points)
    ! Note that gsMap_lnd_gdc2glo goes from (1:endg)
    nglob_x = lni !  decompMod module variables
    nglob_y = lnj !  decompMod module variables
    call get_proc_bounds(bounds)
    allocate(gsMap_lnd_gdc2glo(1:bounds%endg))
    do n = procinfo%begg,procinfo%endg
       gsMap_lnd_gdc2glo(n-procinfo%begg+1) = ldecomp%gdc2glo(n)
    enddo

    deallocate(clumpcnt)

!hh!    ! Set gsMap_lnd_gdc2glo (the global index here includes mask=0 or ocean points)
!hh!
!hh!    call get_proc_bounds(beg, end)
!hh!    allocate(gindex(beg:end))
!hh!    do n = beg,end
!hh!       gindex(n) = ldecomp%gdc2glo(n)
!hh!    enddo
!hh!    lsize = end-beg+1
!hh!    gsize = lni * lnj
!hh!    call mct_gsMap_init(gsMap_lnd_gdc2glo, gindex, mpicom, comp_id, lsize, gsize)
!hh!    deallocate(gindex)

    ! Diagnostic output
    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points               = ',lni
       write(iulog,*)'   latitude points                = ',lnj
       write(iulog,*)'   total number of land gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process             = ',clump_pproc
       write(iulog,*)
    end if
    call shr_sys_flush(iulog)

  end subroutine decompInit_lnd

  !------------------------------------------------------------------------------
!hh!  subroutine decompInit_clumps(glcmask)
  subroutine decompInit_clumps(lni,lnj,glcmask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure.  This assumes each pe has the same number of clumps
    ! set by clump_pproc
    !
    ! !USES:
    use subgridMod, only : subgrid_get_gcellinfo
    use decompMod      , only : bounds_type, clumps, nclumps, procinfo
    use decompMod      , only : get_proc_global, get_proc_bounds
    use decompMod      , only : numg, numl, numc, nump, numCohort
    use decompMod      , only : gsMap_lnd_gdc2glo
    use spmdMod        , only : MPI_INTEGER, MPI_SUM
    !
    ! !ARGUMENTS:
    implicit none
    integer                 , intent(in) :: lni,lnj ! land domain global size
    integer , pointer, optional   :: glcmask(:)  ! glc mask
    !integer , pointer, optional   :: num_tunits_per_grd(:)  ! Number of topounits per grid
    !
    ! !LOCAL VARIABLES:
    integer           :: ln,an             ! indices
    integer           :: i,g,l,k           ! indices
    integer           :: cid,pid           ! indices
    integer           :: n,m,np            ! indices
    integer           :: anumg             ! lnd num gridcells
    integer           :: icells            ! temporary
    integer           :: begg, endg        ! temporary
    integer           :: itunits            ! temporary
    integer           :: ilunits            ! temporary
    integer           :: icols              ! temporary
    integer           :: ipfts              ! temporary
    integer           :: icohorts           ! temporary
    integer           :: ier                ! error code
    integer           :: npmin,npmax,npint ! do loop values for printing
    integer           :: clmin,clmax       ! do loop values for printing
    type(bounds_type) :: bounds            ! bounds
    integer :: glev, tlev, llev, clev, plev, hlev  ! order of subgrid levels in the allvec arrays
    integer :: nlev               ! number of subgrid levels
    integer, allocatable :: allvecg(:,:)  ! temporary vector "global"
    integer, allocatable :: allvecl(:,:)  ! temporary vector "local"
    character(len=32), parameter :: subname = 'decompInit_clumps'
    !------------------------------------------------------------------------------
    
    !--- assign order of subgrid levels in allvecl and allvecg arrays ---
    nlev=6  ! number of subgrid levels
    glev=1  ! gridcell
    tlev=2  ! topounit
    llev=3  ! landunit
    clev=4  ! column
    plev=5  ! pft/patch
    hlev=6  ! cohort

    !--- assign gridcells to clumps (and thus pes) ---
    call get_proc_bounds(bounds)
    begg = bounds%begg; endg = bounds%endg

    allocate(allvecl(nclumps,nlev))   ! local  clumps [gcells,topounits,lunits,cols,pfts,cohs]
    allocate(allvecg(nclumps,nlev))   ! global clumps [gcells,topounits,lunits,cols,pfts,cohs]

    ! Determine the number of gridcells, topounits, landunits, columns, pfts, and cohorts 
    ! on this processor 
    ! Determine number of topounits, landunits, columns and pfts for each global
    ! gridcell index (an) that is associated with the local gridcell index (ln)
    ! More detail: an is the row-major order 1d-index into the global ixj grid.

    itunits=0
    ilunits=0
    icols=0
    ipfts=0
    icohorts=0 

    allvecg= 0
    allvecl= 0
    ! Loop through the gridcells on this proc
    do anumg = begg,endg
       ! an is the row-major order 1d-index into the global ixj grid.
       an  = gsMap_lnd_gdc2glo(anumg - begg + 1)
       cid = lcid(an)
       ln  = anumg
       if(max_topounits > 1) then
          if (present(glcmask)) then
             call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                 ncohorts=icohorts, glcmask=glcmask(ln), num_tunits_per_grd= ldomain%num_tunits_per_grd(ln))
          else
             call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                 ncohorts=icohorts, num_tunits_per_grd= ldomain%num_tunits_per_grd(ln) )
          endif
       else 
          if (present(glcmask)) then
             call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                 ncohorts=icohorts, glcmask=glcmask(ln))
          else
             call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                 ncohorts=icohorts )
          endif
       endif
       
       allvecl(cid,glev) = allvecl(cid,glev) + 1           ! number of gridcells for local clump cid
       allvecl(cid,tlev) = allvecl(cid,tlev) + itunits     ! number of topographic units for local clump cid
       allvecl(cid,llev) = allvecl(cid,llev) + ilunits     ! number of landunits for local clump cid
       allvecl(cid,clev) = allvecl(cid,clev) + icols       ! number of columns for local clump cid
       allvecl(cid,plev) = allvecl(cid,plev) + ipfts       ! number of pfts for local clump cid 
       allvecl(cid,hlev) = allvecl(cid,hlev) + icohorts    ! number of cohorts for local clump cid 
    enddo
    call mpi_allreduce(allvecl,allvecg,size(allvecg),MPI_INTEGER,MPI_SUM,mpicom,ier)

    ! Determine overall  total gridcells, landunits, columns and pfts and distribute
    ! gridcells over clumps

    numg = 0
    numt = 0
    numl = 0
    numc = 0
    nump = 0
    numCohort = 0

    do cid = 1,nclumps
       icells      = allvecg(cid,glev)  ! number of all clump cid gridcells (over all processors)
       itunits  = allvecg(cid,tlev)  ! number of all clump cid topounits (over all processors)
       ilunits     = allvecg(cid,llev)  ! number of all clump cid landunits (over all processors)
       icols       = allvecg(cid,clev)  ! number of all clump cid columns (over all processors)
       ipfts       = allvecg(cid,plev)  ! number of all clump cid pfts (over all processors)
       icohorts    = allvecg(cid,hlev)  ! number of all clump cid cohorts (over all processors)

       !--- overall total ---
       numg = numg + icells         ! total number of gridcells
       numt = numt + itunits     ! total number of landunits
       numl = numl + ilunits        ! total number of landunits
       numc = numc + icols          ! total number of columns
       nump = nump + ipfts          ! total number of pfts
       numCohort = numCohort + icohorts       ! total number of cohorts

       !--- give gridcell to cid ---
       clumps(cid)%ntunits  = clumps(cid)%ntunits  + itunits  
       clumps(cid)%nlunits     = clumps(cid)%nlunits  + ilunits  
       clumps(cid)%ncols       = clumps(cid)%ncols    + icols
       clumps(cid)%npfts       = clumps(cid)%npfts    + ipfts
       clumps(cid)%nCohorts    = clumps(cid)%nCohorts + icohorts

       do m = 1,nclumps
          if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
              (clumps(m)%owner == clumps(cid)%owner .and. m > cid)) then
             clumps(m)%begt = clumps(m)%begt + itunits
             clumps(m)%begl = clumps(m)%begl + ilunits
             clumps(m)%begc = clumps(m)%begc + icols
             clumps(m)%begp = clumps(m)%begp + ipfts
             clumps(m)%begCohort = clumps(m)%begCohort + icohorts
          endif

          if ((clumps(m)%owner >  clumps(cid)%owner) .or. &
              (clumps(m)%owner == clumps(cid)%owner .and. m >= cid)) then
             clumps(m)%endt = clumps(m)%endt + itunits
             clumps(m)%endl = clumps(m)%endl + ilunits
             clumps(m)%endc = clumps(m)%endc + icols
             clumps(m)%endp = clumps(m)%endp + ipfts
             clumps(m)%endCohort = clumps(m)%endCohort + icohorts
          endif
       enddo

       !--- give gridcell to the proc that owns the cid ---
       !--- increment the beg and end indices ---
       if (iam == clumps(cid)%owner) then
          procinfo%ntunits  = procinfo%ntunits  + itunits
          procinfo%nlunits     = procinfo%nlunits  + ilunits
          procinfo%ncols       = procinfo%ncols    + icols
          procinfo%npfts       = procinfo%npfts    + ipfts
          procinfo%nCohorts    = procinfo%nCohorts + icohorts
       endif

       if (iam >  clumps(cid)%owner) then
          procinfo%begt = procinfo%begt + itunits
          procinfo%begl = procinfo%begl + ilunits
          procinfo%begc = procinfo%begc + icols
          procinfo%begp = procinfo%begp + ipfts
          procinfo%begCohort = procinfo%begCohort + icohorts
       endif

       if (iam >= clumps(cid)%owner) then
          procinfo%endt = procinfo%endt + itunits
          procinfo%endl = procinfo%endl + ilunits
          procinfo%endc = procinfo%endc + icols
          procinfo%endp = procinfo%endp + ipfts
          procinfo%endCohort = procinfo%endCohort + icohorts
       endif
    enddo

    do n = 1,nclumps
       if (clumps(n)%ncells      /= allvecg(n,glev) .or. &
           clumps(n)%ntunits  /= allvecg(n,tlev) .or. &
           clumps(n)%nlunits     /= allvecg(n,llev) .or. &
           clumps(n)%ncols       /= allvecg(n,clev) .or. &
           clumps(n)%npfts       /= allvecg(n,plev) .or. &
           clumps(n)%nCohorts    /= allvecg(n,hlev)) then

               write(iulog,*) 'decompInit_glcp(): allvecg error ncells ',iam,n,clumps(n)%ncells ,allvecg(n,glev)
               write(iulog,*) 'decompInit_glcp(): allvecg error topounits ',iam,n,clumps(n)%ntunits,allvecg(n,tlev)
               write(iulog,*) 'decompInit_glcp(): allvecg error lunits ',iam,n,clumps(n)%nlunits,allvecg(n,llev)
               write(iulog,*) 'decompInit_glcp(): allvecg error ncols  ',iam,n,clumps(n)%ncols  ,allvecg(n,clev)
               write(iulog,*) 'decompInit_glcp(): allvecg error pfts   ',iam,n,clumps(n)%npfts  ,allvecg(n,plev)
               write(iulog,*) 'decompInit_glcp(): allvecg error cohorts ',iam,n,clumps(n)%nCohorts ,allvecg(n,hlev)

               call endrun(msg=errMsg(__FILE__, __LINE__))

       endif
    enddo

    deallocate(allvecg,allvecl)
    deallocate(lcid)

    ! Diagnostic output

    call get_proc_global(ng=numg, nl=numl, nc=numc, np=nump, nCohorts=numCohort)
    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)'   total number of topounits = ',numt
       write(iulog,*)'   total number of landunits = ',numl
       write(iulog,*)'   total number of columns   = ',numc
       write(iulog,*)'   total number of pfts      = ',nump
       write(iulog,*)'   total number of cohorts   = ',numCohort
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
       write(iulog,*)
    end if

    ! Write out clump and proc info, one pe at a time, 
    ! barrier to control pes overwriting each other on stdout

    call shr_sys_flush(iulog)
    call mpi_barrier(mpicom,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)

       if (iam == pid) then
          write(iulog,*)
          write(iulog,*)'proc= ',pid,&
               ' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
          write(iulog,*)'proc= ',pid,&
               ' beg topounit= ',procinfo%begt, &
               ' end topounit= ',procinfo%endt,                   &
               ' total topounits per proc= ',procinfo%ntunits
          write(iulog,*)'proc= ',pid,&
               ' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
          write(iulog,*)'proc= ',pid,&
               ' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
          write(iulog,*)'proc= ',pid,&
               ' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
          write(iulog,*)'proc= ',pid,&
               ' beg coh     = ',procinfo%begCohort, &
               ' end coh     = ',procinfo%endCohort,                   &
               ' total coh per proc     = ',procinfo%nCohorts
          write(iulog,*)'proc= ',pid,' nclumps = ',procinfo%nclumps

          clmin = 1
          clmax = procinfo%nclumps
          if (dbug == 1) then
            clmax = 1
          elseif (dbug == 0) then
            clmax = -1
          endif
          do n = clmin,clmax
             cid = procinfo%cid(n)
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg gridcell= ',clumps(cid)%begg, &
                  ' end gridcell= ',clumps(cid)%endg, &
                  ' total gridcells per clump= ',clumps(cid)%ncells
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg topounit= ',clumps(cid)%begt, &
                  ' end topounit= ',clumps(cid)%endt, &
                  ' total topounits per clump = ',clumps(cid)%ntunits
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg landunit= ',clumps(cid)%begl, &
                  ' end landunit= ',clumps(cid)%endl, &
                  ' total landunits per clump = ',clumps(cid)%nlunits
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg column  = ',clumps(cid)%begc, &
                  ' end column  = ',clumps(cid)%endc, &
                  ' total columns per clump  = ',clumps(cid)%ncols
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg pft     = ',clumps(cid)%begp, &
                  ' end pft     = ',clumps(cid)%endp, &
                  ' total pfts per clump     = ',clumps(cid)%npfts
             write(iulog,*)'proc= ',pid,' clump no = ',n, &
                  ' clump id= ',procinfo%cid(n),    &
                  ' beg cohort     = ',clumps(cid)%begCohort, &
                  ' end cohort     = ',clumps(cid)%endCohort, &
                  ' total cohorts per clump     = ',clumps(cid)%nCohorts
          end do
       end if
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom,ier)
    end do

  end subroutine decompInit_clumps

  !------------------------------------------------------------------------------
  subroutine decompInit_gtlcp(lni,lnj,glcmask)
    !
    ! !DESCRIPTION:
    ! Determine gsMaps for topounits, landunits, columns, pfts and cohorts
    !
    ! !USES:
    use spmdMod
    use elm_varctl             , only : use_fates
    use subgridMod,       only : subgrid_get_gcellinfo
    use decompMod              , only : bounds_type, get_proc_global, get_proc_bounds
    use decompMod              , only : gsMap_lnd_gdc2glo
    use decompMod              , only : gsMap_gce_gdc2glo, gsMap_top_gdc2glo, gsMap_lun_gdc2glo, gsMap_col_gdc2glo, gsMap_patch_gdc2glo,gsMap_cohort_gdc2glo
    use decompMod              , only : procinfo, clump_type, clumps, get_proc_global
    use TopounitType           , only : top_pp                
    use LandunitType           , only : lun_pp
    use ColumnType             , only : col_pp
    use VegetationType         , only : veg_pp                
    use FatesInterfaceTypesMod , only : fates_maxElementsPerSite
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: lni,lnj ! land domain global size
    integer , pointer, optional   :: glcmask(:)  ! glc mask
    !integer , pointer, optional   :: num_tunits_per_grd(:)  ! Number of topounits per grid
    !
    ! !LOCAL VARIABLES:
    integer              :: gi,ti,li,ci,pi,coi ! indices
    integer              :: i,g,k,l,n,np       ! indices
    integer              :: cid,pid            ! indices
    integer              :: numg               ! total number of gridcells across all processors
    integer              :: numt               ! total number of topounits across all processors
    integer              :: numl               ! total number of landunits across all processors
    integer              :: numc               ! total number of columns across all processors
    integer              :: nump               ! total number of pfts across all processors
    integer              :: numCohort          ! ED cohorts
    integer              :: itunits         ! temporary
    integer              :: ilunits            ! temporary
    integer              :: icols              ! temporary
    integer              :: ipfts              ! temporary
    integer              :: icohorts           ! temporary
    integer              :: ier                ! error code
    integer, pointer     :: gcount(:)
    integer, pointer     :: tcount(:)
    integer, pointer     :: lcount(:)
    integer, pointer     :: ccount(:)
    integer, pointer     :: pcount(:)
    integer, pointer     :: coCount(:)
    type(bounds_type)    :: bounds
    integer, allocatable :: ioff(:)
    integer, allocatable :: gridcells_per_pe(:) ! needed for gindex at all levels
    integer, allocatable :: gridcell_offsets(:) ! needed for gindex at all levels
    integer, allocatable :: index_gridcells(:)  ! needed for gindex at all levels
    integer, allocatable :: start_global(:)
    integer, allocatable :: start(:) 
    integer, allocatable :: index_lndgridcells(:)
    integer              :: count
    integer              :: temp
    integer              :: lsize_g, lsize_t, lsize_l, lsize_c, lsize_p, lsize_cohort
    integer              :: gsize
    character(len=32), parameter :: subname = 'decompInit_gtlcp'
    !------------------------------------------------------------------------------

    ! Get processor bounds

    call get_proc_bounds(bounds)
    call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump, nCohorts=numCohort)

    lsize_g = bounds%endg
    lsize_t = bounds%endt
    lsize_l = bounds%endl
    lsize_c = bounds%endc
    lsize_p = bounds%endp
    lsize_cohort = bounds%endCohort
    gsize = nglob_x * nglob_y

    ! allocate module variables in decompMod.F90
    allocate(gsMap_gce_gdc2glo(lsize_g))
    allocate(gsMap_top_gdc2glo(lsize_t))
    allocate(gsMap_lun_gdc2glo(lsize_l))
    allocate(gsMap_col_gdc2glo(lsize_c))
    allocate(gsMap_patch_gdc2glo(lsize_p))
    allocate(gsMap_cohort_gdc2glo(lsize_cohort))

    ! Determine counts
    allocate(gcount(lsize_g))  ; gcount(:) = 0
    allocate(tcount(lsize_t))  ; tcount(:) = 0
    allocate(lcount(lsize_g))  ; lcount(:) = 0
    allocate(ccount(lsize_g))  ; ccount(:) = 0
    allocate(pcount(lsize_g))  ; pcount(:) = 0
    allocate(coCount(lsize_g)) ; coCount(:) = 0
    do gi = 1,lsize_g
       if(max_topounits > 1) then
          if (present(glcmask)) then
             call subgrid_get_gcellinfo (gi, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                 ncohorts=icohorts, glcmask=glcmask(gi), num_tunits_per_grd= ldomain%num_tunits_per_grd(gi) )
          else
             call subgrid_get_gcellinfo (gi, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                  ncohorts=icohorts, num_tunits_per_grd= ldomain%num_tunits_per_grd(gi) )
          endif
       else
          if (present(glcmask)) then
             call subgrid_get_gcellinfo (gi, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                 ncohorts=icohorts, glcmask=glcmask(gi))
          else
             call subgrid_get_gcellinfo (gi, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                  ncohorts=icohorts )
          endif
       endif
       gcount(gi)  = 1          ! number of gridcells for local gridcell index gi
       tcount(gi)  = itunits ! number of topounits for local gridcell index gi
       lcount(gi)  = ilunits    ! number of landunits for local gridcell index gi
       ccount(gi)  = icols      ! number of columns for local gridcell index gi
       pcount(gi)  = ipfts      ! number of pfts for local gridcell index gi
       coCount(gi) = icohorts   ! number of ED cohorts for local gricell index gi
    enddo

    ! ---------------------------------------
    ! Arrays needed to determine gindex_xxx(:)
    ! ---------------------------------------

    allocate(ioff(lsize_g))

    if (masterproc) then
       allocate (gridcells_per_pe(0:npes-1))
    else
       allocate(gridcells_per_pe(0))
    endif
    call mpi_gather(lsize_g, 1, MPI_INTEGER, gridcells_per_pe, 1, MPI_INTEGER, 0, mpicom, ier)

    if (masterproc) then
       allocate(gridcell_offsets(0:npes-1))
       gridcell_offsets(0) = 0
       do n = 1 ,npes-1
          gridcell_offsets(n) = gridcell_offsets(n-1) + gridcells_per_pe(n-1)
       end do
    else
       allocate(gridcell_offsets(0))
    end if

    if (masterproc) then
       allocate(start_global(numg)) ! number of landunits in a gridcell
    else
       allocate(start_global(0))
    end if

    allocate(start(lsize_g))

    ! ---------------------------------------
    ! Gridcell gindex (compressed, no ocean points)
    ! ---------------------------------------

    ! gstart_global the global index of all of the land points in task order 
    call mpi_gatherv(gsMap_lnd_gdc2glo, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)

    if (masterproc) then
       ! Create a global size index_gridcells that will have 0 for all ocean points
       ! Fill the location of each land point with the gatherv location of that land point
       allocate(index_gridcells(gsize))
       index_gridcells(:) = 0
       do n = 1,numg
          ! if n = 3, start_global(3)=100, index_gridcells(100)=3
          ! n is the task order location - so for global index 100 - the task order location is 3
          index_gridcells(start_global(n)) = n
       end do

       ! Create a land-only global index based on the original global index ordering
       ! Count is the running global land index
       allocate(index_lndgridcells(numg))
       count = 0
       do n = 1,gsize
          if (index_gridcells(n) > 0) then
             count = count + 1
             ! e.g. n=20, count=4 and index_gridcells(20)=100, then start_global(100)=4
             start_global(index_gridcells(n)) = count  
             index_lndgridcells(count) = index_gridcells(n)
          end if
       end do
       deallocate(index_gridcells)
    end if

    ! Determine gsMap_gce_gdc2glo
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, gsMap_gce_gdc2glo, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)
    deallocate(gcount)

    ! ---------------------------------------
    ! topunit gindex
    ! ---------------------------------------

    start(:) = 0 
    call mpi_gatherv(tcount, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
    if (masterproc) then
       count = 1
       do n = 1,numg
          temp = start_global(index_lndgridcells(n))
          start_global(index_lndgridcells(n)) = count
          count = count + temp
       end do
    endif
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)

    ioff(:) = 0
    do ti = 1,lsize_t
       gi = top_pp%gridcell(ti)
       gsMap_top_gdc2glo(ti) = start(gi) + ioff(gi)
       ioff(gi)  = ioff(gi) + 1
    enddo
    deallocate(tcount)

    ! ---------------------------------------
    ! Landunit gindex
    ! ---------------------------------------

    start(:) = 0 
    call mpi_gatherv(lcount, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
    if (masterproc) then
       count = 1
       do n = 1,numg
          temp = start_global(index_lndgridcells(n))
          start_global(index_lndgridcells(n)) = count
          count = count + temp
       end do
    endif
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)

    ioff(:) = 0
    do li = 1,lsize_l
       gi = lun_pp%gridcell(li)
       gsMap_lun_gdc2glo(li) = start(gi) + ioff(gi)
       ioff(gi)  = ioff(gi) + 1
    enddo
    deallocate(lcount)

    ! ---------------------------------------
    ! Column gindex
    ! ---------------------------------------

    start(:) = 0
    call mpi_gatherv(ccount, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
    if (masterproc) then
       count = 1
       do n = 1,numg
          temp = start_global(index_lndgridcells(n))
          start_global(index_lndgridcells(n)) = count
          count = count + temp
       end do
    endif
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)

    ioff(:) = 0
    do ci = 1,lsize_c
       gi = col_pp%gridcell(ci)
       gsMap_col_gdc2glo(ci) = start(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1
    enddo
    deallocate(ccount)

    ! ---------------------------------------
    ! PATCH gindex
    ! ---------------------------------------

    start(:) = 0
    call mpi_gatherv(pcount, lsize_g, MPI_INTEGER, start_global, &
         gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
    if (masterproc) then
       count = 1
       do n = 1,numg
          temp = start_global(index_lndgridcells(n))
          start_global(index_lndgridcells(n)) = count
          count = count + temp
       end do
    endif
    call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
         lsize_g, MPI_INTEGER, 0, mpicom, ier)

    ioff(:) = 0
    do pi = 1,lsize_p
       gi = veg_pp%gridcell(pi)
       gsMap_patch_gdc2glo(pi) = start(gi) + ioff(gi)
       ioff(gi) = ioff(gi) + 1
    enddo
    deallocate(pcount)

    ! ---------------------------------------
    ! FATES gindex for the cohort/element vector
    ! ---------------------------------------

    if ( use_fates ) then
       start(:) = 0
       call mpi_gatherv(coCount, lsize_g, MPI_INTEGER, start_global, &
            gridcells_per_pe, gridcell_offsets, MPI_INTEGER, 0, mpicom, ier)
       if (masterproc) then
          count = 1
          do n = 1,numg
             temp = start_global(index_lndgridcells(n))
             start_global(index_lndgridcells(n)) = count
             count = count + temp
          end do
       endif
       call mpi_scatterv(start_global, gridcells_per_pe, gridcell_offsets, MPI_INTEGER, start, &
            lsize_g, MPI_INTEGER, 0, mpicom, ier)

       ioff(:) = 0
       gi = 1
       do coi = 1, lsize_cohort
          gsMap_cohort_gdc2glo(coi) = start(gi) + ioff(gi)
          ioff(gi) = ioff(gi) + 1
          if ( mod(coi, fates_maxElementsPerSite ) == 0 ) then
             gi = gi + 1
          end if
       enddo
       deallocate(coCount)
    endif

    ! ---------------------------------------
    ! Deallocate memory
    ! ---------------------------------------

    deallocate(ioff)
    deallocate(gridcells_per_pe)
    deallocate(gridcell_offsets)
    deallocate(start)
    deallocate(start_global)
    if (allocated(index_lndgridcells)) deallocate(index_lndgridcells)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points          = ',lni
       write(iulog,*)'   latitude points           = ',lnj
       write(iulog,*)'   total number of gridcells = ',numg
       write(iulog,*)'   total number of topounits = ',numt
       write(iulog,*)'   total number of landunits = ',numl
       write(iulog,*)'   total number of columns   = ',numc
       write(iulog,*)'   total number of pfts      = ',nump
       write(iulog,*)'   total number of cohorts   = ',numCohort
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process        = ',clump_pproc
    end if

    call shr_sys_flush(iulog)

  end subroutine decompInit_gtlcp

  !------------------------------------------------------------------------------
  subroutine decompInit_lnd_using_gp(lni, lnj, cellsOnCell, ncells_loc, maxEdges, &
       amask)
    !
    ! !DESCRIPTION:
    ! This subroutine initializes the land surface decomposition into a clump
    ! data structure using graph partitioning approach.  This assumes each pe
    ! has the same number of clumps set by clump_pproc.
    !
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
    use elm_varctl, only : nsegspc
#ifdef USE_PETSC_LIB
    use petscsys
    use petscvec
    use petscmat
    use petscdm
#endif
    !
    ! !ARGUMENTS:
    implicit none
    !
    !
    integer , intent(in) :: amask(:)
    integer , intent(in) :: lni,lnj                     ! domain global size
    integer , intent(in) :: cellsOnCell(:,:)
    integer , intent(in) :: ncells_loc
    integer , intent(in) :: maxEdges
    !
    ! !LOCAL VARIABLES:
    integer            :: lns                           ! global domain size
    integer            :: ln,lj                         ! indices
    integer            :: numg                          ! number of land gridcells
    integer            :: cid,pid                       ! indices
    integer            :: n,i,m                         ! indices
    integer            :: ier                           ! error code
    integer            :: beg,end,lsize,gsize           ! used for gsmap init
    integer, pointer   :: gindex(:)                     ! global index for gsmap init
    integer            :: icell, iedge                  ! indices
    integer            :: offset                        ! temporary
    integer            :: cell_id_offset                ! temporary
    integer            :: count                         ! temporary
    integer            :: num_rows, num_cols            ! temporary
    integer            :: istart, iend                  ! temporary
    integer            :: ncells_tot                    ! total number of grid cells
    integer            :: ncells_owned                  ! number of grid cells owned by a processor after domain decomposition
    integer            :: ncells_per_clump              ! number of grid cells per clump
    integer            :: remainder                     ! temporary
    integer            :: cowner                        ! clump owner
    integer, pointer   :: i_index(:)                    ! temporary
    integer, pointer   :: j_index(:)                    ! temporary
    integer, pointer   :: local_conn_offsets(:)         ! temporary
    integer, pointer   :: local_conn(:)                 ! temporary
    integer, pointer   :: clump_ncells(:)               ! temporary
    integer, pointer   :: clump_begg(:)                 ! temporary
    integer, pointer   :: clump_endg(:)                 ! temporary
    integer, pointer   :: local_clump_info(:)           ! temporary
    integer, pointer   :: global_clump_info(:)          ! temporary
    integer, pointer   :: thread_count(:)               ! temporary
    integer, pointer   :: int_array(:)                  ! temporary
    integer, pointer   :: ncells_count(:)               ! temporary
#ifdef USE_PETSC_LIB
    PetscReal, pointer :: real_ptr(:)                   ! temporary
    PetscInt, pointer  :: int_ptr(:)                    ! temporary
    PetscBool          :: success                       ! temporary
    VecScatter         :: scatter                       ! temporary
    Vec                :: ids_old                       ! grid cell IDs before domain decomposition
    Vec                :: ids_new                       ! grid cell IDs after domain decomposition
    Vec                :: lcid_aft_decomp               ! local clump ID for grid cells owned after domain decomposition
    Vec                :: lcid_aft_decomp_for_all_procs ! local clump ID for all grid cells after domain decomposition
    Mat                :: Dual_mat                      ! dual matrix
    Mat                :: Dual_aij                      ! dual matrix in aij format
    Mat                :: Dual_aij_loc                  ! local part of the dual matrix in aij format
    MatPartitioning    :: part                          ! partitioning of dual matrix
    IS                 :: is_new_owner, is_new_id       ! index set that stores mpi rank and id of grid cell after domain decomposition
    IS                 :: is_from, is_to                ! temporary
    PetscErrorCode     :: ierr                          ! get error code from PETSc
#endif
    character(len=255) :: subname = 'decompInit_lnd_using_gp'
    !------------------------------------------------------------------------------

#ifndef USE_PETSC_LIB

    call endrun(msg='ERROR ' // trim(subname) //': Graph partitioning requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

#else

    lns = lni * lnj

    !--- set and verify nclumps ---
    if (clump_pproc > 0) then
       nclumps = clump_pproc * npes
       if (nclumps < npes) then
          write(iulog,*) trim(subname) // '(): Number of gridcell clumps= ',nclumps, &
               ' is less than the number of processes = ', npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    else
       write(iulog,*)trim(subname) // '(): clump_pproc= ',clump_pproc,'  must be greater than 0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! allocate and initialize procinfo and clumps
    ! beg and end indices initialized for simple addition of cells later

    allocate(procinfo%cid(clump_pproc), stat=ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname) // '(): allocation error for procinfo%cid'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    procinfo%nclumps   = clump_pproc
    procinfo%cid(:)    = -1
    procinfo%ncells    = 0
    procinfo%nlunits   = 0
    procinfo%ncols     = 0
    procinfo%npfts     = 0
    procinfo%nCohorts  = 0
    procinfo%begg      = 1
    procinfo%begl      = 1
    procinfo%begc      = 1
    procinfo%begp      = 1
    procinfo%begCohort = 1
    procinfo%endg      = 0
    procinfo%endl      = 0
    procinfo%endc      = 0
    procinfo%endp      = 0
    procinfo%endCohort = 0

    allocate(clumps(nclumps), stat=ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname) // '(): allocation error for clumps'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    clumps(:)%owner     = -1
    clumps(:)%ncells    = 0
    clumps(:)%nlunits   = 0
    clumps(:)%ncols     = 0
    clumps(:)%npfts     = 0
    clumps(:)%nCohorts  = 0
    clumps(:)%begg      = 1
    clumps(:)%begl      = 1
    clumps(:)%begc      = 1
    clumps(:)%begp      = 1
    clumps(:)%begCohort = 1
    clumps(:)%endg      = 0
    clumps(:)%endl      = 0
    clumps(:)%endc      = 0
    clumps(:)%endp      = 0
    clumps(:)%endCohort = 0

    ! assign clumps to proc round robin
    cid = 0
    do n = 1,nclumps
       pid = mod(n-1,npes)
       if (pid < 0 .or. pid > npes-1) then
          write(iulog,*) trim(subname) // '(): round robin pid error ',n,pid,npes
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
       clumps(n)%owner = pid
       if (iam == pid) then
          cid = cid + 1
          if (cid < 1 .or. cid > clump_pproc) then
             write(iulog,*) trim(subname) // '(): round robin pid error ',n,pid,npes
             call endrun(msg=errMsg(_FILE__, __LINE__))
          endif
          procinfo%cid(cid) = n
       endif
    enddo

    ! count total active land gridcells
    numg = 0
    do ln = 1,lns
       if (amask(ln) == 1) then
          numg = numg + 1
       endif
    enddo

    if (npes > numg) then
       write(iulog,*) trim(subname) // '(): Number of processes exceeds number ', &
            'of land grid cells',npes,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (nclumps > numg) then
       write(iulog,*) trim(subname) // '(): Number of clumps exceeds number ', &
            'of land grid cells',nclumps,numg
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (numg /= lns) then
       write(iulog,*) trim(subname) // '(): Only implimented for numg == lns '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Determine the cell id offset on each processor
    cell_id_offset = 0
    call MPI_Scan(ncells_loc, cell_id_offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)
    cell_id_offset = cell_id_offset - ncells_loc

    ! Determine the total number of grid cells
    call MPI_Allreduce(ncells_loc, ncells_tot, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)

    ! Create a Dual matrix.
    call MatCreateAIJ(mpicom , & ! comm
         ncells_loc          , & ! m
         PETSC_DECIDE        , & ! n
         PETSC_DETERMINE     , & ! M
         ncells_tot          , & ! N
         maxEdges            , & ! d_nz
         PETSC_NULL_INTEGER  , & ! d_nnz
         maxEdges            , & ! o_nz
         PETSC_NULL_INTEGER  , & ! o_nnz
         Dual_aij            , & ! Mat
         ierr);CHKERRQ(ierr)

    !
    ! If cell_1 and cell_2 are neighbors, then set:
    !
    !  Dual_aij(cell_1, cell_2) = 1
    !  Dual_aij(cell_2, cell_1) = 1
    !
    do icell = 1, ncells_loc
       do iedge = 1, maxEdges

          if (cellsOnCell(iedge, icell) > 0) then

             call MatSetValue(Dual_aij        , &
                  icell+cell_id_offset-1      , &
                  cellsOnCell(iedge, icell)-1 , &
                  1.d0                        , &
                  INSERT_VALUES               , &
                  ierr);CHKERRQ(ierr)

             call MatSetValue(Dual_aij        , &
                  cellsOnCell(iedge, icell)-1 , &
                  icell+cell_id_offset-1      , &
                  1.d0                        , &
                  INSERT_VALUES               , &
                  ierr);CHKERRQ(ierr)
          end if
       end do
    end do

    call MatAssemblyBegin( Dual_aij,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(   Dual_aij,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

    !
    ! Create sparse matrix representing an adjacency list from Dual matrix
    !

    if ( npes > 1 ) then

       call MatMPIAIJGetLocalMat(Dual_aij , &
            MAT_INITIAL_MATRIX            , &
            Dual_aij_loc                  , &
            ierr);CHKERRQ(ierr)

       call MatGetRowIJF90(Dual_aij_loc   , &
            0                             , &
            PETSC_FALSE                   , &
            PETSC_FALSE                   , &
            num_rows                      , &
            i_index                       , &
            j_index                       , &
            success                       , &
            ierr);CHKERRQ(ierr)
    else

       call MatGetRowIJF90(Dual_aij       , &
            0                             , &
            PETSC_FALSE                   , &
            PETSC_FALSE                   , &
            num_rows                      , &
            i_index                       , &
            j_index                       , &
            success                       , &
            ierr);CHKERRQ(ierr)
    endif

    count = 0
    do icell = 1,num_rows
       istart   = i_index(icell)
       iend     = i_index(icell + 1) - 1
       num_cols = iend - istart + 1
       count    = count + num_cols
    enddo

    allocate (local_conn         (count      ))
    allocate (local_conn_offsets (num_rows+1 ))

    local_conn_offsets (1:num_rows+1) = i_index (1:num_rows+1)
    local_conn         (1:count)      = j_index (1:count     )

    call MatCreateMPIAdj(mpicom , &
         ncells_loc             , &
         ncells_tot             , &
         local_conn_offsets     , &
         local_conn             , &
         PETSC_NULL_INTEGER     , &
         Dual_mat               , &
         ierr); CHKERRQ(ierr)

    if ( npes > 1 ) then
       call MatRestoreRowIJF90(Dual_aij_loc , &
            0                               , &
            PETSC_FALSE                     , &
            PETSC_FALSE                     , &
            num_rows                        , &
            i_index                         , &
            j_index                         , &
            success                         , &
            ierr);CHKERRQ(ierr)
    else
       call MatGetRowIJF90(Dual_aij , &
            0                       , &
            PETSC_FALSE             , &
            PETSC_FALSE             , &
            num_rows                , &
            i_index                 , &
            j_index                 , &
            success                 , &
            ierr);CHKERRQ(ierr)
    endif
    call MatDestroy(Dual_aij,ierr);CHKERRQ(ierr)

    !
    ! Use graph partitioning to decompose the mesh
    !

    call MatPartitioningCreate(mpicom, part, ierr);CHKERRQ(ierr)

    call MatPartitioningSetAdjacency(part, Dual_mat, ierr);CHKERRQ(ierr)

    call MatPartitioningSetFromOptions(part, ierr);CHKERRQ(ierr)

    ! Now perform graph partioning.
    !   - After the call to MatPartitioningApply(),
    !     is_new_owner has entries that corresponds to the rank of
    !     processor which owns each grid cell
    call MatPartitioningApply(part, is_new_owner, ierr);CHKERRQ(ierr)

    ! Free up memory
    call MatDestroy(Dual_mat, ierr); CHKERRQ(ierr);
    call MatPartitioningDestroy(part, ierr);CHKERRQ(ierr)
    deallocate(local_conn        )
    deallocate(local_conn_offsets)

    ! Determine the number of grid cells owned after graph partitioning
    allocate(ncells_count(npes))
    call ISPartitioningCount(is_new_owner, npes, ncells_count, ierr);CHKERRQ(ierr)
    ncells_owned = ncells_count(iam + 1)
    deallocate(ncells_count)

    ! Determine the new ids of grid cells
    call ISPartitioningToNumbering(is_new_owner, is_new_id, ierr);CHKERRQ(ierr)

    !
    ! Compute information for each processor
    !
    procinfo%ncells = ncells_owned

    offset = 0
    call MPI_Scan(ncells_owned, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)
    procinfo%begg = offset + 1 - ncells_owned

    offset = 0
    call MPI_scan(ncells_owned, offset, 1, MPIU_INTEGER, MPI_SUM, mpicom, ierr)
    procinfo%endg = offset

    !
    ! Compute information about all clumps on all processors
    !

    ncells_per_clump = procinfo%ncells/clump_pproc
    remainder        = procinfo%ncells - ncells_per_clump*clump_pproc

    allocate (clump_ncells      (clump_pproc            ))
    allocate (clump_begg        (clump_pproc            ))
    allocate (clump_endg        (clump_pproc            ))

    allocate (local_clump_info  (0:3*clump_pproc-1      ))
    allocate (global_clump_info (0:3*clump_pproc*npes-1 ))

    allocate (thread_count      (0:npes-1               ))

    clump_ncells(:) = ncells_per_clump

    offset = procinfo%begg
    do m = 1,clump_pproc
       if (m-1 < remainder) clump_ncells(m) = clump_ncells(m) + 1

       clump_begg(m) = offset
       clump_endg(m) = offset + clump_ncells(m) - 1
       offset        = offset + clump_ncells(m)

       local_clump_info((m-1)*3 + 0) = clump_ncells(m)
       local_clump_info((m-1)*3 + 1) = clump_begg  (m)
       local_clump_info((m-1)*3 + 2) = clump_endg  (m)
    end do

    call MPI_Allgather(local_clump_info, 3*clump_pproc, MPI_INTEGER, &
         global_clump_info, 3*clump_pproc, MPI_INTEGER, mpicom, ier)

    count = 0
    thread_count(:) = 0
    do m = 1, nclumps
       cowner = clumps(m)%owner

       clumps(m)%ncells = global_clump_info(cowner*3 + thread_count(cowner)*3    )
       clumps(m)%begg   = global_clump_info(cowner*3 + thread_count(cowner)*3 + 1)
       clumps(m)%endg   = global_clump_info(cowner*3 + thread_count(cowner)*3 + 2)

       thread_count(cowner) = thread_count(cowner) + 1
    enddo

    deallocate (clump_ncells      )
    deallocate (clump_begg        )
    deallocate (clump_endg        )
    deallocate (local_clump_info  )
    deallocate (global_clump_info )
    deallocate (thread_count      )

    !
    ! Determine the natural ids of the grid cells that each processor owns after
    ! domain decomposition. This information will be used to set up ldecomp%gdc2glo(:).
    !
    call VecCreateMPI(mpicom, ncells_loc  , PETSC_DETERMINE, ids_old, ierr);CHKERRQ(ierr)
    call VecCreateMPI(mpicom, ncells_owned, PETSC_DETERMINE, ids_new, ierr);CHKERRQ(ierr)

    ! Create a vector containing old ids of grid cells (i.e. ids before
    ! domain decomposition)
    call VecGetArrayF90(ids_old, real_ptr, ierr)
    do i = 1, ncells_loc
       real_ptr(i) = i + cell_id_offset
    end do
    call VecRestoreArrayF90(ids_old, real_ptr, ierr)

    ! Create an IS to scatter data stored in ids_old
    allocate(int_array(ncells_loc))
    do m = 1, ncells_loc
       int_array(m) = m - 1 + cell_id_offset
    end do

    call ISCreateGeneral(mpicom, ncells_loc, int_array, &
         PETSC_COPY_VALUES, is_from, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    ! Create a VectorScatter
    call VecScatterCreate(ids_old, is_from, ids_new, is_new_id, &
         scatter, ierr); CHKERRQ(ierr);

    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_new_id, ierr); CHKERRQ(ierr)

    ! Scatter data to get new ids of grid cells (i.e. ids after
    ! domain decomposition)
    call VecScatterBegin(scatter, ids_old, ids_new, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(scatter, ids_old, ids_new, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

    ! Set ldecomp

    call get_proc_bounds(beg, end)
    allocate(ldecomp%gdc2glo(beg:end), stat=ier)
     if (ier /= 0) then
       write(iulog,*) trim(subname) // '(): allocation error1 for ldecomp, etc'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ldecomp%gdc2glo(:) = 0

    call VecGetArrayF90(ids_new, real_ptr, ierr); CHKERRQ(ierr);
    do m = beg, end
       ldecomp%gdc2glo(m) = INT(real_ptr(m-beg+1))
    end do
    call VecRestoreArrayF90(ids_new, real_ptr, ierr); CHKERRQ(ierr);

    call VecDestroy(ids_old, ierr); CHKERRQ(ierr);
    call VecDestroy(ids_new, ierr); CHKERRQ(ierr);

    !
    ! For each grid cell, identify the processor that owns
    ! the grid cell after domain decomposition.
    !
    ! NOTE: This is not a scalable approach as each processor
    !       stores the information about the entire domain.
    !       But, decompInit_clumps() needs this information.
    !

    call VecCreateMPI(mpicom           , &
         ncells_loc                    , &
         PETSC_DETERMINE               , &
         lcid_aft_decomp               , &
         ierr); CHKERRQ(ierr)

    call VecCreateMPI(mpicom           , &
         numg                          , &
         PETSC_DETERMINE               , &
         lcid_aft_decomp_for_all_procs , &
         ierr); CHKERRQ(ierr)

    call ISGetIndicesF90(is_new_owner   , int_ptr , ierr); CHKERRQ(ierr)
    call VecGetArrayF90( lcid_aft_decomp, real_ptr, ierr); CHKERRQ(ierr)
    real_ptr(:) = real(int_ptr(:) + 1._r8)
    call ISRestoreIndicesF90(is_new_owner   , int_ptr , ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90( lcid_aft_decomp, real_ptr, ierr); CHKERRQ(ierr)

    allocate(int_array(numg))
    do m = 1, numg
       int_array(m) = m - 1
    end do
    call ISCreateGeneral(mpicom, numg, int_array, &
         PETSC_COPY_VALUES, is_from, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    allocate(int_array(numg))
    do m = 1, numg
       int_array(m) = m - 1 + iam*numg
    end do
    call ISCreateGeneral(mpicom, numg, int_array, &
         PETSC_COPY_VALUES, is_to, ierr);CHKERRQ(ierr);
    deallocate(int_array)

    call VecScatterCreate(lcid_aft_decomp, is_from, &
         lcid_aft_decomp_for_all_procs, is_to, scatter, ierr);CHKERRQ(ierr)

    call ISDestroy(is_from, ierr); CHKERRQ(ierr)
    call ISDestroy(is_to  , ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, lcid_aft_decomp, lcid_aft_decomp_for_all_procs, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(scatter, lcid_aft_decomp, lcid_aft_decomp_for_all_procs, &
         INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)

    ! Set lcid

    allocate(lcid(lns))
    lcid(:) = 0

    call VecGetArrayF90(lcid_aft_decomp_for_all_procs, real_ptr, ierr); CHKERRQ(ierr)
    lcid(:) = INT(real_ptr(:))
    call VecRestoreArrayF90(lcid_aft_decomp_for_all_procs, real_ptr, ierr); CHKERRQ(ierr)

    call ISDestroy(is_new_owner, ierr);CHKERRQ(ierr)

    ! Set gsMap_lnd_gdc2glo (the global index here includes mask=0 or ocean points)

    call get_proc_bounds(beg, end)
    allocate(gindex(beg:end))
    do n = beg,end
       gindex(n) = ldecomp%gdc2glo(n)
    enddo
    lsize = end-beg+1
    gsize = lni * lnj

    deallocate(gindex)

    ! Diagnostic output

    if (masterproc) then
       write(iulog,*)' Surface Grid Characteristics'
       write(iulog,*)'   longitude points               = ',lni
       write(iulog,*)'   latitude points                = ',lnj
       write(iulog,*)'   total number of land gridcells = ',numg
       write(iulog,*)' Decomposition Characteristics'
       write(iulog,*)'   clumps per process             = ',clump_pproc
       write(iulog,*)
    end if

    call shr_sys_flush(iulog)

#endif

  end subroutine decompInit_lnd_using_gp

  !------------------------------------------------------------------------------
!hh!  subroutine decompInit_ghosts(glcmask)
  subroutine decompInit_ghosts(lni, lnj, glcmask)
    !
    ! !DESCRIPTION:
    ! On each proc, determine the number of following ghost/halo subgrid quantities:
    !  - topounits
    !  - landunits,
    !  - columns,
    !  - PFTs, and
    !  - cohorts.
    !
    ! !USES:
    use elm_varctl           , only : lateral_connectivity
    use subgridMod           , only : subgrid_get_gcellinfo
#ifdef USE_PETSC_LIB
    use domainLateralMod     , only : ldomain_lateral
    use UnstructuredGridType , only : ScatterDataG2L
#endif
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: lni,lnj ! land domain global size
    integer , pointer, optional  :: glcmask(:)             ! glc mask
    !integer , pointer, optional  :: num_tunits_per_grd(:)             ! glc mask
    !
    ! !LOCAL VARIABLES:
    integer                      :: begg,endg              ! begin/end indices for grid
    integer                      :: anumg                  ! lnd num gridcells
    integer                      :: ln                     ! temporary
    integer                      :: nblocks                ! block size for PETSc vector
    integer                      :: itunits             ! temporary
    integer                      :: ilunits                ! temporary
    integer                      :: icols                  ! temporary
    integer                      :: ipfts                  ! temporary
    integer                      :: icohorts               ! temporary
    integer                      :: ighost                 ! temporary
    integer                      :: ighost_beg, ighost_end ! temporary
    real(r8), pointer            :: data_send(:)
    real(r8), pointer            :: data_recv(:)
    integer                      :: ndata_send
    integer                      :: ndata_recv
    character(len=32), parameter :: subname = 'decompInit_ghosts'

    if (.not.lateral_connectivity) then

       ! No ghost cells
       procinfo%ncells_ghost    = 0
       procinfo%ntunits_ghost   = 0
       procinfo%nlunits_ghost   = 0
       procinfo%ncols_ghost     = 0
       procinfo%npfts_ghost     = 0
       procinfo%nCohorts_ghost  = 0

       procinfo%begg_ghost      = 0
       procinfo%begt_ghost      = 0
       procinfo%begl_ghost      = 0
       procinfo%begc_ghost      = 0
       procinfo%begp_ghost      = 0
       procinfo%begCohort_ghost = 0
       procinfo%endg_ghost      = 0
       procinfo%endt_ghost      = 0
       procinfo%endl_ghost      = 0
       procinfo%endc_ghost      = 0
       procinfo%endp_ghost      = 0
       procinfo%endCohort_ghost = 0

       ! All = local (as no ghost cells)
       procinfo%ncells_all      = procinfo%ncells
       procinfo%ntunits_all     = procinfo%ntunits
       procinfo%nlunits_all     = procinfo%nlunits
       procinfo%ncols_all       = procinfo%ncols
       procinfo%npfts_all       = procinfo%npfts
       procinfo%nCohorts_all    = procinfo%nCohorts

       procinfo%begg_all        = procinfo%begg
       procinfo%begt_all        = procinfo%begt
       procinfo%begl_all        = procinfo%begl
       procinfo%begc_all        = procinfo%begc
       procinfo%begp_all        = procinfo%begp
       procinfo%begCohort_all   = procinfo%begCohort
       procinfo%endg_all        = procinfo%endg
       procinfo%endt_all        = procinfo%endt
       procinfo%endl_all        = procinfo%endl
       procinfo%endc_all        = procinfo%endc
       procinfo%endp_all        = procinfo%endp
       procinfo%endCohort_all   = procinfo%endCohort

    else

#ifndef USE_PETSC_LIB

    call endrun(msg='ERROR ' // trim(subname) //': decompInit_ghosts requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

#else
       call get_proc_bounds(begg, endg)

       ! Approach:
       ! 1) For a global PETSc vector, save the number of subgrid
       !    quantities for each grid cell.
       ! 2) Scatter the global PETSc vector to a local PETSc vector
       ! 3) Finally count the number of subgrid quantities for all
       !    ghost grid cells in the local PETSc vector

       nblocks = 5 ! topo + lun + col + pft + cohort

       ndata_send = nblocks*ldomain_lateral%ugrid%ngrid_local
       ndata_recv = nblocks*ldomain_lateral%ugrid%ngrid_ghosted

       allocate(data_send(ndata_send))
       allocate(data_recv(ndata_recv))

       data_send(:) = 0.d0

       ! Save information about number of subgrid categories for
       ! local grid cells

       do anumg = begg,endg
          ln  = anumg
          if(max_topounits > 1) then
             if (present(glcmask)) then
                call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                     ncohorts=icohorts, glcmask=glcmask(ln), num_tunits_per_grd= ldomain%num_tunits_per_grd(ln) )
             else
                call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                     ncohorts=icohorts, num_tunits_per_grd= ldomain%num_tunits_per_grd(ln) )
             endif
          else
             if (present(glcmask)) then
                call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                     ncohorts=icohorts, glcmask=glcmask(ln))
             else
                call subgrid_get_gcellinfo (ln, ntunits=itunits, nlunits=ilunits, ncols=icols, npfts=ipfts, &
                     ncohorts=icohorts )
             endif
          endif 
          
          data_send((anumg-begg)*nblocks + 1) = itunits
          data_send((anumg-begg)*nblocks + 2) = ilunits
          data_send((anumg-begg)*nblocks + 3) = icols
          data_send((anumg-begg)*nblocks + 4) = ipfts
          data_send((anumg-begg)*nblocks + 5) = icohorts

       enddo

       ! Scatter: Global-to-Local
       call ScatterDataG2L(ldomain_lateral%ugrid, &
            nblocks, ndata_send, data_send, ndata_recv, data_recv)

       ! Get number of ghost quantites at all subgrid categories
       procinfo%ncells_ghost    = ldomain_lateral%ugrid%ngrid_ghost
       procinfo%ntunits_ghost   = 0
       procinfo%nlunits_ghost   = 0
       procinfo%ncols_ghost     = 0
       procinfo%npfts_ghost     = 0
       procinfo%nCohorts_ghost  = 0

       ighost_beg = ldomain_lateral%ugrid%ngrid_local   + 1
       ighost_end = ldomain_lateral%ugrid%ngrid_ghosted

       do ighost = ighost_beg, ighost_end
          procinfo%ntunits_ghost  = procinfo%ntunits_ghost  + data_recv((ighost-1)*nblocks + 1)
          procinfo%nlunits_ghost  = procinfo%nlunits_ghost  + data_recv((ighost-1)*nblocks + 2)
          procinfo%ncols_ghost    = procinfo%ncols_ghost    + data_recv((ighost-1)*nblocks + 3)
          procinfo%npfts_ghost    = procinfo%npfts_ghost    + data_recv((ighost-1)*nblocks + 4)
          procinfo%ncohorts_ghost = procinfo%ncohorts_ghost + data_recv((ighost-1)*nblocks + 5)
       enddo

       ! Set 'begin' index for subgrid categories
       procinfo%begg_all        = procinfo%begg
       procinfo%begt_all        = procinfo%begt
       procinfo%begl_all        = procinfo%begl
       procinfo%begc_all        = procinfo%begc
       procinfo%begp_all        = procinfo%begp
       procinfo%begCohort_all   = procinfo%begCohort

       ! Set 'end' index for subgrid categories
       procinfo%endg_all        = procinfo%endg      + procinfo%ncells_ghost
       procinfo%endt_all        = procinfo%endt      + procinfo%ntunits_ghost
       procinfo%endl_all        = procinfo%endl      + procinfo%nlunits_ghost
       procinfo%endc_all        = procinfo%endc      + procinfo%ncols_ghost
       procinfo%endp_all        = procinfo%endp      + procinfo%npfts_ghost
       procinfo%endCohort_all   = procinfo%endCohort + procinfo%nCohorts_ghost

#endif

    endif

  end subroutine decompInit_ghosts

end module decompInitMod
