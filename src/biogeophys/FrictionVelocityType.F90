module FrictionVelocityType

  !------------------------------------------------------------------------------
  !
  ! !USES:
  use shr_sys_mod    , only : shr_sys_flush
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use elm_varctl     , only : use_cn, use_fates, iulog
  use elm_varpar     , only : nlevcan, nlevsno, nlevgrnd, nlevsoi
  use elm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc
  use LandunitType   , only : lun_pp
  use ColumnType     , only : col_pp
  use VegetationType , only : veg_pp
  use CanopyStateType, only : canopystate_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  type, public :: frictionvel_type
     real(r8), public :: zetamaxstable = -999._r8  ! Max value zeta ("height" used in Monin-Obukhov theory) can go to under stable conditions

     ! Roughness length/resistance for friction velocity calculation

     real(r8), pointer :: forc_hgt_u_patch (:)   ! patch wind forcing height (10m+z0m+d) (m)
     real(r8), pointer :: forc_hgt_t_patch (:)   ! patch temperature forcing height (10m+z0m+d) (m)
     real(r8), pointer :: forc_hgt_q_patch (:)   ! patch specific humidity forcing height (10m+z0m+d) (m)
     real(r8), pointer :: u10_patch        (:)   ! patch 10-m wind (m/s) (for dust model)
     real(r8), pointer :: u10_elm_patch    (:)   ! patch 10-m wind (m/s) (for elm_map2gcell)
     real(r8), pointer :: va_patch         (:)   ! patch atmospheric wind speed plus convective velocity (m/s)
     real(r8), pointer :: vds_patch        (:)   ! patch deposition velocity term (m/s) (for dry dep SO4, NH4NO3)
     real(r8), pointer :: fv_patch         (:)   ! patch friction velocity (m/s) (for dust model)
     real(r8), pointer :: rb1_patch        (:)   ! patch aerodynamical resistance (s/m) (for dry deposition of chemical tracers)
     real(r8), pointer :: ram1_patch       (:)   ! patch aerodynamical resistance (s/m)
     real(r8), pointer :: z0m_patch        (:)   ! patch momentum roughness length (m)
     real(r8), pointer :: z0mv_patch       (:)   ! patch roughness length over vegetation, momentum [m]
     real(r8), pointer :: z0hv_patch       (:)   ! patch roughness length over vegetation, sensible heat [m]
     real(r8), pointer :: z0qv_patch       (:)   ! patch roughness length over vegetation, latent heat [m]
     real(r8), pointer :: z0mg_col         (:)   ! col roughness length over ground, momentum  [m]
     real(r8), pointer :: z0hg_col         (:)   ! col roughness length over ground, sensible heat [m]
     real(r8), pointer :: z0qg_col         (:)   ! col roughness length over ground, latent heat [m]
     real(r8), pointer, public :: z0m_actual_patch (:)   ! patch roughness length actually used in flux calculations, momentum [m]

   contains

     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: SetActualRoughnessLengths ! Set roughness lengths actually used in flux calculations
     procedure, public  :: MoninObukIni           ! Initialization of the Monin-Obukhov length
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold
     procedure, private :: ReadNamelist

  end type frictionvel_type
  !------------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
!hh!  subroutine Init(this, bounds)
  subroutine Init(this, bounds, NLFilename)

    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: NLFilename  ! file name of namelist file

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)
    call this%ReadNamelist(NLFilename)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%forc_hgt_u_patch (begp:endp)) ; this%forc_hgt_u_patch (:)   = spval
    allocate(this%forc_hgt_t_patch (begp:endp)) ; this%forc_hgt_t_patch (:)   = spval
    allocate(this%forc_hgt_q_patch (begp:endp)) ; this%forc_hgt_q_patch (:)   = spval
    allocate(this%u10_patch        (begp:endp)) ; this%u10_patch        (:)   = spval
    allocate(this%u10_elm_patch    (begp:endp)) ; this%u10_elm_patch    (:)   = spval
    allocate(this%va_patch         (begp:endp)) ; this%va_patch         (:)   = spval
    allocate(this%vds_patch        (begp:endp)) ; this%vds_patch        (:)   = spval
    allocate(this%fv_patch         (begp:endp)) ; this%fv_patch         (:)   = spval
    allocate(this%rb1_patch        (begp:endp)) ; this%rb1_patch        (:)   = spval
    allocate(this%ram1_patch       (begp:endp)) ; this%ram1_patch       (:)   = spval
    allocate(this%z0m_patch        (begp:endp)) ; this%z0m_patch        (:)   = spval
    allocate(this%z0mv_patch       (begp:endp)) ; this%z0mv_patch       (:)   = spval
    allocate(this%z0hv_patch       (begp:endp)) ; this%z0hv_patch       (:)   = spval
    allocate(this%z0qv_patch       (begp:endp)) ; this%z0qv_patch       (:)   = spval
    allocate(this%z0mg_col         (begc:endc)) ; this%z0mg_col         (:)   = spval
    allocate(this%z0qg_col         (begc:endc)) ; this%z0qg_col         (:)   = spval
    allocate(this%z0hg_col         (begc:endc)) ; this%z0hg_col         (:)   = spval
    allocate(this%z0m_actual_patch (begp:endp)) ; this%z0m_actual_patch (:)   = spval

  end subroutine InitAllocate
  !-----------------------------------------------------------------------


  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    this%z0mg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0MG', units='m', &
         avgflag='A', long_name='roughness length over ground, momentum', &
         ptr_col=this%z0mg_col, default='inactive')

    this%z0hg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0HG', units='m', &
         avgflag='A', long_name='roughness length over ground, sensible heat', &
         ptr_col=this%z0hg_col, default='inactive')

    this%z0qg_col(begc:endc) = spval
    call hist_addfld1d (fname='Z0QG', units='m', &
         avgflag='A', long_name='roughness length over ground, latent heat', &
         ptr_col=this%z0qg_col, default='inactive')

    this%va_patch(begp:endp) = spval
    call hist_addfld1d (fname='VA', units='m/s', &
         avgflag='A', long_name='atmospheric wind speed plus convective velocity', &
         ptr_patch=this%va_patch, default='inactive')

    this%u10_elm_patch(begp:endp) = spval
    call hist_addfld1d (fname='U10', units='m/s', &
         avgflag='A', long_name='10-m wind', &
         ptr_patch=this%u10_elm_patch)

    if (use_cn) then
       this%u10_patch(begp:endp) = spval
       call hist_addfld1d (fname='U10_DUST', units='m/s', &
            avgflag='A', long_name='10-m wind for dust model', &
            ptr_patch=this%u10_patch, default='inactive')
    end if

    if (use_cn) then
       this%ram1_patch(begp:endp) = spval
       call hist_addfld1d (fname='RAM1', units='s/m', &
            avgflag='A', long_name='aerodynamical resistance ', &
            ptr_patch=this%ram1_patch, default='inactive')
    end if

    if (use_cn) then
       this%fv_patch(begp:endp) = spval
       call hist_addfld1d (fname='FV', units='m/s', &
            avgflag='A', long_name='friction velocity for dust model', &
            ptr_patch=this%fv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0hv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0HV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, sensible heat', &
            ptr_patch=this%z0hv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0m_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0M', units='m', &
            avgflag='A', long_name='momentum roughness length', &
            ptr_patch=this%z0m_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0mv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0MV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, momentum', &
            ptr_patch=this%z0mv_patch, default='inactive')
    end if

    if (use_cn) then
       this%z0qv_patch(begp:endp) = spval
       call hist_addfld1d (fname='Z0QV', units='m', &
            avgflag='A', long_name='roughness length over vegetation, latent heat', &
            ptr_patch=this%z0qv_patch, default='inactive')
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l                         ! indices
    !-----------------------------------------------------------------------

    ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
    ! since this is not initialized before first call to VegStructUpdate,
    ! and it is required to set the upper bound for canopy top height.
    ! Changed 3/21/08, KO: still needed but don't have sufficient information
    ! to set this properly (e.g., patch-level displacement height and roughness
    ! length). So leave at 30m.

    if (use_cn .or. use_fates) then
       do p = bounds%begp, bounds%endp
          this%forc_hgt_u_patch(p) = 30._r8
       end do
    end if

    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%lakpoi(l)) then !lake
          this%z0mg_col(c) = 0.0004_r8
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use spmdMod    , only : masterproc
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(frictionvel_type) :: this
    type(bounds_type) , intent(in)    :: bounds
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='Z0MG', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground momentum roughness length', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%z0mg_col)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine ReadNamelist( this, NLFilename )
    !
    ! !DESCRIPTION:
    ! Read the namelist for Friction Velocity
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use elm_varctl     , only : iulog
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use abortutils     , only : endrun
    !
    ! !ARGUMENTS:
    class(frictionvel_type), intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'FrictionVelocityReadNamelist'
    character(len=*), parameter :: nmlname = 'friction_velocity'
    !-----------------------------------------------------------------------
    real(r8) :: zetamaxstable
    namelist /friction_velocity/ zetamaxstable

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    zetamaxstable = 0.5_r8

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=friction_velocity, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(__FILE__, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (zetamaxstable, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=friction_velocity)
       write(iulog,*) ' '
    end if

    this%zetamaxstable = zetamaxstable

  end subroutine ReadNamelist

  !-----------------------------------------------------------------------
  subroutine SetActualRoughnessLengths(this, bounds, &
       num_nolakeurbanp, filter_nolakeurbanp, &
       num_urbanp, filter_urbanp, &
       num_lakep, filter_lakep, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! Set roughness lengths actually used in flux calculations
    !
    ! !ARGUMENTS:
    class(frictionvel_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_nolakeurbanp        ! number of points in filter_exposedvegp
    integer                 , intent(in)    :: filter_nolakeurbanp(:)  ! patch filter for non-snow-covered veg
    integer                 , intent(in)    :: num_urbanp              ! number of points in filter_urbanp
    integer                 , intent(in)    :: filter_urbanp(:)        ! patch filter for urban
    integer                 , intent(in)    :: num_lakep               ! number of points in filter_lakep
    integer                 , intent(in)    :: filter_lakep(:)         ! patch filter for lake
    class(canopystate_type) , intent(in)    :: canopystate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p, c, l

    character(len=*), parameter :: subname = 'SetActualRoughnessLengths'
    !-----------------------------------------------------------------------

    associate( &
         z_0_town   => lun_pp%z_0_town          , & ! Input:  [real(r8) (:)] momentum roughness length of urban landunit [m]

         z0mv       => this%z0mv_patch       , & ! Input:  [real(r8) (:)] roughness length over vegetation, momentum [m]
         z0mg       => this%z0mg_col         , & ! Input:  [real(r8) (:)] roughness length over ground, momentum [m]
         z0m_actual => this%z0m_actual_patch , & ! Output: [real(r8) (:)] roughness length actually used in flux calculations, momentum [m]
         frac_veg_nosno  => canopystate_vars%frac_veg_nosno_patch  & ! fraction of vegetation not covered by snow [patch]
         )

    do fp = 1, num_nolakeurbanp 
       p = filter_nolakeurbanp(fp)
       c = veg_pp%column(p)

       if (frac_veg_nosno(p) > 0._r8) then
          z0m_actual(p) = z0mv(p)
       else
          z0m_actual(p) = z0mg(c)
       end if
    end do

    do fp = 1, num_urbanp
       p = filter_urbanp(fp)
       l = veg_pp%landunit(p)

       z0m_actual(p) = z_0_town(l)
    end do

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = veg_pp%column(p)

       z0m_actual(p) = z0mg(c)
    end do

    end associate
  end subroutine SetActualRoughnessLengths

  subroutine MoninObukIni (this, ur, thv, dthv, zldis, z0m, um, obu)
    !$acc routine seq
    ! !DESCRIPTION:
    ! Initialization of the Monin-Obukhov length.
    ! The scheme is based on the work of Zeng et al. (1998):
    ! Intercomparison of bulk aerodynamic algorithms for the computation
    ! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
    ! Vol. 11, 2628-2644.
    !
    ! !USES:
    use elm_varcon, only : grav
    !
    ! !ARGUMENTS:
    implicit none
    class(frictionvel_type), intent(in) :: this
    real(r8), intent(in)  :: ur    ! wind speed at reference height [m/s]
    real(r8), intent(in)  :: thv   ! virtual potential temperature (kelvin)
    real(r8), intent(in)  :: dthv  ! diff of vir. poten. temp. between ref. height and surface
    real(r8), intent(in)  :: zldis ! reference height "minus" zero displacement heght [m]
    real(r8), intent(in)  :: z0m   ! roughness length, momentum [m]
    real(r8), intent(out) :: um    ! wind speed including the stability effect [m/s]
    real(r8), intent(out) :: obu   ! monin-obukhov length (m)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wc    ! convective velocity [m/s]
    real(r8) :: rib   ! bulk Richardson number
    real(r8) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: ustar ! friction velocity [m/s]
    !-----------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar=0.06_r8
    wc=0.5_r8
    if (dthv >= 0._r8) then
       um=max(ur,0.1_r8)
    else
       um=sqrt(ur*ur+wc*wc)
    endif

    rib=grav*zldis*dthv/(thv*um*um)

    if (rib >= 0._r8) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1._r8-5._r8*min(rib,0.19_r8))
!hh!       zeta = min(2._r8,max(zeta,0.01_r8 ))
       zeta = min(this%zetamaxstable,max(zeta,0.01_r8 ))
    else                     ! unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta

  end subroutine MoninObukIni

end module FrictionVelocityType
