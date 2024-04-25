module elm_initializeMod

  !-----------------------------------------------------------------------
  ! Performs land model initialization
  !-----------------------------------------------------------------------

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_sys_mod           , only : shr_sys_flush
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use spmdMod               , only : masterproc, mpicom
  use decompMod             , only : bounds_type, get_proc_bounds, get_proc_clumps, get_clump_bounds
  use abortutils            , only : endrun
  use elm_varctl            , only : nsrest, nsrStartup, nsrContinue, nsrBranch, use_fates_sp

  use elm_varctl            , only : create_glacier_mec_landunit, iulog
  use elm_varctl            , only : use_lch4, use_cn, use_cndv, use_voc, use_c13, use_c14, use_fates

  use elm_varsur            , only : wt_lunit, urban_valid, wt_nat_patch, wt_cft, wt_glc_mec, topo_glc_mec,firrig,f_surf,f_grd 
  use elm_varsur            , only : fert_cft
  use elm_varsur            , only : wt_tunit, elv_tunit, slp_tunit,asp_tunit,num_tunit_per_grd
  use perf_mod              , only : t_startf, t_stopf
  !use readParamsMod         , only : readParameters
  use readParamsMod         , only : readSharedParameters, readPrivateParameters
  use ncdio_pio             , only : file_desc_t
  use reweightMod           , only : reweight_wrapup
  use filterMod             , only : allocFilters, filter, filter_inactive_and_active
  use ELMFatesInterfaceMod  , only : ELMFatesGlobals1,ELMFatesGlobals2
  use ELMFatesInterfaceMod  , only : ELMFatesTimesteps
  use dynSubgridControlMod  , only : dynSubgridControl_init
  use CLMFatesParamInterfaceMod, only: FatesReadPFTs
!hh!  use BeTRSimulationELM, only : create_betr_simulation_elm
  !
  !-----------------------------------------
  ! Definition of component types
  !-----------------------------------------
  use GridcellType           , only : grc_pp
  use TopounitType           , only : top_pp
  use TopounitDataType       , only : top_as, top_af, top_es
  use LandunitType           , only : lun_pp
  use ColumnType             , only : col_pp
  use ColumnDataType         , only : col_es
  use VegetationType         , only : veg_pp
  use VegetationDataType     , only : veg_es

  use elm_instMod
  use WaterBudgetMod         , only : WaterBudget_Reset
  use CNPBudgetMod           , only : CNPBudget_Reset
  use elm_varctl             , only : do_budgets
  !
  implicit none
  save
  !
  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
!hh!  subroutine initialize1( )
  subroutine initialize1(dtime)
    !
    ! !DESCRIPTION:
    ! CLM initialization first phase
    !
    ! !USES:
    use elm_varpar           , only: elm_varpar_init
    use elm_varcon           , only: elm_varcon_init
    use landunit_varcon      , only: landunit_varcon_init
    use elm_varctl           , only: fsurdat, version
    use elm_varpar           , only: update_pft_array_bounds
    use controlMod           , only: control_init, control_print, NLFilename
    use ncdio_pio            , only: ncd_pio_init
    use initGridCellsMod     , only: initGridCells
    use SoilTemperatureMod   , only: init_soil_temperature
    use dynSubgridControlMod , only: dynSubgridControl_init
    !
    ! !ARGUMENTS
    integer, intent(in) :: dtime    ! model time step (seconds)
    !
    ! !LOCAL VARIABLES:
    integer           :: ier                     ! error status
    integer           :: i,j,n,k,c,l,g           ! indices
    integer           :: nl                      ! gdc and glo lnd indices
    integer           :: ns, ni, nj              ! global grid sizes
    integer           :: begg, endg              ! processor bounds
    type(bounds_type) :: bounds_proc
    type(bounds_type) :: bounds_clump            ! clump bounds
    integer           :: nclumps                 ! number of clumps on this processor
    integer           :: nc                      ! clump index
    integer ,pointer  :: amask(:)                ! global land mask
    character(len=32) :: subname = 'initialize1' ! subroutine name
    !-----------------------------------------------------------------------

    call t_startf('elm_init1')

    ! ------------------------------------------------------------------------
    ! Initialize run control variables, timestep
    ! ------------------------------------------------------------------------

    if ( masterproc )then
       write(iulog,*) trim(version)
       write(iulog,*)
       write(iulog,*) 'Attempting to initialize the land model .....'
       write(iulog,*)
       call shr_sys_flush(iulog)
    endif

!hh!    call control_init()
    call control_init(dtime)
    call elm_varpar_init()
    call elm_varcon_init()
    call landunit_varcon_init()
    call ncd_pio_init()
    if(use_fates) then
       ! Allow FATES to dictate the number of patches per column.
       ! We still use numcft as dictated by
       ! the host model.
       ! This call will override natpft_size (and its bounds
       ! in the following call) for FATES runs
       call ELMFatesGlobals1()
       call update_pft_array_bounds()
    end if    
    
    call elm_petsc_init()
    call init_soil_temperature()

    if (masterproc) call control_print()

    call dynSubgridControl_init(NLFilename)

    call t_stopf('elm_init1')

  end subroutine initialize1

  !-----------------------------------------------------------------------
!hh!  subroutine initialize2( )
  subroutine initialize2(ni,nj)
    !
    ! !DESCRIPTION:
    ! CLM initialization - second phase
    !
    ! !USES:
    use elm_varcon                    , only : spval
    use elm_varpar                    , only : natpft_lb, natpft_ub, cft_lb, cft_ub, maxpatch_glcmec
    use elm_varpar                    , only : surfpft_lb, surfpft_ub
    use elm_varpar                    , only : nlevsno, numpft, crop_prog, nlevsoi,max_patch_per_col

    use elm_varctl                    , only : fsurdat, fatmlndfrc, flndtopo, fglcmask, noland, version
    use elm_varctl                    , only : finidat, finidat_interp_source, finidat_interp_dest, fsurdat
    use elm_varctl                    , only : use_cn, use_fates

    use elm_varorb                    , only : eccen, mvelpp, lambm0, obliqr
    use landunit_varcon               , only : landunit_varcon_init, max_lunit, numurbl
    use pftvarcon                     , only : pftconrd
    use decompInitMod                 , only : decompInit_clumps, decompInit_gtlcp
    use domainMod                     , only : domain_check, ldomain, domain_init
    use surfrdMod                     , only : surfrd_get_data, surfrd_get_globmask, surfrd_get_grid, surfrd_get_topo, surfrd_get_topo_for_solar_rad
    use controlMod                    , only : NLFilename
    use initGridCellsMod              , only : initGridCells
    use CH4varcon                     , only : CH4conrd
    use UrbanParamsType               , only : UrbanInput
    use shr_orb_mod                   , only : shr_orb_decl
    use seq_drydep_mod                , only : n_drydep
    use accumulMod                    , only : print_accum_fields
    use clm_time_manager              , only : get_step_size_real, get_curr_calday
    use clm_time_manager              , only : get_curr_date, get_nstep, advance_timestep
    use clm_time_manager              , only : timemgr_init, timemgr_restart_io, timemgr_restart

    use DaylengthMod                  , only : InitDaylength
    use dynSubgridDriverMod           , only : dynSubgrid_init
    use fileutils                     , only : getfil
    use initInterpMod                 , only : initInterp
    use subgridWeightsMod             , only : init_subgrid_weights_mod
    use histFileMod                   , only : hist_htapes_build, htapes_fieldlist
    use histFileMod                   , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    use restFileMod                   , only : restFile_getfile, restFile_open, restFile_close
    use restFileMod                   , only : restFile_read, restFile_write
    use ndepStreamMod                 , only : ndep_init, ndep_interp
    use EcosystemDynMod               , only : EcosystemDynInit
    use pdepStreamMod                 , only : pdep_init, pdep_interp
    use DecompCascadeBGCMod           , only : init_decompcascade_bgc
    use DecompCascadeCNMod            , only : init_decompcascade_cn
    use CNDecompCascadeContype        , only : init_decomp_cascade_constants
    use VegetationPropertiesType      , only : veg_vp
    use SoilorderConType              , only : soilorderconInit
    use LakeCon                       , only : LakeConInit
    use SatellitePhenologyMod         , only : SatellitePhenologyInit, readAnnualVegetation, interpMonthlyVeg, SatellitePhenology
    use SnowSnicarMod                 , only : SnowAge_init, SnowOptics_init
    use lnd2atmMod                    , only : lnd2atm_minimal
    use controlMod                    , only : NLFilename
    use elm_varctl                    , only : use_elm_interface, use_pflotran
    use elm_interface_pflotranMod     , only : elm_pf_interface_init !, elm_pf_set_restart_stamp
    use clm_time_manager              , only : is_restart
    use topounit_varcon               , only : max_topounits, has_topounit, topounit_varcon_init    

    use elm_varctl                    , only : use_top_solar_rad
    use elm_varctl                    , only : use_century_decomp
    use elm_varcon                    , only : h2osno_max, bdsno
    use column_varcon                 , only : col_itype_to_icemec_class
    use seq_drydep_mod                , only : drydep_method, DD_XLND
    use landunit_varcon               , only : istice_mec
    use soilorder_varcon              , only : soilorder_conrd
    use surfrdMod                     , only : surfrd_get_grid_conn, surfrd_topounit_data
    use elm_varctl                    , only : lateral_connectivity, domain_decomp_type
    use decompInitMod                 , only : decompInit_lnd_using_gp, decompInit_ghosts
    use domainLateralMod              , only : ldomain_lateral, domainlateral_init
    !
    ! !ARGUMENTS
    implicit none
    !
    ! !ARGUMENTS
    integer, intent(in) :: ni, nj                ! global grid sizes
    !
    ! !LOCAL VARIABLES:
    integer            :: c,g,i,j,k,l,n,p,t,ti,topi ! indices
    integer            :: yr              ! current year (0, ...)
    integer            :: mon             ! current month (1 -> 12)
    integer            :: day             ! current day (1 -> 31)
    integer            :: ncsec           ! current time of day [seconds]
    integer            :: nc              ! clump index
    integer            :: ns              ! indices
    integer            :: nclumps         ! number of clumps on this processor
    integer            :: icemec_class    ! current icemec class (1..maxpatch_glcmec)
    character(len=256) :: fnamer          ! name of netcdf restart file
    character(len=256) :: pnamer          ! full pathname of netcdf restart file
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    real(r8)           :: dtime           ! time step increment (sec)
    integer            :: nstep           ! model time step
    real(r8)           :: calday          ! calendar day for nstep
    real(r8)           :: caldaym1        ! calendar day for nstep-1
    real(r8)           :: declin          ! solar declination angle in radians for nstep
    real(r8)           :: declinm1        ! solar declination angle in radians for nstep-1
    real(r8)           :: eccf            ! earth orbit eccentricity factor
    type(bounds_type)  :: bounds_proc     ! processor bounds
    type(bounds_type)  :: bounds_clump    ! clump bounds
    integer ,pointer   :: amask(:)                ! global land mask
    integer ,pointer   :: cellsOnCell(:,:)        ! grid cell level connectivity
    integer ,pointer   :: edgesOnCell(:,:)        ! index to determine distance between neighbors from dcEdge
    integer ,pointer   :: nEdgesOnCell(:)         ! number of edges
    real(r8), pointer  :: dcEdge(:)               ! distance between centroids of grid cells
    real(r8), pointer  :: dvEdge(:)               ! distance between vertices
    real(r8), pointer  :: areaCell(:)             ! area of grid cells [m^2]
    integer            :: nCells_loc              ! number of grid cell level connectivity saved locally
    integer            :: nEdges_loc              ! number of edge length saved locally
    integer            :: maxEdges                ! max number of edges/neighbors
    integer            :: begg, endg
    integer            :: iun
    integer            :: klen
    integer            :: ioe
    integer            :: ier
    logical            :: lexists
    real(r8), pointer  :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    character(len=32)  :: subname = 'initialize2'
    !----------------------------------------------------------------------

    call t_startf('elm_init2')

!hh!    ! ------------------------------------------------------------------------
!hh!    ! Read in global land grid and land mask (amask)- needed to set decomposition
!hh!    ! ------------------------------------------------------------------------
!hh!
!hh!    ! global memory for amask is allocate in surfrd_get_glomask - must be
!hh!    ! deallocated below
!hh!    if (masterproc) then
!hh!       write(iulog,*) 'Attempting to read global land mask from ',trim(fatmlndfrc)
!hh!       call shr_sys_flush(iulog)
!hh!    endif
!hh!    call surfrd_get_globmask(filename=fatmlndfrc, mask=amask, ni=ni, nj=nj)
!hh!
!hh!    ! Exit early if no valid land points
!hh!    if ( all(amask == 0) )then
!hh!       if (masterproc) write(iulog,*) trim(subname)//': no valid land points do NOT run elm'
!hh!       noland = .true.
!hh!       return
!hh!    end if
!hh!
!hh!
!hh!    ! ------------------------------------------------------------------------
!hh!    ! If specified, read the grid level connectivity
!hh!    ! ------------------------------------------------------------------------
!hh!
!hh!    if (lateral_connectivity) then
!hh!       call surfrd_get_grid_conn(fatmlndfrc, cellsOnCell, edgesOnCell, &
!hh!            nEdgesOnCell, areaCell, dcEdge, dvEdge, &
!hh!            nCells_loc, nEdges_loc, maxEdges)
!hh!    else
!hh!       nullify(cellsOnCell)
!hh!       nCells_loc = 0
!hh!       maxEdges   = 0
!hh!    endif
!hh!
!hh!    ! ------------------------------------------------------------------------
!hh!    ! Determine clm gridcell decomposition and processor bounds for gridcells
!hh!    ! ------------------------------------------------------------------------
!hh!
!hh!    select case (trim(domain_decomp_type))
!hh!    case ("round_robin")
!hh!       call decompInit_lnd(ni, nj, amask)
!hh!       deallocate(amask)
!hh!    case ("graph_partitioning")
!hh!       call decompInit_lnd_using_gp(ni, nj, cellsOnCell, nCells_loc, maxEdges, amask)
!hh!    case default
!hh!       call endrun(msg='ERROR elm_initializeMod: '//&
!hh!            'Unsupported domain_decomp_type = ' // trim(domain_decomp_type))
!hh!    end select
!hh!
!hh!    if (lateral_connectivity) then
!hh!       call domainlateral_init(ldomain_lateral, cellsOnCell, edgesOnCell, &
!hh!            nEdgesOnCell, areaCell, dcEdge, dvEdge, &
!hh!            nCells_loc, nEdges_loc, maxEdges)
!hh!    endif

    ! *** Get JUST gridcell processor bounds ***
    ! Remaining bounds (landunits, columns, patches) will be determined
    ! after the call to decompInit_glcp - so get_proc_bounds is called
    ! twice and the gridcell information is just filled in twice

!hh!    call get_proc_bounds(begg, endg)
    ! Get processor bounds for gridcells
    call get_proc_bounds(bounds_proc)
    begg = bounds_proc%begg; endg = bounds_proc%endg


!hh!    ! ------------------------------------------------------------------------
!hh!    ! Get grid and land fraction (set ldomain)
!hh!    ! ------------------------------------------------------------------------
!hh!
!hh!    if (masterproc) then
!hh!       write(iulog,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
!hh!       call shr_sys_flush(iulog)
!hh!    endif
!hh!    if (create_glacier_mec_landunit) then
!hh!       call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc, fglcmask)
!hh!    else
!hh!       call surfrd_get_grid(begg, endg, ldomain, fatmlndfrc)
!hh!    endif
!hh!    if (masterproc) then
!hh!       call domain_check(ldomain)
!hh!    endif
!hh!    ldomain%mask = 1  !!! TODO - is this needed?
!hh!
!hh!    ! Get topo if appropriate (set ldomain%topo)
!hh!
!hh!    if (flndtopo /= " ") then
!hh!       if (masterproc) then
!hh!          write(iulog,*) 'Attempting to read atm topo from ',trim(flndtopo)
!hh!          call shr_sys_flush(iulog)
!hh!       endif
!hh!
!hh!       call surfrd_get_topo(ldomain, flndtopo)  
!hh!    endif    
!hh!    
!hh!    if (fsurdat /= " " .and. use_top_solar_rad) then
!hh!       if (masterproc) then
!hh!          write(iulog,*) 'Attempting to read topo parameters for TOP solar radiation parameterization from ',trim(fsurdat)
!hh!          call shr_sys_flush(iulog)
!hh!       endif
!hh!       call surfrd_get_topo_for_solar_rad(ldomain, fsurdat)  
!hh!
!hh!    endif
    
    !-------------------------------------------------------------------------
    ! Topounit
    !-------------------------------------------------------------------------
    call topounit_varcon_init(begg, endg,fsurdat,ldomain)  ! Topounits
    !-------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------
    ! Initialize urban model input (initialize urbinp data structure)
    ! This needs to be called BEFORE the call to surfrd_get_data since
    ! that will call surfrd_get_special which in turn calls check_urban

    call UrbanInput(begg, endg, mode='initialize')

    ! Allocate surface grid dynamic memory (just gridcell bounds dependent)

    allocate (wt_lunit     (begg:endg,1:max_topounits, max_lunit           )) 
    allocate (urban_valid  (begg:endg,1:max_topounits                      ))
    allocate (wt_nat_patch (begg:endg,1:max_topounits, surfpft_lb:surfpft_ub ))
    allocate (wt_cft       (begg:endg,1:max_topounits, cft_lb:cft_ub       ))
    allocate (fert_cft     (begg:endg,1:max_topounits, cft_lb:cft_ub       ))
    if (create_glacier_mec_landunit) then
       allocate (wt_glc_mec  (begg:endg,1:max_topounits, maxpatch_glcmec))
       allocate (topo_glc_mec(begg:endg,1:max_topounits, maxpatch_glcmec))
    else
       allocate (wt_glc_mec  (1,1,1))
       allocate (topo_glc_mec(1,1,1))
    endif
    
    allocate (wt_tunit  (begg:endg,1:max_topounits  )) 
    allocate (elv_tunit (begg:endg,1:max_topounits  ))
    allocate (slp_tunit (begg:endg,1:max_topounits  ))
    allocate (asp_tunit (begg:endg,1:max_topounits  ))
    allocate (num_tunit_per_grd (begg:endg))
    allocate (firrig  (begg:endg,1:max_topounits  ))
    allocate (f_surf  (begg:endg,1:max_topounits  ))
    allocate (f_grd  (begg:endg,1:max_topounits  ))

    ! Read list of Patches and their corresponding parameter values
    ! Independent of model resolution, Needs to stay before surfrd_get_data

    call pftconrd()
    call soilorder_conrd()

    ! Read in FATES parameter values early in the call sequence as well
    ! The PFT file, specifically, will dictate how many pfts are used
    ! in fates, and this will influence the amount of memory we
    ! request from the model, which is relevant in set_fates_global_elements()
    if (use_fates) then
       call FatesReadPFTs()
    end if
    
    ! Read surface dataset and set up subgrid weight arrays
    call surfrd_get_data(begg, endg, ldomain, fsurdat)

    if(use_fates) then

       ! Pass various control flags to FATES and setup
       ! FATES allocations
       ! --------------------------------------------------------------------

       call ELMFatesGlobals2()

    end if

    
    ! ------------------------------------------------------------------------
    ! Determine decomposition of subgrid scale topounits, landunits, topounits, columns, patches
    ! ------------------------------------------------------------------------

    if (create_glacier_mec_landunit) then
!hh!       call decompInit_clumps(ldomain%glcmask)
!hh!       call decompInit_ghosts(ldomain%glcmask)
       call decompInit_clumps(ni, nj, ldomain%glcmask)
       call decompInit_ghosts(ni, nj, ldomain%glcmask)
    else
       call decompInit_clumps(ni, nj)
       call decompInit_ghosts(ni, nj)
    endif

    ! *** Get ALL processor bounds - for gridcells, landunit, columns and patches ***

    call get_proc_bounds(bounds_proc)

    ! Allocate memory for subgrid data structures
    ! This is needed here BEFORE the following call to initGridcells
    ! Note that the assumption is made that none of the subgrid initialization
    ! can depend on other elements of the subgrid in the calls below

    ! Initialize the gridcell data types
    call grc_pp%Init (bounds_proc%begg, bounds_proc%endg)
    
    ! Read topounit information from fsurdat
    if (has_topounit) then
         call surfrd_topounit_data(begg, endg, fsurdat)         
    end if
    
    ! Initialize the topographic unit data types
    call top_pp%Init (bounds_proc%begt, bounds_proc%endt) ! topology and physical properties
    call top_as%Init (bounds_proc%begt, bounds_proc%endt) ! atmospheric state variables (forcings)
    call top_af%Init (bounds_proc%begt, bounds_proc%endt) ! atmospheric flux variables (forcings)
    call top_es%Init (bounds_proc%begt, bounds_proc%endt) ! energy state

    ! Initialize the landunit data types
    call lun_pp%Init (bounds_proc%begl, bounds_proc%endl)

    ! Initialize the column data types
    call col_pp%Init (bounds_proc%begc, bounds_proc%endc)

    ! Initialize the vegetation (PFT) data types
    call veg_pp%Init (bounds_proc%begp, bounds_proc%endp)

    ! Initialize the cohort data types (nothing here yet)
    ! ...to be added later...

!hh!    ! Determine the number of active external models.
!hh!    call EMI_Determine_Active_EMs()

    ! Build hierarchy and topological info for derived types
    ! This is needed here for the following call to decompInit_glcp

    nclumps = get_proc_clumps()
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call initGridCells(bounds_clump)
    end do
    !$OMP END PARALLEL DO

    ! Set global seg maps for gridcells, topounits, landlunits, columns and patches
    !if(max_topounits > 1) then 
    !   if (create_glacier_mec_landunit) then
    !      call decompInit_gtlcp(ns, ni, nj, ldomain%glcmask,ldomain%num_tunits_per_grd)
    !   else
    !      call decompInit_gtlcp(ns, ni, nj,ldomain%num_tunits_per_grd)
    !   endif
    !else
    if (create_glacier_mec_landunit) then
       call decompInit_gtlcp(ni, nj, ldomain%glcmask)
    else
       call decompInit_gtlcp(ni, nj)
    endif
    !endif

    ! Set filters

    call t_startf('init_filters')
    call allocFilters()
    call t_stopf('init_filters')
    
!hh!    nclumps = get_proc_clumps()
    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call reweight_wrapup(bounds_clump, &
            ldomain%glcmask(bounds_clump%begg:bounds_clump%endg)*1._r8)
    end do
    !$OMP END PARALLEL DO

    ! ------------------------------------------------------------------------
    ! Remainder of initialization1
    ! ------------------------------------------------------------------------

    ! Set CH4 Model Parameters from namelist.
    ! Need to do before initTimeConst so that it knows whether to
    ! look for several optional parameters on surfdata file.

    if (use_lch4) then
       call CH4conrd()
    end if

    ! Deallocate surface grid dynamic memory for variables that aren't needed elsewhere.
    ! Some things are kept until the end of initialize2; urban_valid is kept through the
    ! end of the run for error checking.

    !deallocate (wt_lunit, wt_cft, wt_glc_mec)
    deallocate (wt_cft, wt_glc_mec)    !wt_lunit not deallocated because it is being used in CanopyHydrologyMod.F90
    deallocate (wt_tunit, elv_tunit, slp_tunit, asp_tunit,num_tunit_per_grd)
    call t_stopf('elm_init1')

    ! initialize glc_topo
    ! TODO - does this belong here?
    do c = bounds_proc%begc, bounds_proc%endc
       l = col_pp%landunit(c)
       g = col_pp%gridcell(c)
       t = col_pp%topounit(c)
       topi = grc_pp%topi(g)
       ti = t - topi + 1

       if (lun_pp%itype(l) == istice_mec) then
          ! For ice_mec landunits, initialize glc_topo based on surface dataset; this
          ! will get overwritten in the run loop by values sent from CISM
          icemec_class = col_itype_to_icemec_class(col_pp%itype(c))
          col_pp%glc_topo(c) = topo_glc_mec(g,ti, icemec_class)
       else
          ! For other landunits, arbitrarily initialize glc_topo to 0 m; for landunits
          ! where this matters, this will get overwritten in the run loop by values sent
          ! from CISM
          col_pp%glc_topo(c) = 0._r8
       end if
    end do

    if (do_budgets) then
       call WaterBudget_Reset('all')
       if (use_cn) then
          call CNPBudget_Reset('all')
       endif
    endif

    ! ------------------------------------------------------------------------
    ! Determine processor bounds and clumps for this processor
    ! ------------------------------------------------------------------------

    call get_proc_bounds(bounds_proc)
    nclumps = get_proc_clumps()

    ! ------------------------------------------------------------------------
    ! Read in shared parameters files
    ! ------------------------------------------------------------------------

    call readSharedParameters()

    ! ------------------------------------------------------------------------
    ! Initialize time manager
    ! ------------------------------------------------------------------------
    if (nsrest == nsrStartup) then  
       call timemgr_init()
    else
       call restFile_getfile(file=fnamer, path=pnamer)
       call restFile_open( flag='read', file=fnamer, ncid=ncid )
       call timemgr_restart_io( ncid=ncid, flag='read' )
       call restFile_close( ncid=ncid )
       call timemgr_restart()
    end if
    
    ! ------------------------------------------------------------------------
    ! Pass model timestep info to FATES
    ! ------------------------------------------------------------------------
    if(use_fates) then
       call ELMFatesTimesteps()
    end if
    
    ! ------------------------------------------------------------------------
    ! Initialize daylength from the previous time step (needed so prev_dayl can be set correctly)
    ! ------------------------------------------------------------------------

    call t_startf('init_orbd')
    calday = get_curr_calday(reuse_day_365_for_day_366=.true.)
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf )
    dtime = get_step_size_real()
    caldaym1 = get_curr_calday(offset=-int(dtime), reuse_day_365_for_day_366=.true.)
    call shr_orb_decl( caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf )
    call t_stopf('init_orbd')
    call InitDaylength(bounds_proc, declin=declin, declinm1=declinm1, obliquity=obliqr)



    ! History file variables
    if (use_cn) then
       call hist_addfld1d (fname='DAYL',  units='s', &
            avgflag='A', long_name='daylength', &
            ptr_gcell=grc_pp%dayl, default='inactive')

       call hist_addfld1d (fname='PREV_DAYL', units='s', &
            avgflag='A', long_name='daylength from previous timestep', &
            ptr_gcell=grc_pp%prev_dayl, default='inactive')
    end if

    ! Initialize component data structures
    ! Note: new logic is in place that sets all the history fields to spval so
    ! there is no guesswork in the initialization to nans of the allocated variables
    ! First put in history calls for subgrid data structures - these cannot appear in the
    ! module for the subgrid data definition due to circular dependencies that are introduced

    data2dptr => col_pp%dz(:,-nlevsno+1:0)
    col_pp%dz(bounds_proc%begc:bounds_proc%endc,:) = spval
    call hist_addfld2d (fname='SNO_Z', units='m', type2d='levsno',  &
         avgflag='A', long_name='Snow layer thicknesses', &
         ptr_col=data2dptr, no_snow_behavior=no_snow_normal, default='inactive')

    col_pp%zii(bounds_proc%begc:bounds_proc%endc) = spval
    call hist_addfld1d (fname='ZII', units='m', &
         avgflag='A', long_name='convective boundary height', &
         ptr_col=col_pp%zii, default='inactive')

    call elm_inst_biogeophys(bounds_proc)


    call SnowOptics_init( ) ! SNICAR optical parameters:

    call SnowAge_init( )    ! SNICAR aging   parameters:

    ! ------------------------------------------------------------------------
    ! Read in private parameters files, this should be preferred for mulitphysics
    ! implementation, jinyun Tang, Feb. 11, 2015
    ! ------------------------------------------------------------------------
    if(use_cn .or. use_fates) then
       call init_decomp_cascade_constants()
    endif
    !read bgc implementation specific parameters when needed
    call readPrivateParameters()

    if (use_cn .or. use_fates) then
!hh!       if (.not. is_active_betr_bgc)then
          if (use_century_decomp) then
           ! Note that init_decompcascade_bgc needs cnstate_vars to be initialized
             call init_decompcascade_bgc(bounds_proc, cnstate_vars, soilstate_vars)
          else
           ! Note that init_decompcascade_cn needs cnstate_vars to be initialized
             call init_decompcascade_cn(bounds_proc, cnstate_vars)
          end if
!hh!       endif
    endif

    ! FATES is instantiated in the following call.  The global is in clm_inst
    call elm_inst_biogeochem(bounds_proc)

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! The time manager needs to be initialized before thes called is made, since
    ! the step size is needed.

    call t_startf('init_accflds')

    call atm2lnd_vars%initAccBuffer(bounds_proc)

    call top_as%InitAccBuffer(bounds_proc)

    call top_af%InitAccBuffer(bounds_proc)

    call veg_es%InitAccBuffer(bounds_proc)

    call canopystate_vars%initAccBuffer(bounds_proc)

    if (crop_prog) then
       call crop_vars%initAccBuffer(bounds_proc)
    end if

    call cnstate_vars%initAccBuffer(bounds_proc)

    call print_accum_fields()

    call t_stopf('init_accflds')

    ! ------------------------------------------------------------------------
    ! Initializate dynamic subgrid weights (for prescribed transient Patches,
    ! and/or dynamic landunits); note that these will be overwritten in a
    ! restart run
    ! ------------------------------------------------------------------------

    call t_startf('init_dyn_subgrid')
    call init_subgrid_weights_mod(bounds_proc)
    call dynSubgrid_init(bounds_proc, glc2lnd_vars, crop_vars)
    call t_stopf('init_dyn_subgrid')

    ! ------------------------------------------------------------------------
    ! Initialize modules (after time-manager initialization in most cases)
    ! ------------------------------------------------------------------------

    if (use_cn .or. use_fates) then
       call EcosystemDynInit(bounds_proc,alm_fates)
    else
       call SatellitePhenologyInit(bounds_proc)
    end if

    if (use_fates_sp) then
       call SatellitePhenologyInit(bounds_proc)
    end if

    if (use_cn .and. n_drydep > 0 .and. drydep_method == DD_XLND) then
       ! Must do this also when drydeposition is used so that estimates of monthly
       ! differences in LAI can be computed
       call SatellitePhenologyInit(bounds_proc)
    end if

    ! ------------------------------------------------------------------------
    ! On restart only - process the history namelist.
    ! ------------------------------------------------------------------------

    ! Later the namelist from the restart file will be used.  This allows basic
    ! checking to make sure you didn't try to change the history namelist on restart.

    if (nsrest == nsrContinue ) then
       call htapes_fieldlist()
    end if

    ! ------------------------------------------------------------------------
    ! Read restart/initial info
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup) then

       if (finidat == ' ') then
          if (finidat_interp_source == ' ') then
             if (masterproc) then
                write(iulog,*)'Using cold start initial conditions '
             end if
          else
             if (masterproc) then
                write(iulog,*)'Interpolating initial conditions from ',trim(finidat_interp_source),&
                     ' and creating new initial conditions ', trim(finidat_interp_dest)
             end if
          end if
       else
          if (trim(finidat) == trim(finidat_interp_dest)) then
             ! Check to see if status file for finidat exists
             klen = len_trim(finidat_interp_dest) - 3 ! remove the .nc
             locfn = finidat_interp_dest(1:klen)//'.status'
             inquire(file=trim(locfn), exist=lexists)
             if (.not. lexists) then
                if (masterproc) then
                   write(iulog,'(a)')' failed to find file '//trim(locfn)
                   write(iulog,'(a)')' this indicates a problem in creating '//trim(finidat_interp_dest)
                   write(iulog,'(a)')' remove '//trim(finidat_interp_dest)//' and try again'
                end if
                call endrun()
             end if
          end if
          if (masterproc) then
             write(iulog,*)'Reading initial conditions from ',trim(finidat)
          end if
          call getfil( finidat, fnamer, 0 )
          call restFile_read(bounds_proc, fnamer,                           &
               atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,  &
               ch4_vars, energyflux_vars, frictionvel_vars, lakestate_vars, &
               photosyns_vars, soilhydrology_vars,                          &
               soilstate_vars, solarabs_vars, surfalb_vars,                 &
!hh!               sedflux_vars, ep_betr, alm_fates, glc2lnd_vars, crop_vars)
               sedflux_vars, alm_fates, glc2lnd_vars, crop_vars)

         call WaterBudget_Reset('all')
         if (use_cn) then
            call CNPBudget_Reset('all')
         endif

       end if

    else if ((nsrest == nsrContinue) .or. (nsrest == nsrBranch)) then
       if (masterproc) then
          write(iulog,*)'Reading restart file ',trim(fnamer)
       end if
       call restFile_read(bounds_proc, fnamer,                           &
            atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,  &
            ch4_vars, energyflux_vars, frictionvel_vars, lakestate_vars ,&
            photosyns_vars, soilhydrology_vars,                          &
            soilstate_vars, solarabs_vars, surfalb_vars,                 &
!hh!            sedflux_vars, ep_betr, alm_fates, glc2lnd_vars, crop_vars)
            sedflux_vars, alm_fates, glc2lnd_vars, crop_vars)

    end if

    ! ------------------------------------------------------------------------
    ! If appropriate, create interpolated initial conditions
    ! ------------------------------------------------------------------------

    if (nsrest == nsrStartup .and. finidat_interp_source /= ' ') then

       ! Check that finidat is not cold start - abort if it is
       if (finidat /= ' ') then
          call endrun(msg='ERROR elm_initializeMod: '//&
               'finidat and finidat_interp_source cannot both be non-blank')
       end if

       ! Determine name if finidat_interp_dest status file
       klen = len_trim(finidat_interp_dest) - 3 ! remove the .nc
       locfn = finidat_interp_dest(1:klen)//'.status'

       ! Remove file if it already exists
       if (masterproc) then
          inquire(file=trim(locfn), exist=lexists)
          if (lexists) then
             open(unit=9876, file=locfn, status='old', iostat=ioe)
             if (ioe == 0) then
                close(9876, status='delete')
             end if
          end if
       end if
       call mpi_barrier(mpicom,ier)

       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1, nclumps
          call get_clump_bounds(nc, bounds_clump)
          call reweight_wrapup(bounds_clump, &
               glc2lnd_vars%icemask_grc(bounds_clump%begg:bounds_clump%endg))
       end do
       !$OMP END PARALLEL DO
!hh!       if(use_betr)then
!hh!         call ep_betr%set_active(bounds_proc, col_pp)
!hh!       endif
       ! Create new template file using cold start
       call restFile_write(bounds_proc, finidat_interp_dest,             &
            atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,  &
            ch4_vars, energyflux_vars, frictionvel_vars, lakestate_vars, &
            photosyns_vars, soilhydrology_vars,          &
            soilstate_vars, solarabs_vars, surfalb_vars, &
!hh!            sedflux_vars, ep_betr, alm_fates, crop_vars)
            sedflux_vars, alm_fates, crop_vars)

       ! Interpolate finidat onto new template file
       call getfil( finidat_interp_source, fnamer,  0 )
       call initInterp(filei=fnamer, fileo=finidat_interp_dest, bounds=bounds_proc)

       ! Read new interpolated conditions file back in
       call restFile_read(bounds_proc, finidat_interp_dest,              &
            atm2lnd_vars, aerosol_vars, canopystate_vars, cnstate_vars,  &
            ch4_vars, energyflux_vars, frictionvel_vars, lakestate_vars, &
            photosyns_vars, soilhydrology_vars,            &
            soilstate_vars, solarabs_vars, surfalb_vars,   &
!hh!            sedflux_vars, ep_betr, alm_fates, glc2lnd_vars, crop_vars)
            sedflux_vars, alm_fates, glc2lnd_vars, crop_vars)

       ! Reset finidat to now be finidat_interp_dest
       ! (to be compatible with routines still using finidat)
       finidat = trim(finidat_interp_dest)

       ! Write out finidat status flag
       call mpi_barrier(mpicom,ier)
       if (masterproc) then
          open (newunit=iun, file=locfn, status='unknown',  iostat=ioe)
          if (ioe /= 0) then
             call endrun(msg='ERROR failed to open file '//trim(locfn))
          end if
          write(iun,'(a)')'Successfully wrote out '//trim(locfn)
          close(iun)
          write(iulog,'(a)')' Successfully wrote finidat status file '//trim(locfn)
       end if
    end if

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call reweight_wrapup(bounds_clump, &
            glc2lnd_vars%icemask_grc(bounds_clump%begg:bounds_clump%endg))
    end do
    !$OMP END PARALLEL DO

!hh!    if(use_betr)then
!hh!      call ep_betr%set_active(bounds_proc, col_pp)
!hh!    endif
    ! ------------------------------------------------------------------------
    ! Initialize nitrogen deposition
    ! ------------------------------------------------------------------------

    if (use_cn .or. use_fates) then
       call t_startf('init_ndep')
       call ndep_init(bounds_proc, NLFilename)
       call ndep_interp(bounds_proc, atm2lnd_vars)
       call t_stopf('init_ndep')
    end if

    ! ------------------------------------------------------------------------
    ! Initialize phosphorus deposition
    ! ------------------------------------------------------------------------

    if (use_cn .or. use_fates) then
       call t_startf('init_pdep')
       call pdep_init(bounds_proc, NLFilename)
       call pdep_interp(bounds_proc, atm2lnd_vars)
       call t_stopf('init_pdep')
    end if


    ! ------------------------------------------------------------------------
    ! Initialize active history fields.
    ! ------------------------------------------------------------------------

    ! This is only done if not a restart run. If a restart run, then this
    ! information has already been obtained from the restart data read above.
    ! Note that routine hist_htapes_build needs time manager information,
    ! so this call must be made after the restart information has been read.

    if (nsrest /= nsrContinue) then
       call hist_htapes_build()
    end if

    ! ------------------------------------------------------------------------
    ! Initialize variables that are associated with accumulated fields.
    ! ------------------------------------------------------------------------

    ! The following is called for both initial and restart runs and must
    ! must be called after the restart file is read
    call atm2lnd_vars%initAccVars(bounds_proc)
    call top_as%InitAccVars(bounds_proc)
    call top_af%InitAccVars(bounds_proc)
    call veg_es%InitAccVars(bounds_proc)
    call canopystate_vars%initAccVars(bounds_proc)
    if (crop_prog) then
       call crop_vars%initAccVars(bounds_proc)
    end if
    call cnstate_vars%initAccVars(bounds_proc)

    !------------------------------------------------------------
    ! Read monthly vegetation
    !------------------------------------------------------------

    ! Even if CN is on, and dry-deposition is active, read CLMSP annual vegetation
    ! to get estimates of monthly LAI

    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
       call readAnnualVegetation(bounds_proc, canopystate_vars)
       if (nsrest == nsrStartup .and. finidat /= ' ') then
          ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will be calculated
          ! This needs to be done even if CN is on!
          call interpMonthlyVeg(bounds_proc, canopystate_vars)
       end if
    elseif ( use_fates_sp ) then
      ! If fates has satellite phenology enabled, get the monthly veg values
      ! prior to the first call to SatellitePhenology()
       call interpMonthlyVeg(bounds_proc, canopystate_vars)
    end if

    !------------------------------------------------------------
    ! Determine gridcell averaged properties to send to atm
    !------------------------------------------------------------

    if (nsrest == nsrStartup) then
       call t_startf('init_map2gc')
       call lnd2atm_minimal(bounds_proc, surfalb_vars, energyflux_vars, lnd2atm_vars)
       call t_stopf('init_map2gc')
    end if

    !------------------------------------------------------------
    ! Initialize sno export state to send to glc
    if (create_glacier_mec_landunit) then
       !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
       do nc = 1,nclumps
          call get_clump_bounds(nc, bounds_clump)

          call t_startf('init_lnd2glc')
          call lnd2glc_vars%update_lnd2glc(bounds_clump,       &
               filter(nc)%num_do_smb_c, filter(nc)%do_smb_c,   &
               init=.true.)
          call t_stopf('init_lnd2glc')
       end do
       !$OMP END PARALLEL DO
    end if

    !------------------------------------------------------------
    ! Deallocate wt_nat_patch
    !------------------------------------------------------------

    ! wt_nat_patch was allocated in initialize1, but needed to be kept around through
    ! initialize2 for some consistency checking; now it can be deallocated

    deallocate(wt_nat_patch)

    ! --------------------------------------------------------------
    ! Initialise the FATES model state structure cold-start
    ! --------------------------------------------------------------

    if ( use_fates .and. .not.is_restart() .and. finidat == ' ' .and. nsrest /= nsrBranch) then
       ! If fates is using satellite phenology mode, make sure to call the SatellitePhenology
       ! procedure prior to init_coldstart which will eventually call leaf_area_profile
       if ( use_fates_sp ) then
          !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
          do nc = 1,nclumps
             call get_clump_bounds(nc, bounds_clump)
             call SatellitePhenology(bounds_clump, &
                  filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp, &
                  waterstate_vars, canopystate_vars)
          end do
          !$OMP END PARALLEL DO
       end if
       call alm_fates%init_coldstart(canopystate_vars, soilstate_vars, frictionvel_vars)
    end if

    ! topo_glc_mec was allocated in initialize1, but needed to be kept around through
    ! initialize2 because it is used to initialize other variables; now it can be
    ! deallocated



    ! topo_glc_mec was allocated in initialize1, but needed to be kept around through
    ! initialize2 because it is used to initialize other variables; now it can be
    ! deallocated

    deallocate(topo_glc_mec)

    !------------------------------------------------------------
    ! initialize clm_bgc_interface_data_type
    call t_startf('init_elm_interface_data & pflotran')
    if (use_elm_interface) then
        call elm_interface_data%Init(bounds_proc)
        ! PFLOTRAN initialization
        if (use_pflotran) then
            call elm_pf_interface_init(bounds_proc)
        end if
    end if
    call t_stopf('init_elm_interface_data & pflotran')
    !------------------------------------------------------------

    !------------------------------------------------------------
    ! Write log output for end of initialization
    !------------------------------------------------------------

    call t_startf('init_wlog')
    if (masterproc) then
       write(iulog,*) 'Successfully initialized the land model'
       if (nsrest == nsrStartup) then
          write(iulog,*) 'begin initial run at: '
       else
          write(iulog,*) 'begin continuation run at:'
       end if
       call get_curr_date(yr, mon, day, ncsec)
       write(iulog,*) '   nstep= ',get_nstep(), ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
       write(iulog,*)
       write(iulog,'(72a1)') ("*",i=1,60)
       write(iulog,*)
    endif
    call t_stopf('init_wlog')

    call t_stopf('elm_init2')

  end subroutine initialize2

!hh!  !-----------------------------------------------------------------------
!hh!  subroutine initialize3( )
!hh!    !
!hh!    ! !DESCRIPTION:
!hh!    ! CLM initialization - third phase
!hh!    !
!hh!    ! !USES:
!hh!    use elm_varpar               , only : nlevsoi, nlevgrnd, nlevsno, max_patch_per_col
!hh!    use landunit_varcon          , only : istsoil, istcrop, istice_mec, istice_mec
!hh!    use landunit_varcon          , only : istice, istdlak, istwet, max_lunit
!hh!    use column_varcon            , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
!hh!    use elm_varctl               , only : use_vsfm, vsfm_use_dynamic_linesearch
!hh!    use elm_varctl               , only : vsfm_include_seepage_bc, vsfm_satfunc_type
!hh!    use elm_varctl               , only : vsfm_lateral_model_type
!hh!    use elm_varctl               , only : use_petsc_thermal_model
!hh!    use elm_varctl               , only : lateral_connectivity
!hh!    use elm_varctl               , only : finidat
!hh!    use decompMod                , only : get_proc_clumps
!hh!    use mpp_varpar               , only : mpp_varpar_init
!hh!    use mpp_varcon               , only : mpp_varcon_init_landunit
!hh!    use mpp_varcon               , only : mpp_varcon_init_column
!hh!    use mpp_varctl               , only : mpp_varctl_init_vsfm
!hh!    use mpp_varctl               , only : mpp_varctl_init_petsc_thermal
!hh!    use mpp_bounds               , only : mpp_bounds_init_proc_bounds
!hh!    use mpp_bounds               , only : mpp_bounds_init_clump
!hh!    use ExternalModelInterfaceMod, only : EMI_Init_EM
!hh!    use ExternalModelConstants   , only : EM_ID_VSFM
!hh!    use ExternalModelConstants   , only : EM_ID_PTM
!hh!
!hh!    implicit none
!hh!
!hh!    type(bounds_type) :: bounds_proc
!hh!    logical           :: restart_vsfm          ! does VSFM need to be restarted
!hh!
!hh!    call t_startf('elm_init3')
!hh!
!hh!    ! Is this a restart run?
!hh!    restart_vsfm = .false.
!hh!    if (nsrest == nsrStartup) then
!hh!       if (finidat == ' ') then
!hh!          restart_vsfm = .false.
!hh!       else
!hh!          restart_vsfm = .true.
!hh!       end if
!hh!    else if ((nsrest == nsrContinue) .or. (nsrest == nsrBranch)) then
!hh!       restart_vsfm = .true.
!hh!    end if
!hh!
!hh!    call mpp_varpar_init (nlevsoi, nlevgrnd, nlevsno, max_patch_per_col)
!hh!
!hh!    call mpp_varcon_init_landunit   (istsoil, istcrop, istice, istice_mec, &
!hh!           istdlak, istwet, max_lunit)
!hh!
!hh!    call mpp_varcon_init_column(icol_roof, icol_sunwall, icol_shadewall, &
!hh!      icol_road_imperv, icol_road_perv)
!hh!
!hh!    call mpp_varctl_init_vsfm(use_vsfm, vsfm_use_dynamic_linesearch, &
!hh!      vsfm_include_seepage_bc, lateral_connectivity, restart_vsfm, &
!hh!      vsfm_satfunc_type, vsfm_lateral_model_type)
!hh!
!hh!    call mpp_varctl_init_petsc_thermal(use_petsc_thermal_model)
!hh!
!hh!    call get_proc_bounds(bounds_proc)
!hh!    call mpp_bounds_init_proc_bounds(bounds_proc%begg    , bounds_proc%endg,     &
!hh!                                     bounds_proc%begg_all, bounds_proc%endg_all, &
!hh!                                     bounds_proc%begc    , bounds_proc%endc,     &
!hh!                                     bounds_proc%begc_all, bounds_proc%endc_all)
!hh!
!hh!    call mpp_bounds_init_clump(get_proc_clumps())
!hh!
!hh!    if (use_vsfm) then
!hh!       call EMI_Init_EM(EM_ID_VSFM)
!hh!    endif
!hh!
!hh!    if (use_petsc_thermal_model) then
!hh!       call EMI_Init_EM(EM_ID_PTM)
!hh!    endif
!hh!
!hh!    call t_stopf('elm_init3')
!hh!
!hh!
!hh!  end subroutine initialize3

  !-----------------------------------------------------------------------
  subroutine elm_petsc_init()
    !
    ! !DESCRIPTION:
    ! Initialize PETSc
    !
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
    use elm_varctl , only : use_vsfm
    use elm_varctl , only : lateral_connectivity
    use elm_varctl , only : use_petsc_thermal_model
#ifdef USE_PETSC_LIB
    use petscsys
#endif
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
#ifdef USE_PETSC_LIB
    PetscErrorCode        :: ierr                  ! get error code from PETSc
#endif

    if ( (.not. use_vsfm)               .and. &
         (.not. lateral_connectivity)   .and. &
         (.not. use_petsc_thermal_model) ) return

#ifdef USE_PETSC_LIB
    ! Initialize PETSc
    PETSC_COMM_WORLD = mpicom
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr);CHKERRQ(ierr)

    PETSC_COMM_SELF  = MPI_COMM_SELF
    PETSC_COMM_WORLD = mpicom
#else
    call endrun(msg='ERROR elm_petsc_init: '//&
         'PETSc required but the code was not compiled using -DUSE_PETSC_LIB')
#endif

  end subroutine elm_petsc_init


end module elm_initializeMod
