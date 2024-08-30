!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2023, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

!! \section{A toy model for the atmosphere}
!!
!! As an atmospheric component we have first coded a toy model. This model is used only
!! for debug purposes. At each timestep it simply update the AO fluxes and it exports
!! them to the ocean component.
module atmosphere_wrf

  !! We call the modules. The toy model is coded into a separate Fortran module that we
  !! do not see here. It is a custom module where the fluxes are hard-coded. In our case
  !! we have represented AO fluxes by simle analytical functions.
  !!  ESMF and NUOPC libraries will be used as well:
  use module_wrf_top
  !use module_ext_esmf, only : ioesmf_create_grid_int()
  use module_domain, only : head_grid, get_ijk_from_grid
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS    => SetServices

  implicit none

  private

  !! We define a \textsf{WRF\_FieldMetadata} for the toy model.
  !! This object that contains important info concerning the field (field names and 
  !! field staggering) has been already discussed in the ocean component.
  type WRF_Metadata
    character(4)                :: fieldName = " "
    character(50)               :: longName = " "
    type(ESMF_StaggerLoc)       :: meshloc = ESMF_STAGGERLOC_CENTER
  end type

  type WRF_Grid
    integer                     :: localPet
    integer                     :: petCount
  end type

  integer :: WRF_numOfImportFields
  integer :: WRF_numOfExportFields
  type(WRF_Metadata), allocatable, save :: WRF_FieldMetadata(:)
  type(WRF_Grid), save :: WrfToEsmf_Grid

  public SetServices

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  !! We are familiar with the creation of the atmospheric component with the 
  !! NUOPC layer. As done for the other components the derive/specialize commands are invoked.
  !! The operations are the same for the ocean component (Advertise, Realize, Advance and Finalize).
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! derive from NUOPC_Model
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! specialize model
    call NUOPC_CompSpecialize(model, specLabel=label_Advertise, &
      specRoutine=Advertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_RealizeProvided, &
      specRoutine=Realize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=Advance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
      specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
  !! The following subroutines register
  !! the members and methods for the the atmospheric toy model object that enters as first
  !! argument. We do not comment extensively such subroutines because we have coded them
  !! very similarly to the SHYFEM component.
  !!
  !! \subsection{Advertise}

  subroutine Advertise(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    type(ESMF_VM)           :: vm
    integer                 :: var
    integer                 :: num
    integer                 :: mpicomtmp

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
      
    call ESMF_VMGet(vm, mpiCommunicator=mpicomtmp, rc=rc)

    call wrf_set_dm_communicator( mpicomtmp )

    call wrf_init(  no_init1=.TRUE. )
    call ESMF_LogWrite("Initialized ATM", ESMF_LOGMSG_INFO, rc=rc)

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    allocate( WRF_FieldMetadata(11) )
    
#undef WITHFIELDS_MOMENTUMFLUX
#undef WITHFIELDS_HEATFLUX
#undef WITHFIELDS_MASSFLUX
#undef WITHFIELDS_SST
    !! We have added one macro to rapidly decouple the atmosphere component from the rest
    !! of the earth system. Disabling the following macros,
    !! will result in an atmosphere component that does not advertise any importable
    !! Fields. We use this if we want to drive the model independently, which is the 
    !! case for our toy model (for instance it does not depends on the ocean state).
    num = 0

#ifdef WITHFIELDS_SST
    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="umap", &
      longName="atmospheric_unmapped_points", &
      meshloc=ESMF_STAGGERLOC_CENTER)

    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="sst", &
      longName="sea_surface_temperature", &
      meshloc=ESMF_STAGGERLOC_CENTER)
#endif

    WRF_numOfImportFields = num

#ifdef WITHFIELDS_MOMENTUMFLUX
    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="pmsl", &
      longName="air_pressure_at_sea_level", &
      meshloc=ESMF_STAGGERLOC_CENTER) 

    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="smes", &
      longName="surface_downward_eastward_stress", &
      meshloc=ESMF_STAGGERLOC_CENTER)

    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="smns", &
      longName="surface_downward_northward_stress", &
      meshloc=ESMF_STAGGERLOC_CENTER)
#endif
#ifdef WITHFIELDS_HEATFLUX
    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="rsns", &
      longName="surface_net_downward_shortwave_flux", &
      meshloc=ESMF_STAGGERLOC_CENTER)

    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="rlns", &
      longName="surface_net_downward_longwave_flux", &
      meshloc=ESMF_STAGGERLOC_CENTER)

    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="stsh", &
      longName="surface_downward_sensible_heat_flux_in_air", &
      meshloc=ESMF_STAGGERLOC_CENTER)

    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="stlh", &
      longName="surface_downward_latent_heat_flux_in_air", &
      meshloc=ESMF_STAGGERLOC_CENTER)
#endif
#ifdef WITHFIELDS_MASSFLUX
    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="evap", &
      longName="evaporation_flux", &
      meshloc=ESMF_STAGGERLOC_CENTER)

    num = num + 1
    WRF_FieldMetadata(num) = WRF_Metadata(fieldName="prec", &
      longName="precipitation_flux", &
      meshloc=ESMF_STAGGERLOC_CENTER)
#endif

    WRF_numOfExportFields = num - WRF_numOfImportFields 

    do var=1,WRF_numOfImportFields

      call NUOPC_Advertise(importState, &
        StandardName=WRF_FieldMetadata(var)%longName, &
        name=WRF_FieldMetadata(var)%fieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

    do var=WRF_numOfImportFields+1, &
      WRF_numOfImportFields+WRF_numOfExportFields

      call NUOPC_Advertise(exportState, &
	StandardName=WRF_FieldMetadata(var)%longName, &
	name=WRF_FieldMetadata(var)%fieldName, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
	line=__LINE__, &
	file=__FILE__)) &
	return  ! bail out

    enddo

  end subroutine

  !-----------------------------------------------------------------------------

  !!
  !! \subsection{Realize}
  !!
  !! Nothing different with respect to the ocean component for the Advertise subroutine. 
  !! Fields that were previously advertised should now be realized. 
  !! Realizing a field means that an \textsf{ESMF\_Field} object is created and 
  !! it is added to the appropriate \textsf{ESMF\_State}, either import or export.  
  subroutine Realize(model, rc)

!    use CONFIG, only : cfgHandler

    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)               :: importState, exportState
    type(ESMF_Field)               :: field
    type(ESMF_Grid)                :: gridIn
    type(ESMF_Grid)                :: gridOut
    double precision, pointer      :: fieldPtr(:,:)
!    integer                        :: fieldtotalLBnd(2)
!    integer                        :: fieldtotalUBnd(2)
!    integer                        :: fieldtotal_count(2)
    integer                        :: var
!    integer                        :: i,j,k,localPet,petCount

    rc = ESMF_SUCCESS

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call WRF_GridGet(gridIn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    call ESMF_LogWrite("  Get Grid ATM", ESMF_LOGMSG_INFO, rc=rc)


!-----------------------------------------------------------------------
!     Get gridded component 
!-----------------------------------------------------------------------

!    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, 
!      msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=FILENAME)) &
!      return
!    call ESMF_VMGetCurrent(vm, rc=rc)
!    
!    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return
      
    !! The field creation has been extensively explained for the ocean component.
    !! We show how to create import fields, even if they are not used in the 
    !! toy model.
    do var=1,WRF_numOfImportFields

      field = ESMF_FieldCreate(grid=gridIn, typekind=ESMF_TYPEKIND_R8, &
        staggerloc=WRF_FieldMetadata(var)%meshloc, &
        name=WRF_FieldMetadata(var)%fieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_Realize(importState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_LogWrite(trim("ATM: Creating Field ")// &
        trim(WRF_FieldMetadata(var)%fieldName), ESMF_LOGMSG_INFO, rc=rc)

    enddo

    !! More important is the creation of export fields where the fluxes 
    !! are stored. We loop over these fields and we call \textsf{ESMF\_FieldCreate} and
    !! \textsf{ESMF\_Realize} as usual.
    do var=WRF_numOfImportFields+1, &
      WRF_numOfImportFields+WRF_numOfExportFields

      field = ESMF_FieldCreate(grid=gridIn, typekind=ESMF_TYPEKIND_R8, &
        staggerloc=WRF_FieldMetadata(var)%meshloc, &
        name=WRF_FieldMetadata(var)%fieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
	line=__LINE__, &
	file=__FILE__)) &
	return  ! bail out

      call NUOPC_Realize(exportState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
	line=__LINE__, &
	file=__FILE__)) &
	return  ! bail out

      !! The model has not been initialized. This is dangerous becouse at the
      !! first timestep the ocean model would run with uninitialized fluxes.
      !! With the following commands we initialize the atmospheric component.
      !! We retrieve the field pointer with \textsf{ESMF\_FieldGet}. We allocate the
      !! arrays in \textsf{toy\_initialize}, we initialize them in \textsf{toy\_run}
      !! giving a zero time. Finally the \textsf{WRF\_FieldSet} command
      !! copy the field pointer to the WRF flux variables.
      call ESMF_FieldGet(field, localDe=0, farrayPtr=fieldPtr, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call WRF_FieldGet(fieldPtr, &
        WRF_FieldMetadata(var)%fieldName, rc)

    enddo

  end subroutine

  !-----------------------------------------------------------------------------

  !!
  !! \subsection{Advance the toy model}
  !!
  !! The atmospheric model advances of one timestep. It means that the atmospheric 
  !! variables and/or fluxes are re-evalauted at the current time
  !! \[
  !! F_{ao} = F_{ao}(t^{n+1})
  !! \]
  !! It is important to note that the we are only updating the component but 
  !! the timestepping scheme stays explicit. 
  !! The flux values infact will go to the ocean component only at the next
  !! AO timestep.
  subroutine Advance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    character(len=160)          :: msgString
    character(len=160)          :: dateString, petString
    double precision            :: timeStepSec
    type(ESMF_Field)            :: field

    double precision, pointer   :: fieldPtr(:,:)
    double precision, pointer   :: unmappedPointsPtr(:,:)    
    integer                     :: fieldtotalLBnd(2)
    integer                     :: fieldtotalUBnd(2)
    integer                     :: fieldtotal_count(2)
    integer                     :: var	
    integer                     :: currYear, currMonth, currDay, &
                                   currHour, currMinute, &
                                   currSecond, currTimeSec
!    character(len=256)          :: timeStr

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
    call ESMF_TraceRegionEnter("ATM:Advance")
#endif

    rc = ESMF_SUCCESS

    !! To update the fields to be exported there is some work to do. But this involves 
    !! the same work already done for the ocean component.
    !! We go into the object hierarcy. We query the model with \textsf{NUOPC\_ModelGet()}
    !! and we print the state to screen.
    !! Then we access the state with \textsf{ESMF\_StateGet()} and we retrieve the fields
    !! Once we have the field object we use \textsf{ESMF\_FieldGet()} to 
    !! extract the pointer to the data array.
    !! Other useful information about the bounds of the field can
    !! be also retrieved and used to bound the loops that updates the field.
    ! query for clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !call ESMF_StatePrint(exportState, rc=rc)

    !! Meanwhile we need also to clock object. Apart from printing to the 
    !! PET file useful information about the situation of the component, we 
    !! retrieve from the clock object,
    !! the \textsf{currTime} object and the \textsf{timeStep} object. With
    !! te timestep we compute the new internal time of the atmospheric component. 
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing ATM from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, &
      timeStep=timeStep, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeIntervalGet(timeStep, s_r8=timeStepSec, rc=rc)

    !! To write the output fields we also write retrieve the current date in two
    !! differen unit. One is used to name the output file in an
    !! ordered fashion. The second one is to check if the current
    !! date in seconds correspond to an output tick. In that case
    !! we write the fields to file. 
    call ESMF_TimeGet(currTime, yy=currYear, mm=currMonth, &
      dd=currDay, h=currHour, m=currMinute, s=currSecond, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_TimeGet(currTime, yy=currYear, s=currTimeSec, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(dateString, & 
      "(I4,'-',I2.2,'-',I2.2,'_',I2.2,':',I2.2,':',I2.2)") & 
      currYear,currMonth, currDay, currHour, currMinute, currSecond
    write(petString, "('.',I1,'.',I1)") &
      WrfToEsmf_Grid%petCount, WrfToEsmf_Grid%localPet

    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="---------------------> to: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    time = time + timeStepSec

    !! The export state must be updated. A loop on the export fields is performed. 
    !! We can retrieve what is needed from the export state.
    !! To say that we want to retrieve a field item (of \textsf{ESMF\_Field} type) we 
    !! add the keyworld \textsf{field}. With \textsf{itemName} we choose
    !! the specific field among the ones that we have created. The next step is to run
    !! the toy model, that is we simply re-evaluate at the new time the arrays of fluxes. The
    !! fields are then assigned to the updated arrays:
    do var=1, WRF_numOfImportFields-1
      call ESMF_StateGet(importState, field=field, &
        itemName="umap", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_FieldGet(field, localDe=0, farrayPtr=unmappedPointsPtr, &
        totalLBound=fieldtotalLBnd, totalUBound=fieldtotalUBnd, &
        totalCount=fieldtotal_count, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!!      if ( currTimeSec>0. .and. mod(currTimeSec,3600).lt.1 ) then
!!        call WRF_FieldWrite(field, &
!!          trim(WRF_FieldMetadata(var)%fieldName)//trim(dateString) &
!!          //trim(".vtk"), rc)
!!        if(ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU, &
!!          line=__LINE__, &
!!          file=__FILE__)) &
!!          return  ! bail out
!!      endif
    enddo

    do var=2, WRF_numOfImportFields

      call ESMF_StateGet(importState, field=field, &
        itemName=WRF_FieldMetadata(var)%fieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

      call ESMF_FieldGet(field, localDe=0, farrayPtr=fieldPtr, &
       totalLBound=fieldtotalLBnd, totalUBound=fieldtotalUBnd, &
       totalCount=fieldtotal_count, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!      if ( currTimeSec>0. .and. mod(currTimeSec,idtatm).lt.1 ) then
      if ( currTimeSec>0. .and. mod(currTimeSec,86400).lt.1 ) then
	call WRF_FieldWrite(field, &
	  trim(WRF_FieldMetadata(var)%fieldName)//trim(dateString) &
	  //trim(petString)//trim(".vtk"), rc)
	if(ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU, &
	  line=__LINE__, &
	  file=__FILE__)) &
	  return  ! bail out
      endif


      call WRF_FieldSet(fieldPtr, &
        WRF_FieldMetadata(var)%fieldName, unmappedPointsPtr, rc)

    enddo

    head_grid%start_subtime = currTime
    head_grid%stop_subtime = currTime + timeStep
  
!    call wrf_timetoa( head_grid%start_subtime, timeStr )
    if( WrfToEsmf_Grid%localPet == 0 ) then
      print *
      print *, 'wrf_component_run: start nuopc timestep at date ', &
        TRIM(dateString)
      print *
    endif

    call wrf_run

    !! The export state must be updated. A loop on the export fields is performed. 
    !! We can retrieve what is needed from the export state.
    !! To say that we want to retrieve a field item (of \textsf{ESMF\_Field} type) we 
    !! add the keyworld \textsf{field}. With \textsf{itemName} we choose
    !! the specific field among the ones that we have created. The next step is to run
    !! the toy model, that is we simply re-evaluate at the new time the arrays of fluxes. The
    !! fields are then assigned to the updated arrays:
    do var=WRF_numOfImportFields+1, &
      WRF_numOfImportFields+WRF_numOfExportFields

      call ESMF_StateGet(exportState, field=field, &
        itemName=WRF_FieldMetadata(var)%fieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call ESMF_FieldGet(field, localDe=0, farrayPtr=fieldPtr, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call WRF_FieldGet(fieldPtr, &
        WRF_FieldMetadata(var)%fieldName, rc)

    enddo


#ifdef NUOPC_TRACE
    call ESMF_TraceRegionExit("ATM:Advance")
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  !!
  !! \subsection{Finalize the toy model}
  !!
  !! We register a method for the finalization of the toy model. We simply clean 
  !! the memory through the deallocations of the data structures.   
  subroutine Finalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    integer              :: var 

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
    call ESMF_TraceRegionEnter("ATM:Finalize")
#endif

    rc = ESMF_SUCCESS

    ! HERE THE MODEL IS FINALIZED: 
    do var=WRF_numOfImportFields+1, &
      WRF_numOfImportFields+WRF_numOfExportFields

!      call toy_finalize(WRF_FieldMetadata(var)%fieldName)

    enddo

    deallocate(WRF_FieldMetadata)

    call wrf_finalize(.true.)

    call ESMF_LogWrite("Finalized ATM", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef NUOPC_TRACE
    call ESMF_TraceRegionExit("ATM:Finalize")
#endif

  end subroutine

  !-----------------------------------------------------------------------------

  !! Eventually we write the regridded fields to files.
  !! This can be helpful for debugging and checking the interpolations.
  !! We can write such a file with ESMF subroutine \textsf{WRF\_FieldWrite}
  !! but this works only with the third party library PARALLELIO (PIO).
  !! Moreover the only format allowed when this manual was written was netcdf
  !! (ugrid). We have preferred to do use vtk format to visualize the data,
  !! as done with the mesh. This lead to only one type of file outputted.
  !! Vtk files can be visualize nicely with Paraview.
  subroutine WRF_FieldWrite(field, filename, rc)

    type(ESMF_Field), intent(in)  :: field
    character(len=*), intent(in)  :: filename
    integer, intent(out)          :: rc

    ! local variables
    double precision, pointer   :: ptr(:,:)
    character(20) :: str
    integer :: ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme, &
               ips, ipe, jps, jpe, kps, kpe
    integer :: i, j, ii

    rc = ESMF_SUCCESS

    !! We recovet the field with the usual call
    call ESMF_FieldGet(field, localDe=0, farrayPtr=ptr, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call get_ijk_from_grid(head_grid, ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      ips, ipe, jps, jpe, kps, kpe)

    !! We spend some word about the vtk format. We use
    !! three digits to represents point coordinates, that means
    !! that you can reach mm of resolution in you mesh but not beyond.
    !! Paraview does not like zeros as float 0.00 and I had to manually
    !! write zero as an integer. For the variable precision we assume
    !! that three digits are also sufficient.
    open(1, file = filename, status = 'unknown')
    write(1,'(a)') "# vtk DataFile Version 3.0"
    write(1,'(a)') "This file was generated by NUOPC"
    write(1,'(a)') "ASCII"
    write(1,'(a)') "DATASET UNSTRUCTURED_GRID"
    write(1,'(a)', advance='no') "POINTS ";
    write(1, *) (ipe-ips+1)*(jpe-jps+1), " double"
    if ( head_grid%xlong(1,1) .ne. head_grid%xlong(2,1) ) then
      do j=jps,jpe!-1
      do i=ips,ipe!-1
        write(str,'(f11.3)') head_grid%xlong(i,j)
        str = trim(adjustl(str))
        do ii = len_trim(str),1,-1
          if (str(ii:ii)/="0") exit
        enddo
        if (str(ii:ii)==".") ii=ii-1
        write(1,'(a,a)', advance="no") str(1:ii), " "
        write(str,'(f11.3)') head_grid%xlat(i,j) 
        str = trim(adjustl(str))
        do ii = len_trim(str),1,-1
          if (str(ii:ii)/="0") exit
        enddo
        if (str(ii:ii)==".") ii=ii-1
        write(1,'(a,a)') str(1:ii), " 0"
      enddo
      enddo
    else  
      do j=jps,jpe!-1
      do i=ips,ipe!-1      
        write(str,'(f11.3)') &
          (real(i)-1) * head_grid%dx !+  0.5*head_grid%dx
        str = trim(adjustl(str))
        do ii = len_trim(str),1,-1
          if (str(ii:ii)/="0") exit
        enddo
        if (str(ii:ii)==".") ii=ii-1
        write(1,'(a,a)', advance="no") str(1:ii), " "
        write(str,'(f11.3)') &
          (real(j)-1) * head_grid%dy !+  0.5*head_grid%dy       
        str = trim(adjustl(str))
        do ii = len_trim(str),1,-1
          if (str(ii:ii)/="0") exit
        enddo
        if (str(ii:ii)==".") ii=ii-1
        write(1,'(a,a)') str(1:ii), " 0"	      
      enddo
      enddo
    endif  
    write(1,'(a)', advance='no') "CELLS "
    write(1, *) (ipe-ips)*(jpe-jps), (ipe-ips)*(jpe-jps)*5      
    !do j=jps,jpe-2
    !do i=ips,ipe-2
    do j=1,(jpe-jps)
    do i=1,(ipe-ips)  
      write(1,"(i6X,i6X,i6X,i6X,i6)") 4, &
       (ipe-ips+1)*(j-1)+i-1, (ipe-ips+1)*(j-1)+i, &
       (ipe-ips+1)*(j)+i, (ipe-ips+1)*(j)+i-1
    end do
    enddo
    write(1,'(a)', advance="no") "CELL_TYPES "
    write(1, *) (ipe-ips)*(jpe-jps)
    do i=1,(ipe-ips)*(jpe-jps)
      write(1,*) 9
    end do
    write(1,'(a)', advance="no") "POINT_DATA "
    write(1, *) (ipe-ips+1)*(jpe-jps+1)
    write(1,'(a)') "SCALARS _NODE_NUM double 1"
    write(1,'(a)') "LOOKUP_TABLE default"
    do j=jps,jpe
    do i=ips,ipe
      write(str,'(f11.4)') ptr(i,j)
      str = trim(adjustl(str))
      do ii = len_trim(str),1,-1
        if (str(ii:ii)/="0") exit
      enddo
      if (str(ii:ii)==".") ii=ii-1
      write(1,'(a,a)', advance="no") str(1:ii), " "
    end do  
    end do
    close(1)

  end subroutine

  !-----------------------------------------------------------------------------

  !! We use this routine to read the WRF grid. This part is copied from RegESM
  !! Earth System model of []
  subroutine WRF_GridGet(WRF_grid, rc)

    type(ESMF_Grid), intent(out)    :: WRF_grid
    integer, intent(out) 	    :: rc

    type(ESMF_CoordSys_Flag)        :: WRF_coordSys  
    integer :: i, j, k, pe
    integer :: ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme, &
               ips, ipe, jps, jpe, kps, kpe
    integer :: petCount 
!    integer :: numprocsX, numprocsY
    integer :: lbnd(3), ubnd(3)
    integer :: minIndex(2), maxIndex(2)
    integer, allocatable :: ipatchStarts(:), jpatchStarts(:)
    integer, allocatable :: ipatchEnds(:), jpatchEnds(:)
    integer, allocatable :: deBlockList(:,:,:)
    real(ESMF_KIND_R8), pointer :: ptrX(:,:), ptrY(:,:)!, ptrA(:,:)
    integer(ESMF_KIND_I4), pointer :: ptrM(:,:)
!    character(ESMF_MAXSTR) :: cname, name
!    type(domain), pointer :: grid
    type(ESMF_VM) :: vm
    type(ESMF_DistGrid) :: distGrid

    rc = ESMF_SUCCESS        
        
!-----------------------------------------------------------------------
!     Get gridded component 
!-----------------------------------------------------------------------

!    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, 
!      msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=FILENAME)) &
!      return
    call ESMF_VMGetCurrent(vm, rc=rc)
    
    call ESMF_VMGet(vm, localPet=WrfToEsmf_Grid%localPet, &
      petCount=WrfToEsmf_Grid%petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    petCount = WrfToEsmf_Grid%petCount                                       

    !! We get the WRF variables related with grid partitioning: they
    !! are defined as follow:
    !! ids, ide, jds, jde, kds, kde => domain extent
    !! ims, ime, jms, jme, kms, kme => memory extent
    !! ips, ipe, jps, jpe, kps, kpe => patch extent
    call get_ijk_from_grid(head_grid, ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      ips, ipe, jps, jpe, kps, kpe)

!    print *, "Conversion to ESMF distributed mesh: done"
!    print *, "          is   ie   js   je   ks   ke"
!    print *, "domain:", ids, ide, jds, jde, kds, kde
!    print *, "memory:", ims, ime, jms, jme, kms, kme
!    print *, "patch :", ips, ipe, jps, jpe, kps, kpe
!    print *, "PETs  :", localPet, petCount
 
    !! We collect the index of all patches that are on other
    !! PET.
    if (.not. allocated(ipatchStarts)) then
      allocate(ipatchStarts(0:petCount-1))
    end if

    if (.not. allocated(jpatchStarts)) then
      allocate(jpatchStarts(0:petCount-1))
    end if

    call ESMF_VMAllGatherV(vm, sendData=(/ ips /), sendCount=1, &
      recvData=ipatchStarts, &
      recvCounts=(/ (1, k = 0, petCount-1) /), &
      recvOffsets=(/ (k, k = 0, petCount-1) /), &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    call ESMF_VMAllGatherV(vm, sendData=(/ jps /), sendCount=1, &
      recvData=jpatchStarts, &
      recvCounts=(/ (1, k = 0, petCount-1) /), &
      recvOffsets=(/ (k, k = 0, petCount-1) /), &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    if (.not. allocated(ipatchEnds)) then
      allocate(ipatchEnds(0:petCount-1))
    end if

    if (.not. allocated(jpatchEnds)) then
      allocate(jpatchEnds(0:petCount-1))
    end if

    call ESMF_VMAllGatherV(vm, sendData=(/ ipe /), &
      sendCount=1, recvData=ipatchEnds, &
      recvCounts=(/ (1, k = 0, petCount-1) /), &
      recvOffsets=(/ (k, k = 0, petCount-1) /), &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    call ESMF_VMAllGatherV(vm, sendData=(/ jpe /), &
      sendCount=1, recvData=jpatchEnds, &
      recvCounts=(/ (1, k = 0, petCount-1) /), &
      recvOffsets=(/ (k, k = 0, petCount-1) /), &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    !! We create an \textsf{ESMF\_DistGrid} from a single logically 
    !! rectangular tile with decomposition specified by \textsf{deBlockList}. 
    !! We fill the last array and then we search for the min/max index.
    if (.not.allocated(deBlockList)) then
      allocate(deBlockList(2,2,petCount))
    end if

    do pe = 1, petCount
      deBlockList(1,1,pe) = ipatchStarts(pe-1)
      deBlockList(1,2,pe) = ipatchEnds(pe-1)
      deBlockList(2,1,pe) = jpatchStarts(pe-1) 
      deBlockList(2,2,pe) = jpatchEnds(pe-1)
    end do

    minIndex = (/ minval(ipatchStarts), minval(jpatchStarts) /)
    maxIndex = (/ maxval(ipatchEnds), maxval(jpatchEnds) /)

    distGrid = ESMF_DistGridCreate(minIndex=minIndex, &
       maxIndex=maxIndex, &
       deBlockList=deBlockList, &
       rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

!-----------------------------------------------------------------------
!     Define component grid (dot and cross points)
!-----------------------------------------------------------------------

!    if (.not. allocated(models(Iatmos)%mesh)) then
!      allocate(models(Iatmos)%mesh(1))
!      models(Iatmos)%mesh(1)%gtype = Icross
!    end if

!-----------------------------------------------------------------------
!     Create ESMF Grid
!-----------------------------------------------------------------------

    !! The coordinate system is asked to the SHYFEM configuration file. Please
    !! note that the user must assure that a single coordinate system
    !! among all the the components. Two choice are possible in SHYFEM :
    !! \texttt{isphe=1} stands for spherical (in degrees) or 
    !! \texttt{isphe=0} stands for cartesian. Please check the SHYFEM manual
    !! for the details.
    if ( head_grid%xlong(1,1) .ne. head_grid%xlong(2,1) ) then
      WRF_coordSys = ESMF_COORDSYS_SPH_DEG
    else
      WRF_coordSys = ESMF_COORDSYS_CART
    endif

    WRF_grid = ESMF_GridCreate(distgrid=distGrid, &
      coordSys=WRF_coordSys, indexflag=ESMF_INDEX_GLOBAL, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

!-----------------------------------------------------------------------
!     Allocate coordinates 
!-----------------------------------------------------------------------

    call ESMF_GridAddCoord(WRF_grid, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

!-----------------------------------------------------------------------
!     Allocate items for masking
!-----------------------------------------------------------------------

    call ESMF_GridAddItem(WRF_grid, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

!-----------------------------------------------------------------------
!     Set mask value for land and ocean 
!-----------------------------------------------------------------------
!
!      models(Iatmos)%isLand = 1
!      models(Iatmos)%isOcean = 0
!
!-----------------------------------------------------------------------
!     Allocate items for grid area 
!-----------------------------------------------------------------------
!
!    call ESMF_GridAddItem(models(Iatmos)%grid, &
!      staggerLoc=ESMF_STAGGERLOC_CENTER, &
!      itemflag=ESMF_GRIDITEM_AREA, &
!      rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
!                             line=__LINE__, file=FILENAME)) return
!
!-----------------------------------------------------------------------
!     Get pointers and set coordinates for the grid 
!-----------------------------------------------------------------------
! 
    call ESMF_GridGetCoord(WRF_grid, localDE=0, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=1, &
      farrayPtr=ptrX, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    call ESMF_GridGetCoord(WRF_grid, localDE=0, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=2, &
      computationalLBound=lbnd, &
      computationalUBound=ubnd, &
      farrayPtr=ptrY, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    call ESMF_GridGetItem(WRF_grid, localDE=0, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, &
      farrayPtr=ptrM, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

!      call ESMF_GridGetItem (models(Iatmos)%grid,                       &
!                             localDE=0,                                 &
!                             staggerLoc=ESMF_STAGGERLOC_CENTER,         &
!                             itemflag=ESMF_GRIDITEM_AREA,               &
!                             farrayPtr=ptrA,                            &
!                             rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
!                             line=__LINE__, file=FILENAME)) return

!-----------------------------------------------------------------------
!     Fill the pointers    
!-----------------------------------------------------------------------
!
!      if (debugLevel > 0) then
!          write(*,30) localPet, 1, adjustl("DAT/ATM/GRD/"//name),       &
!                      lbound(head_grid%xlong, dim=1),                   &
!                      ubound(head_grid%xlong, dim=1),                   &
!                      lbound(head_grid%xlong, dim=2),                   &
!                      ubound(head_grid%xlong, dim=2),                   &
!                      head_grid%bdy_mask(1), head_grid%bdy_mask(2),     &
!                      head_grid%bdy_mask(3), head_grid%bdy_mask(4)
!      end if

    if (WRF_coordSys == ESMF_COORDSYS_SPH_DEG ) then
      do i = lbnd(1), ubnd(1)
        do j = lbnd(2), ubnd(2)
          ptrX(i,j) = head_grid%xlong(i,j)
          ptrY(i,j) = head_grid%xlat(i,j)
          ptrM(i,j) = head_grid%landmask(i,j) 
!         ptrA(i,j) = head_grid%dx*head_grid%dy
        end do
      end do
    else
      do i = lbnd(1), ubnd(1)
        do j = lbnd(2), ubnd(2)
          ptrX(i,j) = (real(i)-1) * head_grid%dx !+  0.5*head_grid%dx
          ptrY(i,j) = (real(j)-1) * head_grid%dy !+  0.5*head_grid%dy     
          ptrM(i,j) = head_grid%landmask(i,j)
!         ptrA(i,j) = head_grid%dx*head_grid%dy
        end do
      end do      
    endif

    if (head_grid%bdy_mask(2)) then
      ptrX(ubnd(1),:) = ptrX(ubnd(1)-1,:)+(ptrX(ubnd(1)-1,:)-ptrX(ubnd(1)-2,:))
      ptrY(ubnd(1),:) = ptrY(ubnd(1)-1,:)+(ptrY(ubnd(1)-1,:)-ptrY(ubnd(1)-2,:))
    end if

    if (head_grid%bdy_mask(4)) then
      ptrX(:,ubnd(2)) = ptrX(:,ubnd(2)-1)+(ptrX(:,ubnd(2)-1)-ptrX(:,ubnd(2)-2))
      ptrY(:,ubnd(2)) = ptrY(:,ubnd(2)-1)+(ptrY(:,ubnd(2)-1)-ptrY(:,ubnd(2)-2))
    end if

!-----------------------------------------------------------------------
!     Nullify pointers 
!-----------------------------------------------------------------------

    if (associated(ptrX)) then
      nullify(ptrX)
    end if
    if (associated(ptrY)) then
      nullify(ptrY)
    end if
!      if (associated(ptrM)) then
!        nullify(ptrM)
!      end if
!      if (associated(ptrA)) then
!        nullify(ptrA)
!      end if
!
!-----------------------------------------------------------------------
!     Deallocate arrays    
!-----------------------------------------------------------------------

    if (allocated(ipatchStarts)) then
      deallocate(ipatchStarts)
    end if
    if (allocated(jpatchStarts)) then
      deallocate(jpatchStarts)
    end if
    if (allocated(ipatchEnds)) then
      deallocate(ipatchEnds)
    end if
    if (allocated(jpatchEnds)) then
      deallocate(jpatchEnds)
    end if
    if (allocated(deBlockList)) then
      deallocate(deBlockList)
    end if

!-----------------------------------------------------------------------
!     Assign grid to gridded component 
!-----------------------------------------------------------------------
!
!      call ESMF_GridCompSet(gcomp, grid=models(Iatmos)%grid, rc=rc)
!      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
!                             line=__LINE__, file=FILENAME)) return      
!
!-----------------------------------------------------------------------
!     Debug: write out component grid in VTK format 
!-----------------------------------------------------------------------

!    call ESMF_GridWriteVTK(WRF_grid, &
!      staggerLoc=ESMF_STAGGERLOC_CENTER, &
!      filename="wrf_grid", &
!      rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine WRF_FieldSet(fieldPtr, fieldName, unmappedPointsPtr, rc)

    double precision, pointer, intent(inout) :: fieldPtr(:,:)
    character(len=*), intent(in)             :: fieldName
    double precision, pointer, intent(in)    :: unmappedPointsPtr(:,:)
    integer, intent(out)                     :: rc

    integer :: ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme, &
   	       ips, ipe, jps, jpe, kps, kpe
    integer :: i, j


    rc = ESMF_SUCCESS

!-----------------------------------------------------------------------
!     Get variables related with grid partitioning
!     ids, ide, jds, jde, kds, kde => domain extent
!     ims, ime, jms, jme, kms, kme => memory extent
!     ips, ipe, jps, jpe, kps, kpe => patch extent
!-----------------------------------------------------------------------
!
    call get_ijk_from_grid(head_grid, ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      ips, ipe, jps, jpe, kps, kpe)

    !! We assign the ocean state. This has to be done with some
    !! care. The land-points should be excluded as
    !! well as ocean points that fall outside the atmospheric grid.
    !! Identifiying the list of points belonging to the first case
    !! is easy thanks to \textsf{head\_grid%landsea\_mask}. The list
    !! of non-masked points that fall outside the ocean boundaries is
    !! contained in the list \textsf{unmappedDstList}. This list is not
    !! available in ESMF 8.6.0. Here we have used the auxiliary field
    !! that flag unmapped points.
    select case (fieldName)
      case ("sst")
        do i = ips, ipe
          do j = jps, jpe
            if (unmappedPointsPtr(i,j)  .ge. 0.99999D0 .and. & 
                head_grid%landmask(i,j) .eq. 0) then
!               head_grid%sst(i,j) = fieldPtr(i,j) +  273.15
	       head_grid%tsk(i,j) = fieldPtr(i,j) +  273.15
	    else
!              head_grid%sst(i,j) = 28.D0 +  273.15
              head_grid%tsk(i,j) = 28.D0 +  273.15	      
            end if
          end do
        end do
      case default
        call ESMF_LogWrite("  ATM unknown field name", &
    	  ESMF_LOGMSG_INFO, rc=rc)
        rc = 1
    end select

  end subroutine

  subroutine WRF_FieldGet(fieldPtr, fieldName, rc)

    double precision, pointer, intent(in) :: fieldPtr(:,:)
    character(len=*), intent(in)          :: fieldName
    integer, intent(out)                  :: rc
 
    real(ESMF_KIND_R8) :: cff, cff1, cff2, cff3   
    real(ESMF_KIND_R8), parameter :: eps=1.0e-10
    integer :: ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme, &
               ips, ipe, jps, jpe, kps, kpe
    integer :: i, j

    rc = ESMF_SUCCESS

!-----------------------------------------------------------------------
!     Get variables related with grid partitioning
!     ids, ide, jds, jde, kds, kde => domain extent
!     ims, ime, jms, jme, kms, kme => memory extent
!     ips, ipe, jps, jpe, kps, kpe => patch extent
!-----------------------------------------------------------------------
!
    call get_ijk_from_grid(head_grid, ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      ips, ipe, jps, jpe, kps, kpe)

    select case (fieldName)
      case ("pmsl")
        do i = ips, ipe 
          do j = jps, jpe
            fieldPtr(i,j) = head_grid%PSFC(i,j) * &
              exp( (9.81*head_grid%ht(i,j)) / (287.0*head_grid%T2(i,j)* & 
              (1.0+0.61*head_grid%Q2(i,j)) + 1e-12 ) ) !lrp: 1e-12 just to desingularize in case of T2=0 (tc testcase)
          end do
        end do	      
      case ("smes")
        do i = ips, ipe
          do j = jps, jpe
            cff1 = 1.0/(head_grid%alt(i,1,j)+eps)
!            cff2 = 2.0/(((head_grid%u_2(i,1,j)+head_grid%u_2(i+1,1,j))**2+&
!              (head_grid%v_2(i,1,j)+head_grid%v_2(i,1,j+1))**2)**0.5+eps)
!            cff3 = 0.5*(head_grid%u_2(i,1,j)+head_grid%u_2(i+1,1,j))*head_grid%cosa(i,j)-&
!              0.5*(head_grid%v_2(i,1,j)+head_grid%v_2(i,1,j+1))*head_grid%sina(i,j)
            cff2 = 2.0/(((head_grid%U10(i,j)+head_grid%U10(i+1,j))**2+&
              (head_grid%V10(i,j)+head_grid%V10(i,j+1))**2)**0.5+eps)
            cff3 = 0.5*(head_grid%U10(i,j)+head_grid%U10(i+1,j))*head_grid%cosa(i,j)-&
              0.5*(head_grid%V10(i,j)+head_grid%V10(i,j+1))*head_grid%sina(i,j)
	    cff = cff1*cff2*(head_grid%UST(i,j)**2)*cff3
	    fieldPtr(i,j) = cff
          end do
	end do
      case ("smns")
        ! following code inherited from COAWST modelling system
        ! for more information please refer to Warner et al., 2010   
        do i = ips, ipe
          do j = jps, jpe
            cff1 = 1.0/(head_grid%alt(i,1,j)+eps)
!            cff2 = 2.0/(((head_grid%u_2(i,1,j)+head_grid%u_2(i+1,1,j))**2+&
!              (head_grid%v_2(i,1,j)+head_grid%v_2(i,1,j+1))**2)**0.5+eps)
!            cff3 = 0.5*(head_grid%v_2(i,1,j)+head_grid%v_2(i,1,j+1))*head_grid%cosa(i,j)+&
!              0.5*(head_grid%u_2(i,1,j)+head_grid%u_2(i+1,1,j))*head_grid%sina(i,j)
            cff2 = 2.0/(((head_grid%U10(i,j)+head_grid%U10(i+1,j))**2+&
              (head_grid%V10(i,j)+head_grid%V10(i,j+1))**2)**0.5+eps)
            cff3 = 0.5*(head_grid%V10(i,j)+head_grid%V10(i,j+1))*head_grid%cosa(i,j)+&
              0.5*(head_grid%U10(i,j)+head_grid%U10(i+1,j))*head_grid%sina(i,j)
            cff = cff1*cff2*(head_grid%UST(i,j)**2)*cff3	  
            fieldPtr(i,j) = cff
          end do
        end do
      case ("stsh")
        do i = ips, ipe
          do j = jps, jpe	
            fieldPtr(i,j) = head_grid%HFX(i,j)
	  enddo
	enddo
      case ("stlh")
        do i = ips, ipe
          do j = jps, jpe
            fieldPtr(i,j) = head_grid%LH(i,j)
          enddo
        enddo
      case ("rsns")
        do i = ips, ipe
          do j = jps, jpe
            fieldPtr(i,j) = head_grid%GSW(i,j)
          end do
        end do	      
      case ("rlns")
        do i = ips, ipe
          do j = jps, jpe
            fieldPtr(i,j) = - ( head_grid%GLW(i,j) - &
             (STBOLT*head_grid%EMISS(i,j)*head_grid%SST(i,j)**4) )
          end do
        end do
      case ("evap")
        do i = ips, ipe
          do j = jps, jpe
            fieldPtr(i,j) = head_grid%QFX(i,j)
          end do
        end do
      case ("prec")
        do i = ips, ipe
          do j = jps, jpe
	    fieldPtr(i,j) = (head_grid%RAINCV(i,j)+head_grid%RAINNCV(i,j))/&
	      head_grid%DT * 86400.D0
          end do
        end do
      case default
        call ESMF_LogWrite("  ATM unknown field name", &
          ESMF_LOGMSG_INFO, rc=rc) 
        rc = 1 
    end select

  end subroutine


end module
