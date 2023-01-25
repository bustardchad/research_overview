!!****if* source/Simulation/SimulationMain/WindTunnel/Grid_applyBCEdge
!!  Updated by Chad -- May 22, 2019 to correct errors with ram pressure
!!      Also correcting initial velocity when useRamPressure == .true. in
!!      Simulation_initBlock.F90
!! NAME
!!  Grid_applyBCEdge
!!
!! SYNOPSIS
!!
!!  Grid_applyBCEdge(integer(IN)              :: bcType, 
!!                   integer(IN)              :: bcDir,
!!                   integer(IN)              :: guard,
!!                   integer(IN)              :: var,
!!                   real(INOUT),dimension(:) :: dataRow(2*guard),
!!                   integer(IN)              :: face,
!!                   integer(IN)              :: gridDataStruct,
!!                   integer(IN),OPTIONAL     :: blockHandle,
!!                   real(in),OPTIONAL    :: secondCoord,
!!                   real(in),OPTIONAL    :: thirdCoord)
!!  
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Applies the boundary conditions to a given vector.
!!  This routine applies the boundary conditions on a given face (lowerface
!!  or upperface) of a given vector. 
!!     If (face=LOW)dataRow(1:guard) = boundary values
!!     If (face=HIGH) dataRow(guard+1:2*guard) = boundary values
!!  The reason why information about direction and variable is included in
!!  this interface is because velocities need to be treated specially
!!  for REFLECTING boundary conditions. 
!!  This implementation is specific to the WindTunnel problem which 
!!  requires inflow and outflow boundary conditions.
!!
!!
!!  
!! ARGUMENTS 
!!
!!
!!  bcType -   the type of boundary condition being applied to this face
!!              -  USER_DEFINED: This routine does its own thing.
!!              -  other types:  This routine reproduces the actions of the Grid unit's
!!                 own implementation (at least for cases that may occur in the
!!                 WindTunnel simulation.)     
!!  bcDir -    can take on values IAXIS,JAXIS or KAXIS. This is needed
!!             for handling the reflective boundary conditions. If bcDir=IAXIS,
!!             and boundary conditions are reflective, and X velocity is
!!             treated differently from all other variables. similarly if bcDir
!!             is JAXIS, then Y velocity is different.
!!  guard -    number of guardcells 
!!  var   -    The variable on which boundary conditions are applied
!!             It is used with bcDir for reflective boundary conditions
!!             to correctly handle velocities. This argument is redundant 
!!             for all other variables.
!!  dataRow -  storage for the data being operated upon.
!!  face    -  can take values LOW and HIGH, defined in constants.h
!!             to indicate whether to apply boundary on lowerface or 
!!             upperface
!!  gridDataStruct : integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (unk or work for PM) (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!
!!   blockHandle - the identity of the block under consideration
!!  secondCoord,thirdCoord - scalar coordinate values in the coordinate
!!                         directions perpendicular to the sweep direction.
!!                         This is not needed for simple boundary condition types
!!                         such as REFLECTIVE or OUTFLOW, but is provided so that
!!                         more complex boundary conditions can make use of it.
!!  NOTES 
!!            This routine exists in the simulation directory of 
!!            the WindTunnel problem, and therefore replaces
!!            the default implementation in the Grid unit.
!!            
!!
!!***


subroutine Grid_applyBCEdge(bcType,bcDir,guard,var,dataRow,face,&
     gridDataStruct, blockHandle, secondCoord, thirdCoord)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getCellCoords, Grid_putPointData
  use Driver_interface !, ONLY : Driver_abortFlash
!  use Driver_data, ONLY :dr_myPE                   ! only needed for print* in error case
  implicit none

# include "constants.h"
# include "Flash.h"

  integer, intent(in):: bcType
  integer,intent(IN) :: var,guard,face,bcDir,gridDataStruct
  real,dimension(:),intent(INOUT)::dataRow
  integer,intent(IN),OPTIONAL :: blockHandle
  real,intent(IN),OPTIONAL :: secondCoord,thirdCoord
  real :: kine,jcoord,delyCoord, error,kcoord,delzcoord,error2
  integer :: i     !loop counter
  integer :: k, sign, sizeY,istat,sizeZ
  real,pointer :: solnData(:,:,:,:)
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  logical :: gcell = .true.
  real,allocatable,dimension(:) :: yCoord,yCoordf,zCoord,zCoordf
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real :: meanmass, kb, ai, sim_pInflow, radSalemkpc, rhoSalem
  real :: simTime, simTimeGyr, rhoIn, velx,vely,velz,presZone,pcrZone
  meanmass = 1.67E-24 ! g
  kb = 1.380658E-16 ! Boltzmann constant
  ai = sqrt(6.022141E23*kb*tempGal)


!  if (var.eq.MAG_FACE_VAR.and. &
!       (gridDataStruct.eq.FACEX.or.&
!        gridDataStruct.eq.FACEY.or.&
!        gridDataStruct.eq.FACEZ)) then
!     dataRow(i)=0.
!  end if
!  jcoord=secondCoord
!  jcoord =(jcoord *64) +256 
!  size=64
call Grid_getBlkIndexLimits(blockHandle,blkLimits,blkLimitsGC)
call Grid_getBlkPtr(blockHandle, solnData, CENTER) !use blockhandle instead
call Grid_getBlkPtr(blockHandle,facexData,FACEX)
call Grid_getBlkPtr(blockHandle,faceyData,FACEY)
call Grid_getBlkPtr(blockHandle,facezData,FACEZ)

sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1

allocate(yCoordf(sizeY+1),stat=istat)
call Grid_getCellCoords(JAXIS, blockHandle, FACES, gcell, yCoordf, sizeY+1)
allocate(yCoord(sizeY),stat=istat)
call Grid_getCellCoords(JAXIS, blockHandle, CENTER, gcell, ycoord, sizeY)!GW

allocate(zCoordf(sizeZ+1),stat=istat)
call Grid_getCellCoords(KAXIS, blockHandle, FACES, gcell, zCoordf, sizeZ+1)
allocate(zCoord(sizeZ),stat=istat)
call Grid_getCellCoords(KAXIS, blockHandle, CENTER, gcell, zcoord, sizeZ)!GW

!print*,'blkLimits(LOW,IAXIS) blkLimits(HIGH,IAXIS)',blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS) 
!print*,'blkLimitsGC(LOW,IAXIS) blkLimitsGC(HIGH,IAXIS)',blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS) 

call Driver_getSimTime(simTime)
simTimeGyr = simTime/(3.155E16)
radSalemkpc = -75.06*(1.-simTimeGyr)**2 + 263.61*(1.-simTimeGyr) + 60.1
!rhoSalem = 2.0*10**(-24)*0.46*(1.0+(radSalemkpc/0.35)**2.0)**(-3.0*0.559/2.0)
rhoSalem = (2.E-24)*0.46*(1.0+(radSalemkpc/0.35)**2.0)**(-3.0*0.559/2.0)
!rhoIn = rhoIGM + 1.e-28*(-0.65*(2.-simTimeGyr)**3 + 3.33*(2.-simTimeGyr)**2 - 6*(2.-simTimeGyr) + 4.54)
rhoIn = rhoIGM + 1.e-28*(-5.42*(1.-simTimeGyr)**3 + 10.86*(1.-simTimeGyr)**2 - 8.56*(1.-simTimeGyr) + 3.97)
if (SalemRam == .true.) then
       ! rhoIn = rhoIGM + 1.e-28*(-7.79*(1.-simTimeGyr)**3 + 12.89*(1.-simTimeGyr)**2 - 7.31*(1.-simTimeGyr) + 1.73)
       ! rhoIn = rhoIGM + rhoSalem
        rhoIn = rhoSalem  ! Changed Nov 29 -- Chad
endif
sim_pInflow = rhoIn*ai*ai  ! not the problem
presZone = sim_pInflow
pcrZone = 0.0

!print *, "rhoIn: ", rhoIn
!print *, "rhoSalem ", rhoSalem
!print *, radSalemkpc
!sim_pInflow = rhoIn*ai*ai
!if ((gridDataStruct.eq.FACEX.or.&
!        gridDataStruct.eq.FACEY.or.&
!        gridDataStruct.eq.FACEZ)) then
            
          !  call Grid_getBlkPtr(blockHandle,facexData,FACEX)
          !  call Grid_getBlkPtr(blockHandle,faceyData,FACEY)
 if (gridDataStruct.eq.FACEY) then  !var.eq.MAG_FACE_VAR.or.var.eq.MAGI_FACE_VAR gw         
          delyCoord=yCoordf(blkLimitsGC(LOW,JAXIS)+1)- yCoordf(blkLimitsGC(LOW, JAXIS))
          !print*,"delycoord",delycoord,yCoordf(blkLimitsGC(LOW,JAXIS)+1)
          jcoord=secondCoord -yCoordf(blkLimitsGC(LOW, JAXIS))
          jcoord=jcoord/delycoord +blkLimitsGC(LOW,JAXIS)  !jugad GW
          error= secondCoord - yCoordf(jcoord)
          !print*,"face",jcoord,error 
else
         
         delyCoord=yCoord(blkLimitsGC(LOW,JAXIS)+1)- ycoord(blkLimitsGC(LOW, JAXIS))
         jcoord=secondCoord -ycoord(blkLimitsGC(LOW, JAXIS))
         jcoord=jcoord/delycoord +blkLimitsGC(LOW,JAXIS)
         error= secondCoord - yCoord(jcoord)
 !        print*,"center",jcoord,error  
endif

if (gridDataStruct.eq.FACEZ) then !var.eq.MAG_FACE_VAR.or.var.eq.MAGI_FACE_VAR gw         
          delzCoord=zCoordf(blkLimitsGC(LOW,KAXIS)+1)- zCoordf(blkLimitsGC(LOW, KAXIS))
         !print*,"delycoord",delycoord,yCoordf(blkLimitsGC(LOW,JAXIS)+1)
          kcoord=thirdCoord -zCoordf(blkLimitsGC(LOW, KAXIS))
          kcoord=kcoord/delzcoord +blkLimitsGC(LOW,KAXIS)  !jugad GW
          error2= thirdCoord - zCoordf(kcoord)
         !print*,"face",jcoord,error 
else
        
         delzCoord=zCoord(blkLimitsGC(LOW,KAXIS)+1)- zcoord(blkLimitsGC(LOW, KAXIS))
         kcoord=thirdCoord -zcoord(blkLimitsGC(LOW, KAXIS))
         kcoord=kcoord/delzcoord + blkLimitsGC(LOW,KAXIS)
         error2= thirdCoord - zCoord(kcoord)
  !      print*,"center",kcoord,error2  
endif
!kcoord = blkLimitsGC(LOW, KAXIS)

!print*,'jcoord',jcoord
!print*,'secondCoord ycalculated',secondCoord,yCoord(jcoord)
  sign = 1
  if (gridDataStruct==CENTER) then
#ifdef VELX_VAR
     if ((bcDir==IAXIS).and.(var==VELX_VAR))sign=-1
#endif
#ifdef VELY_VAR
     if((bcDir==JAXIS).and.(var==VELY_VAR))sign=-1
#endif
#ifdef VELZ_VAR
     if((bcDir==KAXIS).and.(var==VELZ_VAR))sign=-1
#endif
  end if

  if((bcType)==PERIODIC) return

  select case (gridDataStruct)
  case(CENTER)
     if(face==LOW) then
        select case (bcType)
        case(REFLECTING)
           k = 2*guard+1
           do i = 1,guard
              dataRow(i)= dataRow(k-i)*sign
           end do
        case(OUTFLOW)
           do i = 1,guard
              dataRow(i)= dataRow(guard+1)
           end do
        case(USER_DEFINED) !use this for ram pressure inflow
         !  print *, "using User Defined BC"
         !  do i = 1,guard !set outflow as default
         !     dataRow(i)= dataRow(guard+1)
         !  end do
           select case(var)
           case(GAMC_VAR)
              dataRow(1:guard)=(sim_gamma_gas*presZone + sim_gamma_cr*pcrZone) / &
                  (presZone + pcrZone)
           case(GAME_VAR)
              dataRow(1:guard)=(sim_gamma_gas*(sim_gamma_cr-1.)*presZone + & 
                  sim_gamma_cr*(sim_gamma_gas-1.)*pcrZone) / &
                  ((sim_gamma_cr-1.)*presZone + (sim_gamma_gas-1.)*pcrZone)
         !  case(DENS_VAR) !not using rhoUp in this version, just outflow
         !     do i = 1,guard
         !       dataRow(i)=dataRow(guard+1)
         !     end do
           case(DENS_VAR) 
             ! dataRow(1:guard)=rhoIGM
              dataRow(1:guard)=rhoIn
           case(PRES_VAR) !
              dataRow(1:guard)=sim_pInflow
           case(VELX_VAR)
             ! dataRow(1:guard)=velUp*(1.-exp(-simTime/tauRam))
              !velx = (1.-exp(-simTime/tauRam))*1.E5*(30.03*(2.-simTimeGyr)**2 - 148.68*(2.-simTimeGyr) + 328.26)
             ! velx = (1.-exp(-simTime/tauRam))*1.E5*(22.13*(1.-simTimeGyr)**2.0 - 177.67*(1.-simTimeGyr) + 265.69)
              velx = (1.-exp(-simTime/tauRam))*1.E5*(-187.909*(1.0-simTimeGyr)**3.0 + 545.606*(1.0-simTimeGyr)**2.0 - 565.672*(1.0-simTimeGyr) + 220.114)
              dataRow(1:guard)=velx
           case(VELY_VAR)
             ! vely = (1.-exp(-simTime/tauRam))*1.E5*(1127.35*(1.-simTimeGyr)**3.0 + 1728.84*(1.-simTimeGyr)**2.0 - 840.58*(1.-simTimeGyr) + 218.56)
              vely = (1.-exp(-simTime/tauRam))*1.E5*(75.755*(1.0-simTimeGyr)**2.0 - 178.239*(1.0-simTimeGyr) + 158.829)
           !   dataRow(1:guard)=vely
              dataRow(1:guard)=vely
           case(VELZ_VAR)
             ! velz = (1.-exp(-simTime/tauRam))*velUp
             ! velz = (1.-exp(-simTime/tauRam))*1.E5*(38.82*(1.-simTimeGyr)**2.0 - 123.64*(1.-simTimeGyr) + 162.03)
              velz = (1.-exp(-simTime/tauRam))*1.E5*(-739.545*(1.0-simTimeGyr)**4.0 + 2122.302*(1.0-simTimeGyr)**3.0 - 2119.255*(1.0-simTimeGyr)**2.0 + 764.728*(1.0-simTimeGyr) + 192.382)  
              dataRow(1:guard)=velz
           case(MAGX_VAR)
             ! dataRow(1:guard)=0.0
              dataRow(1:guard)= BHalo ! Chad -- Feb 14, 2019
           case(MAGY_VAR)
              dataRow(1:guard)=0.0
           case(MAGZ_VAR)
              dataRow(1:guard)=0.0
           case(ENER_VAR)
              kine = 0.5 *(velx**2 + vely**2 + velz**2) !vortex frame GW
              dataRow(1:guard) = max( sim_pInflow / ((sim_gamma_gas-1.) * rhoIn) + kine, 1.E-60)
           case(EINT_VAR)
              dataRow(1:guard) = max( sim_pInflow / ((sim_gamma_gas-1.) * rhoIn), 1.E-60)
           case(SPECIES_BEGIN) ! igm? 
              dataRow(1:guard)=1.0e0-sim_smallFrac
           case(SPECIES_BEGIN+1) ! ism?
              dataRow(1:guard)=sim_smallFrac
           case(SPECIES_BEGIN+2)  ! mtl?
              dataRow(1:guard)= 0.01
           end select
        case default
           call Driver_abortFlash("unsupported boundary condition on Lower Face")
        end select
        
     else
        
        select case (bcType)
        case(REFLECTING)
           k = 2*guard+1
           do i = 1,guard
              dataRow(k-i)= dataRow(i)*sign
           end do
           
        case(OUTFLOW)
           do i = 1,guard
              dataRow(guard+i)= dataRow(guard)
           end do
        case(USER_DEFINED)
          ! print*,'boundary is',bcType,dr_myPE,face
           !call Driver_abortFlash("Simulation WindTunnel does not support USER_DEFINED boundary on Upper Face")
           !solnData(MAGX_VAR,i,j,k) = bx_zone
           !solnData(MAGY_VAR,i,j,k) = by_zone
           !solnData(MAGZ_VAR,i,j,k) = bz_zone
        !   do i = 1,guard
        !      dataRow(guard+i)= dataRow(guard)
        !   end do

           select case(var)
           case(GAMC_VAR)
              dataRow(guard+1:2*guard)=solnData(GAMC_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!sim_gamma
           case(GAME_VAR)
              dataRow(guard+1:2*guard)=solnData(GAME_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!sim_gamma
           case(DENS_VAR)
 !             print *,"printing density"
 !             print *,solnData(DENS_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)
 !             do i = 1,guard
 !                if (solnData(DENS_VAR,blkLimits(LOW,IAXIS)+i-1,jcoord,kcoord).gt.0.0) then
               dataRow(guard+1:2*guard)=solnData(DENS_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!sim_rhoAmbient removing GC1/23/2017
 !                else
                    
           case(PRES_VAR)
              dataRow(guard+1:2*guard)=solnData(PRES_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!sim_pAmbient
           case(VELX_VAR)
              dataRow(guard+1:2*guard)=solnData(VELX_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!-1.8536 !vortex frame gw
           case(VELY_VAR)
              dataRow(guard+1:2*guard)=solnData(VELY_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!0.0
           case(VELZ_VAR)
              dataRow(guard+1:2*guard)=solnData(VELZ_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!0.0
           case(MAGX_VAR)
              !call grid get cell coords from posn. 
              !dataRow(guard+1:2*guard)=solnData(MAGX_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!dataRow(guard-128) !periodic mag field gw
              !dataRow(guard+1:2*guard)=solnData(MAGX_VAR,blkLimitsGC(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)
              dataRow(guard+1:2*guard)=solnData(MAGX_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord) !clarify with SH  
           case(MAGY_VAR)
              dataRow(guard+1:2*guard)=solnData(MAGY_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!dataRow(guard-128) !periodic mag field gw
             !dataRow(guard+1:2*guard)=solnData(MAGX_VAR,blkLimitsGC(LOW,IAXIS):blkLimitsGC(LOW,IAXIS)+guard-1,jcoord,kcoord) 
           case(MAGZ_VAR)
              dataRow(guard+1:2*guard)=solnData(MAGZ_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!0.0
           case(ENER_VAR)
              kine = 0.5 * ((-1.8536)**2) !vortex frame GW
              dataRow(guard+1:2*guard)=solnData(ENER_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!max( sim_pAmbient / ((sim_gamma-1.) * sim_rhoAmbient) + kine, sim_smallP)
           case(EINT_VAR)
              dataRow(guard+1:2*guard)=solnData(EINT_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!max( sim_pAmbient / ((sim_gamma-1.) * sim_rhoAmbient), sim_smallP)
           case(SPECIES_BEGIN)
              do i = 1,guard
                 dataRow(guard+i) = 0.0
              end do
!              dataRow(guard+1:2*guard)=solnData(SPECIES_BEGIN,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!1.0e0-sim_smallX
!           case(SPECIES_BEGIN+1)
!              dataRow(guard+1:2*guard)=solnData(SPECIES_BEGIN+1,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!sim_smallX
!            print *,bcType,dr_myPE
           end select
        case default
!           print*,'boundary is',bcType,dr_myPE
           call Driver_abortFlash("unsupported boundary condition on Upper Face")
        end select
     end if
  case(FACEX)
     if(face==LOW) then !new condition 
       if (bcType.eq.USER_DEFINED) then
          if (var.eq.MAGI_FACE_VAR) then
          ! facexData(MAG_FACE_VAR,i,j,k)=bx_face ! switching blkLimitsGC to blkLimits
          ! dataRow(guard+1:2*guard+1)=facexData(MAGI_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard,jcoord,kcoord)!0.0!dataRow(guard-128)!periodic mag field gw
           dataRow(1:guard)= BHalo
          endif       
          if (var.eq.MAG_FACE_VAR) then
           !facexData(MAG_FACE_VAR,i,j,k)=bx_face  switching blkLimitsGC to blkLimits
          ! dataRow(guard+1:2*guard+1)=facexData(MAG_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard,jcoord,kcoord)!0.0!dataRow(guard-128)!periodic mag field gw
           dataRow(1:guard)= BHalo
          endif
       endif
     endif
  case(FACEY)
    if(face==LOW) then !new condition
     if (bcType.eq.USER_DEFINED) then
        if (var.eq.MAGI_FACE_VAR) then
           dataRow(1:guard)=0.0
        endif
        if (var.eq.MAG_FACE_VAR) then
           dataRow(1:guard)=0.0
        endif
     endif
    endif
  case(FACEZ)
    if(face==LOW) then !new condition
     if (bcType.eq.USER_DEFINED) then
        if (var.eq.MAGI_FACE_VAR) then
           dataRow(1:guard)= 0.0
        endif
        if (var.eq.MAG_FACE_VAR) then
           dataRow(1:guard)= 0.0
        endif
     endif
    endif
!  case(FACEY)
!    if(face==HIGH) then !new condition
!     if (bcType.eq.USER_DEFINED) then
!        if (var.eq.MAGI_FACE_VAR) then
!           !faceyData(MAG_FACE_VAR,i,j,k)=by_face
!           dataRow(guard+1:2*guard)=faceyData(MAGI_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!0.0!dataRow(guard-128)!periodic mag field gw
!          ! dataRow(guard+1:2*guard+1)=faceyData(MAGI_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard,jcoord,kcoord)!0.0!dataRow(guard-128)!periodic mag field gw
!        endif
!        if (var.eq.MAG_FACE_VAR) then
!           !faceyData(MAG_FACE_VAR,i,j,k)=by_face
!           dataRow(guard+1:2*guard)=faceyData(MAG_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!0.0!dataRow(guard-128)!periodic mag field gw
!          ! dataRow(guard+1:2*guard+1)=faceyData(MAG_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard,jcoord,kcoord)!0.0!dataRow(guard-128)!periodic mag field gw
!        endif
!     endif
!    endif 
!  case(FACEZ)
!    if(face==HIGH) then !new condition
!     if (bcType.eq.USER_DEFINED) then
!        if (var.eq.MAGI_FACE_VAR) then
!           dataRow(guard+1:2*guard)=facezData(MAGI_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!0.
!          ! dataRow(guard+1:2*guard+1)=facezData(MAGI_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard,jcoord,kcoord)!0.
!        endif
!        if (var.eq.MAG_FACE_VAR) then
!           dataRow(guard+1:2*guard)=facezData(MAG_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard-1,jcoord,kcoord)!0.
!           !dataRow(guard+1:2*guard+1)=facezData(MAG_FACE_VAR,blkLimits(LOW,IAXIS):blkLimits(LOW,IAXIS)+guard,jcoord,kcoord)!0.
!        endif
!     endif
!    endif 
  end select
  call Grid_releaseBlkPtr(blockHandle,solnData,CENTER)
  call Grid_releaseBlkPtr(blockHandle,facexData,FACEX)
  call Grid_releaseBlkPtr(blockHandle,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockHandle,facezData,FACEZ)
deallocate(yCoord)
deallocate(yCoordf)
deallocate(zCoord)
deallocate(zCoordf)
  return
end subroutine Grid_applyBCEdge
