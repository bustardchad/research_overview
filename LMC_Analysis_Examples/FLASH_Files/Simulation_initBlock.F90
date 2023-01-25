!!****if* source/Simulation/SimulationMain/magneto/CRStream/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID, 
!!                       integer(IN) :: myPE)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the soundwave problem in Uhlig 2012
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  myPE -             my processor number
!!
!!***

subroutine Simulation_initBlock(blockID, myPE)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, &
                             Grid_getDeltas, &
                             Grid_getBlkPtr, &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID
  integer, intent(in) :: myPE

  integer :: i, j, k, n, istat, l
  real :: xx, yy, zz, del(MDIM)
  real :: wave
  real :: KuzminPhi, kb, meanmass, a_scale, b_scale
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
          eintZone, enerZone, ekinZone, ecrZone, gamcZone, gameZone, pcrZone
  real,allocatable, dimension(:) :: x, y, z, xL, yL, zL, xR, yR, zR
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Az,Ax,Ay
#else
  real, allocatable, dimension(:,:,:) :: Az !,Ax,Ay
#endif

  real,pointer,dimension(:,:,:,:) :: solnData, facexData, faceyData, facezData
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis
  logical :: gcell = .true.
  real :: rad, M, G, R, a, ai, Na
  real :: vphi, costheta, sintheta, rho20kpc, Rs, xparam, radSphere, radSphere1, radSphere2
  real :: rhoPoint1, rhoPoint2, dphidz1, dphidz2, dpdx, dpdy, dpdr
  real, allocatable, dimension(:) :: zVec
  integer :: cnt, point

  ! for Tonnesen + Stone 2014 B fields:
  real :: a_zf

  real, parameter :: kpc  = 3.0856775807000E+21
  real, parameter :: msun = 1.9889225000000E+33
  real, parameter :: pi = 3.14159265359

  ! dump some output to stdout listing the paramters
  if (myPE == MASTER_PE) then
     
1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
     
  endif
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getDeltas(blockID, del)
  if (NDIM == 3) then
    sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
    allocate(z(sizeZ), zL(sizeZ), zR(sizeZ))
    call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, sizeZ)
    call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, gcell, zL, sizeZ)
    call Grid_getCellCoords(KAXIS, blockId, RIGHT_EDGE, gcell, zR, sizeZ)
  endif
  if (NDIM >= 2) then
    sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
    allocate(y(sizeY), yL(sizeY), yR(sizeY))
    call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, sizeY)
    call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, gcell, yL, sizeY)
    call Grid_getCellCoords(JAXIS, blockId, RIGHT_EDGE, gcell, yR, sizeY)
  endif
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  allocate(x(sizeX), xL(sizeX), xR(sizeX))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xL, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xR, sizeX)

  call Grid_getBlkPtr(blockID,solnData,CENTER)
#if NFACE_VARS > 0 
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     if (NDIM >= 2) call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif
#ifndef FIXEDBLOCKSIZE
  if (NDIM == 2) then
     allocate(Az(sizeX+1,sizeY+1,1),stat=istat)
  elseif (NDIM == 3) then
     allocate(Az(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
  endif
#endif

  allocate(zVec(1020))
!------------------------------------------------------------------------------
 ! print *, "checkpoint 1"
  G = 6.67259E-8
  M = galaxyMass * msun
  !meanmass = 1.021E-24 ! g
  meanmass = 1.67E-24 ! g
  kb = 1.380658E-16 ! Boltzmann constant
  Na = 6.022140E23
  a_scale = (a_kpc)*kpc ! 1.7 kpc, same as in gravity file
  b_scale = (b_kpc)*kpc ! 0.34 kpc, same as in gravity file
  ai = sqrt(sim_gamma_gas*kb*tempGal/meanmass)
  R = galaxyRadius
 ! a = sqrt((ai**2)*rhoGal*exp(-G*M/(2*(ai**2)*R))/rhoR)

! Loop over cells in the block.  
  !populate array of Z coordinates to use for numerical integration of dp/dz
  do cnt = 0, 1000
    ! zVec(cnt+1) = ((10000. - cnt)/10000.)*(5.0*sim_zMax)
     zVec(cnt+1) = ((1000. - cnt)/1000.)*(1.50*sim_zMax)
  enddo
        
 ! print *, "checkpoint 2"
  Rs = 3.0*3.0856E21

!  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)+1
!     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)+1
!        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)+1
!                        
!           xx  = x(i)
!           yy  = y(j)
!           zz  = z(k)
!
!           rad = xx**2 + yy**2 ! + zz**2
!           rad = sqrt(rad)
!           
!           if (fieldLoop == .true.) then ! TOR from Tonnesen + Stone 2014
!                if ((rad < cutoff*kpc) .and. (abs(zz) < 5.0*b_scale)) then
!                        ! Az(i,j,k) = 0.000005*(sim_fieldLoopRadius-rad)
!                       ! a_zf = a_init*(-abs(zz)/b_scale + 1.)**(4.0)
!                        ! Made changes to this March 3, 2019
!
!                        if (useTS2014 == .true.) then
!                               ! a_zf = a_init*(-abs(zz)/(cutoff*kpc) + 1.)**(80.0)
!                                a_zf = a_init*(cosh(zz/b_scale))**(-2.) ! squared because it goes into square root
!                                Az(i,j,k) = sqrt(a_zf)*exp(-6.0*rad/(cutoff*kpc))*(-6.0*sin(2.5*rad/(cutoff*kpc)) - &
!                                        2.5*cos(2.5*rad/(cutoff*kpc)))/(6.0**2.0 + 2.5**2.0)
!                        else
!                               ! Az(i,j,k) = -BToroidal*(cutoff*kpc-rad)*sqrt((-abs(zz)/(cutoff*kpc) + 1.)**(80.0))
!                                Az(i,j,k) = -BToroidal*(cutoff*kpc-rad)*(cosh(zz/b_scale))**(-1.)
!                        end if                        
!                        Ay(i,j,k) = 0.0
!                        Ax(i,j,k) = 0.0
!                else
!                        Az(i,j,k) = 0.0
!                        Ay(i,j,k) = 0.0
!                        Ax(i,j,k) = 0.0
!                end if
!           end if
!        enddo
!     enddo
!  enddo

  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           xx  = x(i)
           yy  = y(j)
           zz  = z(k)

           rad = xx**2 + yy**2 ! + zz**2
           rad = sqrt(rad)
           
           KuzminPhi = G*M/sqrt(rad**2 + (a_scale + sqrt(zz**2 + b_scale**2))**2) - G*M/(a_scale+b_scale)
           if (useRotation == .true.) then
              KuzminPhi = KuzminPhi + 0.5*(vphi**2)
           end if

           if (useHSE == .true.) then
              rhoZone = rhoGal*exp(KuzminPhi/(ai**2)) !only true if isothermal
              presZone = pre_factor*rhoZone*ai**2 ! assuming isothermal
           else if ((useHSE == .false.) .and. (useSechProfile == .true.)) then
              !rhoZone = rhoGal
              rhoZone = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.),cutrho) 
              if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                rho20kpc = 5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh((20.0*kpc)/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.)
                rhoZone = max((rhoZone*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoZone = cutrho
                end if
               ! print *, rhoZone
              end if
              presZone = 0.0
              
              point = 1
              do cnt = 1, 1000
                if ((zVec(cnt+1) < abs(zz)) .and. (zVec(cnt) > abs(zz))) then
                        point = cnt ! point to stop integration
                    !    print *, point
                end if
              enddo
            !  print *, "checkpoint 3, point = ", point
            !  if (rhoZone > cutrho) then
            !    presZone = Na*kb*rhoZone*tempGal ! or using cs^2 = gamma P / rho
            !  else
            !    presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
            !  end if
              do l = 1, point
                rhoPoint1 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l)/b_scale))**(-1.), cutrho)
                rhoPoint2 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l+1)/b_scale))**(-1.), cutrho)
                if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                        rhoPoint1 = max((rhoPoint1*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                        
                        rhoPoint2 = max((rhoPoint2*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                end if
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoPoint1 = cutrho
                        rhoPoint2 = cutrho
                end if
                radSphere1 = sqrt(xx**2 + yy**2 + zVec(l)**2)
                radSphere2 = sqrt(xx**2 + yy**2 + zVec(l+1)**2)
                dphidz1 = G*M*zVec(l)*(1+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere1/Rs) / (radSphere1**3.0) &
                                - 1.0/(Rs*(radSphere1**2.0)*(1.0+radSphere1/Rs))))
                dphidz2 = G*M*zVec(l+1)*(1+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere2/Rs) / (radSphere2**3.0) &
                                - 1.0/(Rs*(radSphere2**2.0)*(1.0+radSphere2/Rs))))
                presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(1.50*sim_zMax/1000.)
               ! dphidz1 = G*M*zVec(l)*(1+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5)
               ! dphidz2 = G*M*zVec(l+1)*(1+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5)
               ! presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(15.0*sim_zMax/10000.)
              enddo
             ! print *, "checkpoint 4"

              if (rhoZone <= cutrho) then
                presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
              end if
           else
               rhoZone  = rhoR
               presZone = presR
           end if

           ecrZone  = ecr0
           pcrZone  = ecrZone * rhoZone * (sim_gamma_cr - 1.)
                
           if ((useRamPressure == .true.) .and. (rhoZone <= cutrho)) then
               ! velxZone = 1.E5*(22.13 - 177.67 + 265.69)
               ! velyZone = 1.E5*(1127.35 + 1728.84 - 840.58 + 218.56)
               ! velzZone = 1.E5*(38.82 - 123.64 + 162.03)
                velxZone = 1.E5*(-187.909 + 545.606 - 565.672 + 220.114)
                velyZone = 1.E5*(75.755 - 178.239 + 158.829)
                velzZone = 1.E5*(-739.545 + 2122.302 - 2119.255 + 764.728 + 192.382)
           else
                velxZone = 0.0
                velyZone = 0.0
                velzZone = 0.0
           end if

           ekinZone = 0.5 * (velxZone**2 + velyZone**2 + velzZone**2)
           eintZone = presZone / (sim_gamma_gas-1.)
           eintZone = eintZone / rhoZone + ecrZone
           enerZone = eintZone + ekinZone
           gamcZone = (sim_gamma_gas*presZone + sim_gamma_cr*pcrZone) / (presZone + pcrZone)
           gameZone = (sim_gamma_gas*(sim_gamma_cr-1.)*presZone + sim_gamma_cr*(sim_gamma_gas-1.)*pcrZone) / &
                      ((sim_gamma_cr-1.)*presZone + (sim_gamma_gas-1.)*pcrZone)


           solnData(STMP_VAR,i,j,k) = presZone/(Na*kb*rhoZone)          
           
           presZone = presZone + pcrZone
           
           solnData(CLOO_VAR,i,j,k) = 0.0
           solnData(DENS_VAR,i,j,k) = rhoZone
           solnData(DENO_VAR,i,j,k) = rhoZone
           solnData(PRES_VAR,i,j,k) = presZone
           Az(i,j,k) = presZone
           solnData(VELX_VAR,i,j,k) = velxZone
           solnData(VELY_VAR,i,j,k) = velyZone
           solnData(VELZ_VAR,i,j,k) = velzZone
           solnData(ENER_VAR,i,j,k) = enerZone
           solnData(EINT_VAR,i,j,k) = eintZone
           solnData(GAMC_VAR,i,j,k) = gamcZone
           solnData(GAME_VAR,i,j,k) = gameZone
           solnData(MASS_SCALARS_BEGIN:MASS_SCALARS_END,i,j,k) = ecrZone
           if (rad .gt. 0.0) then
                   costheta = xx/rad
                   sintheta = yy/rad
           else
                   costheta = 0.0
                   sintheta = 0.0
           end if
          ! if ((fieldLoop == .true.) .and. (rad < cutoff*kpc) .and. (abs(zz) < 5.0*b_scale)) then
           if ((fieldLoop == .true.) .and. (rad < (cutoff+cutlength)*kpc) .and. (abs(zz) < 5.0*b_scale)) then
               ! solnData(MAGX_VAR,i,j,k)= BToroidal*sintheta
               ! solnData(MAGY_VAR,i,j,k)= -BToroidal*costheta

               ! print *, Az(i,j,k), (Az(i,j+1,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i+1,j,k)), del(JAXIS)
               ! solnData(MAGX_VAR,i,j,k)=  .5*3.0856E21*(Az(i,j+1,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i+1,j,k))/del(JAXIS)
               ! solnData(MAGY_VAR,i,j,k)= -.5*3.0856E21*(Az(i+1,j,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i,j+1,k))/del(IAXIS)
                solnData(MAGX_VAR,i,j,k) = a_init*(cosh(zz/b_scale))**(-1.) * exp(-6.0*rad/(cutoff*kpc))*sintheta*sin(2.5*rad/(cutoff*kpc))
                solnData(MAGY_VAR,i,j,k) = -a_init*(cosh(zz/b_scale))**(-1.) * exp(-6.0*rad/(cutoff*kpc))*costheta*sin(2.5*rad/(cutoff*kpc))
               ! print *, solnData(MAGX_VAR,i,j,k), solnData(MAGY_VAR,i,j,k)
              !  solnData(MAGY_VAR,i,j,k)= -.5*3.0856E21*(Az(i+1,j,k)-Az(i,j,k) + Az(i+1,j+1,k)-Az(i,j+1,k))/del(IAXIS)
                solnData(MAGZ_VAR,i,j,k) = 0.
                solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                     solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
               ! print *, solnData(MAGX_VAR,i,j,k)
                solnData(DIVB_VAR,i,j,k) = 0.
           else ! no field loop
                solnData(MAGX_VAR,i,j,k)= BHalo
                solnData(MAGY_VAR,i,j,k)= 0.0
                solnData(MAGZ_VAR,i,j,k)= 0.0
           end if
           if (rhoZone > 1.E-26) then  ! for painting cells, in the ISM
                solnData(IGM_SPEC,i,j,k) = sim_smallFrac
                solnData(ISM_SPEC,i,j,k) = 1.0-sim_smallFrac
                solnData(MTL_SPEC,i,j,k) = 0.3  ! metallicity
           else ! in the halo
                solnData(IGM_SPEC,i,j,k) = 1.0-sim_smallFrac
                solnData(ISM_SPEC,i,j,k) = sim_smallFrac
                solnData(MTL_SPEC,i,j,k) = 0.01
           end if
          ! xx  = x(i+1)
           xx  = x(i) + del(IAXIS)
           yy  = y(j)
           zz  = z(k)

          ! print *, "checkpoint 5, xx = ", xx
           rad = xx**2 + yy**2 ! + zz**2
           rad = sqrt(rad)
          

           KuzminPhi = G*M/sqrt(rad**2 + (a_scale + sqrt(zz**2 + b_scale**2))**2) - G*M/(a_scale+b_scale)
           if (useRotation == .true.) then
              KuzminPhi = KuzminPhi + 0.5*(vphi**2)
           end if

           if (useHSE == .true.) then
              rhoZone = rhoGal*exp(KuzminPhi/(ai**2)) !only true if isothermal
              presZone = pre_factor*rhoZone*ai**2 ! assuming isothermal
           else if ((useHSE == .false.) .and. (useSechProfile == .true.)) then
              !rhoZone = rhoGal
              rhoZone = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.),cutrho) 
              if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                rho20kpc = 5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh((20.0*kpc)/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.)
                rhoZone = max((rhoZone*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoZone = cutrho
                end if
               ! print *, rhoZone
              end if
              presZone = 0.0
              
           !   if (rhoZone > cutrho) then
           !     presZone = Na*kb*rhoZone*tempGal ! or using cs^2 = gamma P / rho
           !   else
           !     presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
           !   end if
              point = 1
              do cnt = 1, 1000
                if ((zVec(cnt+1) < abs(zz)) .and. (zVec(cnt) > abs(zz))) then
                        point = cnt ! point to stop integration
                    !    print *, point
                end if
              enddo
              do l = 1, point
                rhoPoint1 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l)/b_scale))**(-1.), cutrho)
                rhoPoint2 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l+1)/b_scale))**(-1.), cutrho)
                if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                        rhoPoint1 = max((rhoPoint1*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                        rhoPoint2 = max((rhoPoint2*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                end if
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoPoint1 = cutrho
                        rhoPoint2 = cutrho
                end if
                radSphere1 = sqrt(xx**2 + yy**2 + zVec(l)**2)
                radSphere2 = sqrt(xx**2 + yy**2 + zVec(l+1)**2)
                dphidz1 = G*M*zVec(l)*(1.0+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere1/Rs) / (radSphere1**3.0) &
                                - 1.0/(Rs*(radSphere1**2.0)*(1.0+radSphere1/Rs))))
                dphidz2 = G*M*zVec(l+1)*(1.0+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere2/Rs) / (radSphere2**3.0) &
                                - 1.0/(Rs*(radSphere2**2.0)*(1.0+radSphere2/Rs))))
                presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(1.50*sim_zMax/1000.)
               ! dphidz1 = G*M*zVec(l)*(1+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5) &
               !         + 
               ! dphidz2 = G*M*zVec(l+1)*(1+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5)
               ! presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(15.0*sim_zMax/10000.)
              enddo
              if (rhoZone <= cutrho) then
                presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
              end if

           else
               rhoZone  = rhoR
               presZone = presR
           end if

           ecrZone  = ecr0
           pcrZone  = ecrZone * rhoZone * (sim_gamma_cr - 1.)

           if ((useRamPressure == .true.) .and. (rhoZone <= cutrho)) then
               ! velxZone = 1.E5*(22.13 - 177.67 + 265.69)
               ! velyZone = 1.E5*(1127.35 + 1728.84 - 840.58 + 218.56)
               ! velzZone = 1.E5*(38.82 - 123.64 + 162.03)
                velxZone = 1.E5*(-187.909 + 545.606 - 565.672 + 220.114)
                velyZone = 1.E5*(75.755 - 178.239 + 158.829)
                velzZone = 1.E5*(-739.545 + 2122.302 - 2119.255 + 764.728 + 192.382)
           else
                velxZone = 0.0
                velyZone = 0.0
                velzZone = 0.0
           end if
        
          ! velxZone = 0.0
          ! velyZone = 0.0
          ! velzZone = 0.0

           ekinZone = 0.5 * (velxZone**2 + velyZone**2 + velzZone**2)
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone + ecrZone
           enerZone = eintZone + ekinZone
           gamcZone = (sim_gamma_gas*presZone + sim_gamma_cr*pcrZone) / (presZone + pcrZone)
           gameZone = (sim_gamma_gas*(sim_gamma_cr-1.)*presZone + sim_gamma_cr*(sim_gamma_gas-1.)*pcrZone) / &
                      ((sim_gamma_cr-1.)*presZone + (sim_gamma_gas-1.)*pcrZone)

           presZone = presZone + pcrZone
           
          ! solnData(PRES_VAR,i+1,j,k) = presZone
           Az(i+1,j,k) = presZone
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! xx  = x(i-1)
           xx  = x(i) - del(IAXIS)
           yy  = y(j)
           zz  = z(k)

          ! print *, "checkpoint 6, xx = ", xx
           rad = xx**2 + yy**2 ! + zz**2
           rad = sqrt(rad)
          
           KuzminPhi = G*M/sqrt(rad**2 + (a_scale + sqrt(zz**2 + b_scale**2))**2) - G*M/(a_scale+b_scale)
           if (useRotation == .true.) then
              KuzminPhi = KuzminPhi + 0.5*(vphi**2)
           end if

           if (useHSE == .true.) then
              rhoZone = rhoGal*exp(KuzminPhi/(ai**2)) !only true if isothermal
              presZone = pre_factor*rhoZone*ai**2 ! assuming isothermal
           else if ((useHSE == .false.) .and. (useSechProfile == .true.)) then
              rhoZone = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.),cutrho) 
              if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                rho20kpc = 5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh((20.0*kpc)/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.)
                rhoZone = max((rhoZone*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoZone = cutrho
                end if
              end if
              presZone = 0.0
              
              point = 1
              do cnt = 1, 1000
                if ((zVec(cnt+1) < abs(zz)) .and. (zVec(cnt) > abs(zz))) then
                        point = cnt ! point to stop integration
                end if
              enddo
              do l = 1, point
                rhoPoint1 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l)/b_scale))**(-1.), cutrho)
                rhoPoint2 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l+1)/b_scale))**(-1.), cutrho)
                if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                        rhoPoint1 = max((rhoPoint1*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                        rhoPoint2 = max((rhoPoint2*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                end if
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoPoint1 = cutrho
                        rhoPoint2 = cutrho
                end if
                radSphere1 = sqrt(xx**2 + yy**2 + zVec(l)**2)
                radSphere2 = sqrt(xx**2 + yy**2 + zVec(l+1)**2)
                dphidz1 = G*M*zVec(l)*(1.0+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere1/Rs) / (radSphere1**3.0) &
                                - 1.0/(Rs*(radSphere1**2.0)*(1.0+radSphere1/Rs))))
                dphidz2 = G*M*zVec(l+1)*(1.0+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere2/Rs) / (radSphere2**3.0) &
                                - 1.0/(Rs*(radSphere2**2.0)*(1.0+radSphere2/Rs))))
                presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(1.50*sim_zMax/1000.)
               ! dphidz1 = G*M*zVec(l)*(1+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5)
               ! dphidz2 = G*M*zVec(l+1)*(1+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5)
               ! presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(15.0*sim_zMax/10000.)
              enddo
              
              if (rhoZone <= cutrho) then
                presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
              end if

           else
               rhoZone  = rhoR
               presZone = presR
           end if

           ecrZone  = ecr0
           pcrZone  = ecrZone * rhoZone * (sim_gamma_cr - 1.)

           if ((useRamPressure == .true.) .and. (rhoZone <= cutrho)) then
               ! velxZone = 1.E5*(22.13 - 177.67 + 265.69)
               ! velyZone = 1.E5*(1127.35 + 1728.84 - 840.58 + 218.56)
               ! velzZone = 1.E5*(38.82 - 123.64 + 162.03)
                velxZone = 1.E5*(-187.909 + 545.606 - 565.672 + 220.114)
                velyZone = 1.E5*(75.755 - 178.239 + 158.829)
                velzZone = 1.E5*(-739.545 + 2122.302 - 2119.255 + 764.728 + 192.382)
           else
                velxZone = 0.0
                velyZone = 0.0
                velzZone = 0.0
           end if
          ! velxZone = 0.0
          ! velyZone = 0.0
          ! velzZone = 0.0

           ekinZone = 0.5 * (velxZone**2 + velyZone**2 + velzZone**2)
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone + ecrZone
           enerZone = eintZone + ekinZone
           gamcZone = (sim_gamma_gas*presZone + sim_gamma_cr*pcrZone) / (presZone + pcrZone)
           gameZone = (sim_gamma_gas*(sim_gamma_cr-1.)*presZone + sim_gamma_cr*(sim_gamma_gas-1.)*pcrZone) / &
                      ((sim_gamma_cr-1.)*presZone + (sim_gamma_gas-1.)*pcrZone)

           presZone = presZone + pcrZone
           
          ! solnData(PRES_VAR,i-1,j,k) = presZone

           if (i .gt. 1) then
                Az(i-1,j,k) = presZone
           end if
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

                 
           xx  = x(i)
          ! yy  = y(j+1)
           yy  = y(j) + del(JAXIS)
           zz  = z(k)

          ! print *, "checkpoint 7, yy = ", yy
           rad = xx**2 + yy**2 ! + zz**2
           rad = sqrt(rad)
          
           KuzminPhi = G*M/sqrt(rad**2 + (a_scale + sqrt(zz**2 + b_scale**2))**2) - G*M/(a_scale+b_scale)
           if (useRotation == .true.) then
              KuzminPhi = KuzminPhi + 0.5*(vphi**2)
           end if

           if (useHSE == .true.) then
              rhoZone = rhoGal*exp(KuzminPhi/(ai**2)) !only true if isothermal
              presZone = pre_factor*rhoZone*ai**2 ! assuming isothermal
           else if ((useHSE == .false.) .and. (useSechProfile == .true.)) then
              rhoZone = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.),cutrho) 
              if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                rho20kpc = 5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh((20.0*kpc)/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.)
                rhoZone = max((rhoZone*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoZone = cutrho
                end if
              end if
              presZone = 0.0
             ! if (rhoZone > cutrho) then
             !   presZone = Na*kb*rhoZone*tempGal ! or using cs^2 = gamma P / rho
             ! else
             !   presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
             ! end if
              
              point = 1
              do cnt = 1, 1000
                if ((zVec(cnt+1) < abs(zz)) .and. (zVec(cnt) > abs(zz))) then
                        point = cnt ! point to stop integration
                    !    print *, point
                end if
              enddo
              do l = 1, point
                rhoPoint1 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l)/b_scale))**(-1.), cutrho)
                rhoPoint2 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l+1)/b_scale))**(-1.), cutrho)
                if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                        rhoPoint1 = max((rhoPoint1*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                        rhoPoint2 = max((rhoPoint2*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                end if
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoPoint1 = cutrho
                        rhoPoint2 = cutrho
                end if
                radSphere1 = sqrt(xx**2 + yy**2 + zVec(l)**2)
                radSphere2 = sqrt(xx**2 + yy**2 + zVec(l+1)**2)
                dphidz1 = G*M*zVec(l)*(1+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere1/Rs) / (radSphere1**3.0) &
                                - 1.0/(Rs*(radSphere1**2.0)*(1.0+radSphere1/Rs))))
                dphidz2 = G*M*zVec(l+1)*(1+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere2/Rs) / (radSphere2**3.0) &
                                - 1.0/(Rs*(radSphere2**2.0)*(1.0+radSphere2/Rs))))
                presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(1.50*sim_zMax/1000.)
              enddo
              if (rhoZone <= cutrho) then
                presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
              end if

           else
               rhoZone  = rhoR
               presZone = presR
           end if

           ecrZone  = ecr0
           pcrZone  = ecrZone * rhoZone * (sim_gamma_cr - 1.)

        
           if ((useRamPressure == .true.) .and. (rhoZone <= cutrho)) then
               ! velxZone = 1.E5*(22.13 - 177.67 + 265.69)
               ! velyZone = 1.E5*(1127.35 + 1728.84 - 840.58 + 218.56)
               ! velzZone = 1.E5*(38.82 - 123.64 + 162.03)
                velxZone = 1.E5*(-187.909 + 545.606 - 565.672 + 220.114)
                velyZone = 1.E5*(75.755 - 178.239 + 158.829)
                velzZone = 1.E5*(-739.545 + 2122.302 - 2119.255 + 764.728 + 192.382)
           else
                velxZone = 0.0
                velyZone = 0.0
                velzZone = 0.0
           end if
          ! velxZone = 0.0
          ! velyZone = 0.0
          ! velzZone = 0.0

           ekinZone = 0.5 * (velxZone**2 + velyZone**2 + velzZone**2)
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone + ecrZone
           enerZone = eintZone + ekinZone
           gamcZone = (sim_gamma_gas*presZone + sim_gamma_cr*pcrZone) / (presZone + pcrZone)
           gameZone = (sim_gamma_gas*(sim_gamma_cr-1.)*presZone + sim_gamma_cr*(sim_gamma_gas-1.)*pcrZone) / &
                      ((sim_gamma_cr-1.)*presZone + (sim_gamma_gas-1.)*pcrZone)

           presZone = presZone + pcrZone
           
          ! solnData(PRES_VAR,i,j+1,k) = presZone
           Az(i,j+1,k) = presZone
           !!!!!!!!!!!!!!!!!!!!!!!!!


           xx  = x(i)
           !yy  = y(j-1)
           yy  = y(j) - del(JAXIS)
           zz  = z(k)

          ! print *, "checkpoint 8, yy = ", yy
           rad = xx**2 + yy**2 ! + zz**2
           rad = sqrt(rad)
          
           KuzminPhi = G*M/sqrt(rad**2 + (a_scale + sqrt(zz**2 + b_scale**2))**2) - G*M/(a_scale+b_scale)
           if (useRotation == .true.) then
              KuzminPhi = KuzminPhi + 0.5*(vphi**2)
           end if

           if (useHSE == .true.) then
              rhoZone = rhoGal*exp(KuzminPhi/(ai**2)) !only true if isothermal
              presZone = pre_factor*rhoZone*ai**2 ! assuming isothermal
           else if ((useHSE == .false.) .and. (useSechProfile == .true.)) then
              !rhoZone = rhoGal
              rhoZone = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.),cutrho) 
              if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                rho20kpc = 5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh((20.0*kpc)/a_scale))**(-1.) * (cosh(zz/b_scale))**(-1.)
                rhoZone = max((rhoZone*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoZone = cutrho
                end if
               ! print *, rhoZone
              end if
              presZone = 0.0
              
           !   if (rhoZone > cutrho) then
           !     presZone = Na*kb*rhoZone*tempGal ! or using cs^2 = gamma P / rho
           !   else
           !     presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
           !   end if
              point = 1
              do cnt = 1, 1000
                if ((zVec(cnt+1) < abs(zz)) .and. (zVec(cnt) > abs(zz))) then
                        point = cnt ! point to stop integration
                    !    print *, point
                end if
              enddo
              do l = 1, point
                rhoPoint1 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l)/b_scale))**(-1.), cutrho)
                rhoPoint2 = max(5.e8*msun*(0.5**2)/(2*pi*(a_scale**2)*b_scale) &
                        * (cosh(rad/a_scale))**(-1.) * (cosh(zVec(l+1)/b_scale))**(-1.), cutrho)
                if (rad > cutoff*kpc) then !truncating smoothly (Roedigger and Bruggen)
                        rhoPoint1 = max((rhoPoint1*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                        rhoPoint2 = max((rhoPoint2*0.5*(1.0+cos(pi*(rad-(cutoff*kpc))/(cutlength*kpc)))), cutrho)
                end if
                if (rad > (cutoff*kpc + cutlength*kpc)) then
                        rhoPoint1 = cutrho
                        rhoPoint2 = cutrho
                end if
                radSphere1 = sqrt(xx**2 + yy**2 + zVec(l)**2)
                radSphere2 = sqrt(xx**2 + yy**2 + zVec(l+1)**2)
                dphidz1 = G*M*zVec(l)*(1+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere1/Rs) / (radSphere1**3.0) &
                                - 1.0/(Rs*(radSphere1**2.0)*(1.0+radSphere1/Rs))))
                dphidz2 = G*M*zVec(l+1)*(1+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
                        /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5) + &
                        (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+radSphere2/Rs) / (radSphere2**3.0) &
                                - 1.0/(Rs*(radSphere2**2.0)*(1.0+radSphere2/Rs))))
                presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(1.50*sim_zMax/1000.)
               ! dphidz1 = G*M*zVec(l)*(1.0+a_scale/sqrt(zVec(l)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l)**2 + b_scale**2))**2)**(1.5) + &
               !         (zVec(l)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+xparam) / (radSphere**3.0) & 
               !                 - 1.0/(Rs*(radSphere**2.0)*(1.0+xparam))))
               ! dphidz2 = G*M*zVec(l+1)*(1.0+a_scale/sqrt(zVec(l+1)**2 + b_scale**2)) & 
               !         /(rad**2 + (a_scale + sqrt(zVec(l+1)**2 + b_scale**2))**2)**(1.5) + &
               !         (zVec(l+1)*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+xparam) /(radSphere**3.0) &
               !                 - 1.0/(Rs*(radSphere**2.0)*(1.0+xparam))))
               ! presZone = presZone + (rhoPoint1*dphidz1 + rhoPoint2*dphidz2)*0.5*(1.50*sim_zMax/1000.)
              enddo

              if (rhoZone <= cutrho) then
                presZone = Na*kb*rhoZone*1.E6  ! halo temp = 1e6 K
              end if
           else
               rhoZone  = rhoR
               presZone = presR
           end if

           ecrZone  = ecr0
           pcrZone  = ecrZone * rhoZone * (sim_gamma_cr - 1.)

           if ((useRamPressure == .true.) .and. (rhoZone <= cutrho)) then
               ! velxZone = 1.E5*(22.13 - 177.67 + 265.69)
               ! velyZone = 1.E5*(1127.35 + 1728.84 - 840.58 + 218.56)
               ! velzZone = 1.E5*(38.82 - 123.64 + 162.03)
                velxZone = 1.E5*(-187.909 + 545.606 - 565.672 + 220.114)
                velyZone = 1.E5*(75.755 - 178.239 + 158.829)
                velzZone = 1.E5*(-739.545 + 2122.302 - 2119.255 + 764.728 + 192.382)
           else
                velxZone = 0.0
                velyZone = 0.0
                velzZone = 0.0
           end if
          ! velxZone = 0.0
          ! velyZone = 0.0
          ! velzZone = 0.0

           ekinZone = 0.5 * (velxZone**2 + velyZone**2 + velzZone**2)
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone + ecrZone
           enerZone = eintZone + ekinZone
           gamcZone = (sim_gamma_gas*presZone + sim_gamma_cr*pcrZone) / (presZone + pcrZone)
           gameZone = (sim_gamma_gas*(sim_gamma_cr-1.)*presZone + sim_gamma_cr*(sim_gamma_gas-1.)*pcrZone) / &
                      ((sim_gamma_cr-1.)*presZone + (sim_gamma_gas-1.)*pcrZone)

           presZone = presZone + pcrZone
           
          ! solnData(PRES_VAR,i,j-1,k) = presZone
           if (j .gt. 1) then
                Az(i,j-1,k) = presZone
           end if

 
 !       enddo
 !    enddo
 ! enddo


          ! print *, "checkpoint 9, starting to get rotation vphi"
           if (useRotation == .true.) then
                if (Az(i+1,j,k) > 0.0) then
                        dpdx = (Az(i+1,j,k) - Az(i,j,k))/(del(IAXIS))
                else
                        dpdx = (Az(i,j,k) - Az(i-1,j,k))/(del(IAXIS))
                end if
                if (Az(i,j+1,k) > 0.0) then
                        dpdy = (Az(i,j+1,k) - Az(i,j,k))/(del(JAXIS))
                else  
                        dpdy = (Az(i,j,k) - Az(i,j-1,k))/(del(JAXIS))
                end if
                dpdr = sqrt(dpdx**2 + dpdy**2)
                xx  = x(i)
                yy  = y(j)
                zz  = z(k)

                rad = xx**2 + yy**2 ! + zz**2
                rad = sqrt(rad)
                Rs = 3.0*3.0856E21
                radSphere = sqrt(xx**2 + yy**2 + zz**2)
                xparam = radSphere/Rs
               ! if (rad <= cutoff*kpc) then
                if (rad <= (cutoff+cutlength)*kpc) then
                        vphi = sqrt(abs(((rad/solnData(DENS_VAR,i,j,k)) * abs(dpdr)) - &
                                (rad*(G*M*rad/((rad**2.0 + (a_scale+sqrt((b_scale**2.) & 
                                + zz**2.0))**2.)**(1.5))) + rad*(rad*G*4.0*pi*3.0E-24*(Rs**3.0)*(log(1.0+xparam) / (radSphere**3.0) & 
                                - 1.0/(Rs*(radSphere**2.0)*(1.0+xparam)))))))
                else
                        vphi = 0.0
                end if
                if (rad .gt. 0.0) then
                        costheta = xx/rad
                        sintheta = yy/rad
                else
                        costheta = 0.0
                        sintheta = 0.0
                end if


                if ((useRamPressure == .true.) .and. (solnData(DENS_VAR,i,j,k) <= cutrho)) then
                       ! velxZone = 1.E5*(22.13 - 177.67 + 265.69)
                       ! velyZone = 1.E5*(1127.35 + 1728.84 - 840.58 + 218.56)
                       ! velzZone = 1.E5*(38.82 - 123.64 + 162.03)
                        velxZone = 1.E5*(-187.909 + 545.606 - 565.672 + 220.114)
                        velyZone = 1.E5*(75.755 - 178.239 + 158.829)
                        velzZone = 1.E5*(-739.545 + 2122.302 - 2119.255 + 764.728 + 192.382)
                else
                        velxZone = vphi*sintheta
                        velyZone = -vphi*costheta
                        velzZone = 0.0e0
                end if
                solnData(VPHI_VAR,i,j,k) = vphi
                solnData(VELX_VAR,i,j,k) = velxZone
                solnData(VELY_VAR,i,j,k) = velyZone
                solnData(VELZ_VAR,i,j,k) = velzZone

           end if
        enddo
       ! print *, "finished first do loop"
     enddo
    ! print *, "finished second do loop"
  enddo      

!  print *, "end of first loops"

! Initialize face variables for MHD
#ifdef MAGX_VAR
#if NFACE_VARS > 0

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
  do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
  do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
    rad = sqrt(x(i)**2.0 + y(j)**2.0)
    if (rad .gt. 0.0) then
            costheta = x(i)/rad
            sintheta = y(j)/rad
    else
            costheta = 0.0
            sintheta = 0.0
    end if
    ! Insert uniform B field 
    if(sim_killdivb) then
      facexData(MAG_FACE_VAR,i,j,k) = BHalo
      if(NDIM >= 2) faceyData(MAG_FACE_VAR,i,j,k) = 0.0
      if(NDIM == 3) facezData(MAG_FACE_VAR,i,j,k) = 0.0
     ! if ((fieldLoop == .true.) .and. (rad < cutoff*kpc) .and. (abs(z(k)) < 5.0*b_scale)) then
      if ((fieldLoop == .true.) .and. (rad < (cutoff+cutlength)*kpc) .and. (abs(z(k)) < 5.0*b_scale)) then
                facexData(MAG_FACE_VAR,i,j,k)= a_init*(cosh(z(k)/b_scale))**(-1.) * exp(-6.0*rad/(cutoff*kpc))*sintheta*sin(2.5*rad/(cutoff*kpc))
                faceyData(MAG_FACE_VAR,i,j,k)= -a_init*(cosh(z(k)/b_scale))**(-1.) * exp(-6.0*rad/(cutoff*kpc))*costheta*sin(2.5*rad/(cutoff*kpc))
               ! facexData(MAG_FACE_VAR,i,j,k)= -3.0856E21*(Ay(i,j,k+1)-Ay(i,j,k))/del(KAXIS) + (Az(i,j+1,k)-Az(i,j,k))/del(JAXIS)
               ! faceyData(MAG_FACE_VAR,i,j,k)=  3.0856E21*(Ax(i,j,k+1)-Ax(i,j,k))/del(KAXIS) - (Az(i+1,j,k)-Az(i,j,k))/del(IAXIS)
               ! facezData(MAG_FACE_VAR,i,j,k)= -3.0856E21*(Ax(i,j+1,k)-Ax(i,j,k))/del(JAXIS) + (Ay(i+1,j,k)-Ay(i,j,k))/del(IAXIS)
                facezData(MAG_FACE_VAR,i,j,k)= 0.0
      endif

    endif

  enddo !i
  enddo !j
  enddo !k
#endif
#endif

 ! print *, "end of second loops"
 ! do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
 ! do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
 ! do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
  do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
  do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
#if NFACE_VARS > 0
    solnData(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i,j,k)+facexData(MAG_FACE_VAR,i+1,j,k))
    if(NDIM >= 2) solnData(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j,k)+faceyData(MAG_FACE_VAR,i,j+1,k))
    if(NDIM == 3) solnData(MAGZ_VAR,i,j,k) = 0.5*(facezData(MAG_FACE_VAR,i,j,k)+facezData(MAG_FACE_VAR,i,j,k+1))

#if NDIM == 1
    solnData(DIVB_VAR,i,j,k) = 0.
#elif NDIM >= 2
    solnData(DIVB_VAR,i,j,k)= &
      (facexData(MAG_FACE_VAR,i+1,j,  k  ) - facexData(MAG_FACE_VAR,i,j,k))/del(IAXIS) &
    + (faceyData(MAG_FACE_VAR,i,  j+1,k  ) - faceyData(MAG_FACE_VAR,i,j,k))/del(JAXIS)
#if NDIM == 3
    solnData(DIVB_VAR,i,j,k)= solnData(DIVB_VAR,i,j,k) &
    + (facezData(MAG_FACE_VAR,i,  j,  k+1) - facezData(MAG_FACE_VAR,i,j,k))/del(KAXIS)
#endif
#endif

#else !NFACE_VARS == 0
    solnData(DIVB_VAR,i,j,k) = 0.
#endif !NFACE_VARS

    ! Update the magnetic pressure
    solnData(MAGP_VAR,i,j,k) = .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                              solnData(MAGX_VAR:MAGZ_VAR,i,j,k))
  enddo !i
  enddo !j
  enddo !k

 ! print *, "end of third loops"
  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     if (NDIM >= 2) call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  deallocate(x)
  deallocate(xL)
  deallocate(xR)
  if (NDIM >= 2) then
    deallocate(y)
    deallocate(yL)
    deallocate(yR)
  endif
  if (NDIM == 3) then
    deallocate(z)
    deallocate(zL)
    deallocate(zR)
  endif
 
  deallocate(zVec)

 ! print *, "released blk pointer and deallocated everything"
  return
end subroutine Simulation_initBlock










