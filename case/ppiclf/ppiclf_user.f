!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rpi,rmu,rhof,rmass,vmag,dp,rep,ap,volp,  
     >       psip, cds, fqsx, fqsy, fqsz,fbx,fby,fbz,
     >       grad, fmagz, fsafz, flift, rvx, rvy, rvz, zloc

      integer*4 i
!
      rpi  = 4.0*atan(1.0)
      rmu  = 1.8E-5
      rhof = 1.225
      
! evaluate ydot
      do i=1,ppiclf_npart
         
         ! particle mass
         
         dp    = ppiclf_rprop(PPICLF_R_JDP,i)
         ap = rpi*(dp/2.0D0)**2
         volp = rpi/6.0D0*dp**3
         psip = ppiclf_rprop(PPICLF_R_JPSIP,i)
         
         rmass = volp
     >          *(ppiclf_rprop(PPICLF_R_JRHOP,i) - rhof)



         ! Dioguardia 2018

         rvx = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
         rvy = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
         rvz = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
         
         vmag = sqrt(rvx**2 + rvy**2 + rvz**2)

         rep   = vmag*dp*rhof/rmu
         
         cds =    24/rep*((1-psip)/rep + 1)**0.25
     >          + 24/rep*(0.1806*rep**0.6459)*psip**(-(rep**0.08))
     >          + 0.4251/(1 + 6880.95/rep*psip**(5.05)) 

         fqsx = 0.5*rhof*cds*ap*vmag*rvx
         fqsy = 0.5*rhof*cds*ap*vmag*rvy
         fqsz = 0.5*rhof*cds*ap*vmag*rvz

         
         ! Gravity
         fbx  = 0.0
         fby  = 0.0
         fbz  = -9.8*rmass

         ! NEW LIFTING MODEL
         
         grad = ppiclf_rprop(PPICLF_R_JVXZ,i)
         zloc = ppiclf_y(PPICLF_JZ,i)

         fmagz = rpi/8.0D0*rhof*vmag*dp**3*(0.5*abs(grad))
         flift = 190.0*fmagz*exp(-zloc**2/(2*0.4247**2)) 
         !190 - 0.4247 (50% attenuation at 0.5 m)



         ! set ydot for all PPICLF_SLN number of equations

         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fbx+fqsx)/rmass
         ppiclf_ydot(PPICLF_JVY,i) = (fby+fqsy)/rmass
         ppiclf_ydot(PPICLF_JVZ,i) = (fbz+fqsz+flift)/rmass


      enddo 
! evaluate ydot

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
!
! Input:
!
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
!
! Output:
!
      real*8 map  (PPICLF_LRP_PRO)
!
! Internal:
!
      real*8 dp_norm
!

      ! particle volume divided by particle diameter for 2d
!      dp_norm = 1./rprop(PPICLF_R_JDP)
!      dp_norm = 1.0
!      map(PPICLF_P_JPHIP) = dp_norm*rprop(PPICLF_R_JVOLP)
!      map(PPICLF_P_JFX)   = dp_norm*ydotc(PPICLF_JVX)
!      map(PPICLF_P_JFY)   = dp_norm*ydotc(PPICLF_JVY)
!      map(PPICLF_P_JFZ)   =         ydotc(PPICLF_JVZ) 
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
      
      
      
      
      
!----------------------------------------
! THIS IS NOT BEING CALLED
!----------------------------------------






      integer*4 i
      integer*4 j
      real*8 yi    (PPICLF_LRS)
      real*8 rpropi(PPICLF_LRP)
      real*8 yj    (PPICLF_LRS)
      real*8 rpropj(PPICLF_LRP)
!
! Internal:
!
      real*8 ksp,erest
      common /ucollision/ ksp,erest

      real*8 rpi2, rthresh, rxdiff,rad, rydiff, rzdiff, rdiff, rm1, rm2,
     >       rmult, eta, rbot, rn_12x, rn_12y, rn_12z, rdelta12,
     >       rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra

     
! boundaries
      if (j .eq. 0) then

         ! give a bit larger collision threshold for walls
         rextra   = 0.0d0
         rthresh  = 5.0*rpropi(PPICLF_R_JDP)

         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = sqrt(rxdiff**2 + rydiff**2 + rzdiff**2)
         ppiclf_rprop(PPICLF_R_JVXZ,i) = rdiff	 
         if (rdiff .gt. rthresh) return
                 
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
                  
         rv12_mag = -1.0d0*(yi(PPICLF_JVX)*rn_12x +
     >                      yi(PPICLF_JVY)*rn_12y +
     >                      yi(PPICLF_JVZ)*rn_12z)


!         if (rv12_mag .lt. 0.0) then
         rnmag     = 2.0*rv12_mag

!         if (yi(PPICLF_JZ) .lt. 0.5*rpropi(PPICLF_R_JDP)) then 
!         ppiclf_y(PPICLF_JZ,i) = 0.5d0
!         endif
         
!         ppiclf_y(PPICLF_JVX,i) = yi(PPICLF_JVX)
!     >                              + rnmag*rn_12x
!         ppiclf_y(PPICLF_JVY,i) = yi(PPICLF_JVY)
!     >                              + rnmag*rn_12y
!         ppiclf_y(PPICLF_JVZ,i) = yi(PPICLF_JVZ)
!     >                              + rnmag*rn_12z
!         endif

      endif

      return
      end
!-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine ppiclf_user_PostRK3
!
! Ground collision (thesis Eqs. 5.9-5.12). Called from the one-line
! hook added to ppiclf_solve_IntegrateRK3 by patches/.
!
      implicit none
!
      include "PPICLF"
!
      integer*4 j, k, ip, ndum, dp_loc
      real*8 dp_p, z_p, vz_p, t_frc, t_imp, vz_imp
!
      ip = 0
      ndum = PPICLF_NPART*PPICLF_LRS
      do j=3,ndum-3,PPICLF_LRS
         dp_loc = PPICLF_R_JDP-3+j + ip*(PPICLF_LRP-PPICLF_LRS)
         dp_p = 0.5*ppiclf_rprop(dp_loc,1)
         ip = ip + 1
         k = j + 3
         z_p = ppiclf_y(j,1)
         vz_p = ppiclf_y(k,1)
         if ((z_p.lt.dp_p).and.(vz_p.le.0.0)) then
            if (abs(z_p - ppiclf_y1(j)) .gt. 1.0e-6) then
               t_frc = (dp_p-ppiclf_y1(j))/(z_p-ppiclf_y1(j))
            else
               t_frc = 0.0
            endif
            t_imp = ppiclf_dt*t_frc
            vz_imp = ppiclf_y1(k) + t_frc*(vz_p - ppiclf_y1(k))
            ppiclf_y(k,1) = -1.0*vz_imp
            ppiclf_y(j,1) = dp_p + ppiclf_y(k,1)*(ppiclf_dt-t_imp)
         endif
      enddo
      return
      end
