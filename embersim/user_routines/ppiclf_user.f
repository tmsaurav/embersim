!-----------------------------------------------------------------------
      subroutine ppiclf_user_SetYdot
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 rpi, rmu, rhof, rmass, vmag, dp, rep, rphip, rphif,
     >       beta, reps, cds, fqsx, fqsy, fqsz, fbx, fby, fbz, fcx, fcy,
     >       fcz, beta1, beta2, rs,rdum,rhop,taup,fsafx,fsafy,fsafz,
     >       mag_curluc,a1,a2,a3,b1,b2,b3,Rew,Re,Betasaf,alpha,f,Cl,Cld

      integer*4 i
!
      rpi  = 4.0*atan(1.0)
      rmu  = 0.005555555
      rhof = 1.00

! evaluate ydot
      do i=1,ppiclf_npart
         ! particle mass
         rmass = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! Stokes drag
         vmag  = sqrt((ppiclf_rprop(PPICLF_R_JUX,i)
     >                 -ppiclf_y(PPICLF_JVX,i))**2
     >               +(ppiclf_rprop(PPICLF_R_JUY,i)
     >                 -ppiclf_y(PPICLF_JVY,i))**2
     >               +(ppiclf_rprop(PPICLF_R_JUZ,i)
     >                 -ppiclf_y(PPICLF_JVZ,i))**2)

         dp    = ppiclf_rprop(PPICLF_R_JDP,i)
         rep   = vmag*dp*rhof/rmu

         rhop  = ppiclf_rprop(PPICLF_R_JRHOP,i)


         taup = rhop*(dp**2)/(18.0*rmu)

         fqsx = rmass/taup*(1.0 + 1.0/6.0*(rep**(2.0/3.0)))
     >          *(ppiclf_rprop(PPICLF_R_JUX,i)-ppiclf_y(PPICLF_JVX,i))

         fqsy = rmass/taup*(1.0 + 1.0/6.0*(rep**(2.0/3.0)))
     >          *(ppiclf_rprop(PPICLF_R_JUY,i)-ppiclf_y(PPICLF_JVY,i))

         fqsz = rmass/taup*(1.0 + 1.0/6.0*(rep**(2.0/3.0)))
     >          *(ppiclf_rprop(PPICLF_R_JUZ,i)-ppiclf_y(PPICLF_JVZ,i))


!       Begin - Saffman Lift
         mag_curluc = sqrt((ppiclf_rprop(PPICLF_R_JVZY,i)
     >                     -ppiclf_rprop(PPICLF_R_JVYZ,i))**2
     >                     +(ppiclf_rprop(PPICLF_R_JVXZ,i)
     >                     -ppiclf_rprop(PPICLF_R_JVZX,i))**2
     >                     +(ppiclf_rprop(PPICLF_R_JVYX,i)
     >                     -ppiclf_rprop(PPICLF_R_JVXY,i))**2)


         Rew        =    rhof*mag_curluc*(dp**(2))/(rmu)
         Re         =    dp*vmag/rmu
         Betasaf    =    0.5*(Rew/(Re+ 1E-12))
         alpha      =    0.3314*sqrt(Betasaf)
         f          =    (1.0 - alpha)*exp(-0.1*Re) + alpha

         if (Re .lt. 40) then
                Cld = 6.46*f
         else
                Cld = 6.46*0.0524*sqrt(Betasaf*Re)
         endif

         Cl         =    3.0/(2.0*(3.14159265)*sqrt(Rew + 1E-12))*Cld

         a1    = ppiclf_rprop(PPICLF_R_JUX,i)   - ppiclf_y(PPICLF_JVX,i)
         a2    = ppiclf_rprop(PPICLF_R_JUY,i)   - ppiclf_y(PPICLF_JVY,i)
         a3    = ppiclf_rprop(PPICLF_R_JUZ,i)   - ppiclf_y(PPICLF_JVZ,i)

         b1    = ppiclf_rprop(PPICLF_R_JVZY,i)
     >                        - ppiclf_rprop(PPICLF_R_JVYZ,i)
         b2    = ppiclf_rprop(PPICLF_R_JVXZ,i)
     >                        - ppiclf_rprop(PPICLF_R_JVZX,i)
         b3    = ppiclf_rprop(PPICLF_R_JVYX,i)
     >                        - ppiclf_rprop(PPICLF_R_JVXY,i)



         fsafx      = 0.0!   rmass/rhop*rhof*Cl*(a2*b3 - a3*b2)
         fsafy      = 0.0!   rmass/rhop*rhof*Cl*(a3*b1 - a1*b3)
         fsafz      = 0.0!   rmass/rhop*rhof*Cl*(a1*b2 - a2*b1)

         ! Gravity
         fbx  = 0.0
         fby  = 0.0
         fbz  = -9.8*rmass !0.0

         fcx  = 0.0
         fcy  = 0.0
         fcz  = 0.0

         ! set ydot for all PPICLF_SLN number of equations
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ ,i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fbx+fqsx+fcx+fsafx)/rmass
         ppiclf_ydot(PPICLF_JVY,i) = (fby+fqsy+fcy+fsafy)/rmass
         ppiclf_ydot(PPICLF_JVZ,i) = (fbz+fqsz+fcz+fsafz)/rmass

         ppiclf_ydotc(PPICLF_JVX,i) = -fqsx
         ppiclf_ydotc(PPICLF_JVY,i) = -fqsy
         ppiclf_ydotc(PPICLF_JVZ,i) = -fqsz

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
      real*8 dp_norm

      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)
      map(PPICLF_P_JFX)   = ydotc(PPICLF_JVX)
      map(PPICLF_P_JFY)   = ydotc(PPICLF_JVY)
      map(PPICLF_P_JFZ)   = ydotc(PPICLF_JVZ)
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
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
  !
      rpi2  =  9.869604401089358d0

      ! boundaries
      if (j .eq. 0) then

         rad     = ppiclf_rprop(PPICLF_R_JDP,i) /2

	if (yi(PPICLF_JZ) .lt. 0) then
        	rydiff = (yj(PPICLF_JZ)+rad) - (ppiclf_y(PPICLF_JZ,i))
        		if (rydiff .gt. 0) then
				ppiclf_y(PPICLF_JZ,i)  = rad
     & + 0.15*yi(PPICLF_JVX)*yi(PPICLF_JVX) !0.0 ! ppiclf_y(PPICLF_JY,i)
			endif
	endif 


c         if (yi(PPICLF_JZ) .gt. 0) then
c                rydiff = yj(PPICLF_JZ) - (ppiclf_y(PPICLF_JZ,i) + rad)
c                if (rydiff .lt. 0) then
c                    ppiclf_y(PPICLF_JVZ,i) = -ppiclf_y(PPICLF_JVZ,i)
c                    ppiclf_y(PPICLF_JZ,i)  = rad !0.0 !ppiclf_y(PPICLF_JY,i)
c     &   + 2*rydiff
c                endif
c         endif

      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_lagrangian
!
      implicit none
!
      include "PPICLF"

      integer*4 i
      character*8 fmt
      character*7 x1

      do i=1,ppiclf_npart

         if (ppiclf_iprop(6,i) .gt. 1) then
!                   print*, "ppiclf_nid:",ppiclf_nid
!                   print*, "i:", i
!                   print*, "i:", ppiclf_iprop(6,i)
                   fmt = '(I7.7)'
                   write (x1,fmt) ppiclf_iprop(6,i)
                   open(unit=59,file='lagrangian/lagrangian_vel'
     >                  //trim(x1)//'.dat',
     >                  position='append')
                   write(59,*)
     &             ppiclf_y(PPICLF_JX,i),
     &             ppiclf_y(PPICLF_JY,i),
     &             ppiclf_y(PPICLF_JZ,i),
     &             ppiclf_y(PPICLF_JVX,i),
     &             ppiclf_y(PPICLF_JVY,i),
     &             ppiclf_y(PPICLF_JVZ,i),
     &             ppiclf_rprop(PPICLF_R_JUX,i),
     &             ppiclf_rprop(PPICLF_R_JUY,i),
     &             ppiclf_rprop(PPICLF_R_JUZ,i)
                   close(59)
         endif
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_user_fbforce
!
      implicit none
!
      include "PPICLF"

      integer*4 i
      character*8 fmt
      character*7 x1
      real*8      fb_sum, rmass, fb_sum_total, ppiclf_glsum

      external ppiclf_glsum

      rmass = ppiclf_rprop(PPICLF_R_JVOLP,1)
     >          *ppiclf_rprop(PPICLF_R_JRHOP,1)

!      parz = ppiclf_y(PPICLF_JZ,i)

!      print*, "rmass:", parz

      fb_sum    = 0.0
      do i=1,ppiclf_npart
        fb_sum  = fb_sum + ppiclf_ydot(PPICLF_JVX,i)
      enddo

      fb_sum    = fb_sum*rmass

      fb_sum_total = ppiclf_glsum(fb_sum,1)

      if (ppiclf_nid .eq. 0) then
      open(unit=59,file='fb_forcets.dat',position='append')
      write(59,*)
     &         fb_sum_total
      close(59)

      endif

      return
      end
!-----------------------------------------------------------------------
