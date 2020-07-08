!  Copyright (C) 2019, respective authors of MCFM.
!
!  This program is free software: you can redistribute it and/or modify it under
!  the terms of the GNU General Public License as published by the Free Software
!  Foundation, either version 3 of the License, or (at your option) any later
!  version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!  PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program. If not, see <http://www.gnu.org/licenses/>
 
      function scetint(r,wgt)
          use ieee_arithmetic
          use types
          use Scalevar
          use PDFerrors
          use SCET
          use MCFMStorage
          use m_gencuts, only : enable_reweight_user, reweight_user
      implicit none
      real(dp):: scetint
      include 'constants.f'
      include 'mxpart.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'kprocess.f'
      include 'kpart.f'
      include 'energy.f'
      include 'npart.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'taucut.f'
      include 'first.f'
      include 'scalevar.f'
      include 'scale.f'
      include 'facscale.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'x1x2.f'
      include 'xmin.f'
      include 'debug.f'
      integer itrial
      real(dp):: savescale,savefacscale,xmsqvar(6),alphas
      real(dp):: p(mxpart,4),pjet(mxpart,4),r(mxdim),W,xmsq,
     & val,val2,pswt,wgt,z1,z2,flux,BrnRat
      integer j
      logical:: bin,includedipole
      real(dp) :: xjac
      common/bin/bin
      common/BrnRat/BrnRat

      real(dp) :: scet_xmsq

      scetint=0._dp

      W=sqrts**2

      currentPDF = 0

      p(:,:)=0._dp
      pjet(:,:)=0._dp

      call gen_lops(r,p,pswt,*999)
    
      call dotem(npart+2,p,s)

      if (ntau == 0) then
c----reject event if any tau is too small -- important for precision
         call smalltau(p,npart,*999)
      else
! small safety cuts
        call smallnew(p,npart,*999)
      endif

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif

      z1=r(ndim-1)**2
      z2=r(ndim)**2

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

      if ( (xx(1) > one)  .or. (xx(2) > one)
     & .or.(xx(1) < xmin) .or. (xx(2) < xmin)) goto 999

      ! CENTRAL VALUE
      if (dynamicscale) call scaleset(initscale,initfacscale,p)
      xmsq = scet_xmsq(z1,z2,p,.true.)
      xjac = four*sqrt(z1*z2)
      flux = fbGeV2/(two*xx(1)*xx(2)*W)
      scetint = flux*xjac*pswt*xmsq/BrnRat

      call getptildejet(0,pjet)
      val=scetint*wgt 
      val2=val**2 

      if (ieee_is_nan(val) .or. (.not. ieee_is_finite(val))) then
          if (debug) then
              write(6,*) 'Discarded NaN, val=',val
          endif
          goto 999
      endif

      ! SCALE VARIATION WITH CENTRAL PDF
      if (doScalevar .and. bin .and. xmsq /= 0._dp) then
          if (dynamicscale) then
              call scaleset(initscale,initfacscale,p)
          endif

          savescale = scale
          savefacscale = facscale

          do j=1,maxscalevar
              call usescales(savescale*scalevarmult(j),
     &                       savefacscale*facscalevarmult(j))
              scalereweight(j) = scet_xmsq(z1,z2,p,.false.)/xmsq
              
          enddo

          ! restore
          call usescales(savescale, savefacscale)
      endif

      ! PDF VARIATION WITH CENTRAL SCALE

      if (doPDFerrors .and. bin) then
          do j=1,maxPDFsets
              currentPDF = j
              if (doPDFAlphas) then
                  if (dynamicscale) then
                      call scaleset(initscale,initfacscale,p)
                  else
                      call usescales(initscale,savefacscale)
                  endif
                  call updateAlphas(scale)
              endif
              pdfreweight(currentPDF) = (scetint - flux*xjac*pswt*scet_xmsq(z1,z2,p,.false.)/BrnRat)*wgt
          enddo
      endif

      if (bin) then
          includeTaucutgrid(0) = .true.
          call nplotter(pjet,val,val2,0)
      endif

      if (enable_reweight_user) then
          scetint = scetint * reweight_user(pjet)
      endif

      return

 999  continue

      scetint = 0._dp
      end

      function scet_xmsq(z1,z2,p,central)
          use types
          implicit none
          include 'mxpart.f'
          include 'constants.f'
          real(dp) :: scet_xmsq
          real(dp), intent(in) :: z1, z2
          real(dp), intent(in) :: p(mxpart,4)
          logical, intent(in) :: central

          real(dp) :: xjac, QB(2)
          include 'first.f'
          include 'taucut.f'
          include 'kprocess.f'
          include 'kpart.f'
          include 'x1x2.f'
          include 'energy.f'
          integer, save :: iorder
!$omp threadprivate(iorder)

          real(dp) :: xmsq

          xx(1)=-2._dp*p(1,4)/sqrts
          xx(2)=-2._dp*p(2,4)/sqrts

          QB(1)=-two*p(1,4)
          QB(2)=-two*p(2,4)
          
          if (tauboost) then
! prov    ide beam energies in singlet c.o.m. instead if required
            QB(1)=sqrt(QB(1)*QB(2))
            QB(2)=QB(1)
          endif

c--- d    etermine order of calculation on first call
          if (first) then
            first=.false.
            if (kpart==knnlo) then
              iorder=2
            elseif (kpart==ksnlo) then
              iorder=1
            else
              write(6,*) 'Error in scetint: kpart=',kpart
              stop
            endif
          endif

c--- C    alculate the required matrix elements
          if     (kcase==kW_only) then
            call lumxmsq_w(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kZ_only) then
            call lumxmsq_z(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif ((kcase==kggfus0) .or. (kcase==kHigaga)) then
            call lumxmsq_h(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kHi_Zga) then
            call lumxmsq_h_Zga(p,xx,z1,z2,QB,iorder,xmsq)
          elseif ((kcase==kWHbbar) .or. (kcase==kWHgaga)
     &       .or. (kcase==kWH__WW)) then
            call lumxmsq_wh(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif ((kcase==kZHbbar) .or. (kcase==kZHgaga)
     &       .or. (kcase==kZH__WW)) then
            call lumxmsq_zh(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kgamgam) then
             call lumxmsq_gaga(p,xx,z1,z2,QB,iorder,xmsq,central)
!          elseif (kcase==kW_1jet) then
!            call lumxmsq_w1jet(p,xx,z1,z2,QB,iorder,xmsq)
!          elseif (kcase==kZ_1jet) then
!            call lumxmsq_z1jet(p,xx,z1,z2,QB,iorder,xmsq)
!          elseif (kcase==kggfus1) then
!            call lumxmsq_h1jet(p,xx,z1,z2,QB,iorder,xmsq)
          elseif (kcase==kZgamma) then
             call set_anomcoup(p)
             call lumxmsq_zgamma(p,xx,z1,z2,QB,iorder,xmsq,central)
          elseif (kcase==kWgamma) then
             call lumxmsq_wgamma(p,xx,z1,z2,QB,iorder,xmsq,central)
          else
            error stop 'Process not yet available in jettiness formalism'
          endif

          scet_xmsq = xmsq

      end function

