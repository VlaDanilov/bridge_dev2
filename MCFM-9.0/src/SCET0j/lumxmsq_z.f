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
 
      subroutine lumxmsq_z(p,xx,z1,z2,QB,order,xmsq,central)
          use SCET
      implicit none
      include 'types.f'
c----Matrix element for W production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)-->W^-(e^-(p3)+nbar(p4))
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'scet_const.f'
      include 'taucut.f'
      integer:: j,k,ih1,ih2,m,n,order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),fac,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     & msq(-nf:nf,-nf:nf),assemble,msqlo
      real(dp):: Qh, dot, getxmsq_z
      complex(dp):: prop,qqb,qbq
      common/density/ih1,ih2
      include 'cplx.h'

      logical, intent(in) :: central
      real(dp) :: origtaucut

      fac=4._dp*esq**2*xn

      call spinoru(4,p,za,zb)

c--   calculate propagator
      fac=aveqq*fac/s(3,4)**2
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
c---case dbar-u or ubar-d
      qqb=za(2,3)*zb(4,1)
      qbq=za(1,3)*zb(4,2)

      call softqqbis(order,soft1,soft2)
      call hardqq(s(1,2),musq,hard)

      if (order >= 0) then
      call fdist(ih1,xx(1),facscale,beama0)
      call fdist(ih2,xx(2),facscale,beamb0)
      endif
      if (order >= 1) then
      call xbeam1bis(ih1,z1,xx(1),QB(1),beama1)
      call xbeam1bis(ih2,z2,xx(2),QB(2),beamb1)
      endif
      if (order >= 2) then
      call xbeam2bis(ih1,z1,xx(1),QB(1),beama2)
      call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2)
      endif

      xmsq = getxmsq_z(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qqb,qbq,prop)

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              origtaucut = taucut
              do m=1,size(tcutarray)
                  taucut = tcutarray(m)
                  scetreweight(m) = getxmsq_z(p,xx,order,soft1,soft2,hard,
     &                beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qqb,qbq,prop)
              enddo
              taucut = origtaucut
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif


      return
      end


      function getxmsq_z(p,xx,order,soft1,soft2,hard,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qqb,qbq,prop)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'scet_const.f'
      include 'taucut.f'
      real(dp):: getxmsq_z
      integer:: j,k,order
      real(dp):: p(mxpart,4),fac,
     & xx(2),soft1(-1:1),soft2(-1:3),hard(2),
     & beama0(-5:5),beamb0(-5:5),
     & beama1(-5:5,-1:1),beamb1(-5:5,-1:1),
     & beama2(-5:5,-1:3),beamb2(-5:5,-1:3),
     & z1,z2,QB(2),lum0,lum1(-1:1),lum2(-1:3),bit,
     & msq(-nf:nf,-nf:nf),assemble,msqlo, tauc, getdynamictau
      real(dp):: Qh, dot, powc(-5:5,-5:5)
      complex(dp):: prop,qqb,qbq

      if (dynamictau) then
        tauc=getdynamictau(p)
      else
        tauc=taucut
      endif

! compute power corrections if required
      if ((incpowcorr) .or. (onlypowcorr)) then
        Qh=sqrt(two*dot(p,1,2))
        call powcorr_qa(order,tauc,xx(1),xx(2),Qh,beama0,beamb0,powc)
      endif

      getxmsq_z=zip
      do j=-nf,nf
      k=-j
      if (j == 0) cycle

      if (onlypowcorr) then
        bit=zip
      else
        bit=assemble(order,tauc,
     &   beama0(j),beamb0(k),beama1(j,:),beamb1(k,:),
     &   beama2(j,:),beamb2(k,:),soft1,soft2,hard)
      endif

      if (j > 0) then
        msqlo=fac*(
     &     abs((Q(j)*q1+L(j)*l1*prop)*qqb)**2
     &    +abs((Q(j)*q1+R(j)*r1*prop)*qqb)**2
     &    +abs((Q(j)*q1+L(j)*r1*prop)*qbq)**2
     &    +abs((Q(j)*q1+R(j)*l1*prop)*qbq)**2)
      elseif (j < 0) then
        msqlo=fac*(
     &     abs((Q(k)*q1+L(k)*l1*prop)*qbq)**2
     &    +abs((Q(k)*q1+R(k)*r1*prop)*qbq)**2
     &    +abs((Q(k)*q1+L(k)*r1*prop)*qqb)**2
     &    +abs((Q(k)*q1+R(k)*l1*prop)*qqb)**2)
      else
        msqlo=zip
      endif

      if ((incpowcorr) .or. (onlypowcorr)) then
        bit=bit+powc(j,k)+powc(j,0)+powc(0,k)
      endif

      getxmsq_z=getxmsq_z+bit*msqlo

      enddo

      return
      end
      
