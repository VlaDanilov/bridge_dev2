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
 
      subroutine dk2qqb_QQb_gs(p,msqc)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis                                               *
*     November, 2011.                                                  *
*     calculate the subtraction term for radiation in the              *
*     anti-top quark decay for the process                             *
*                                                                      *
*     q(-p1)+qbar(-p2) = nu(p3)+e+(p4)+b(p5)                           *
*                        +bbar(p6)+e-(p7)+nubar(p8)+g(p9)              *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained. B-quark is taken to be either massless or massive  *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      include 'alfacut.f'
      include 'incldip.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     & p(mxpart,4),q(mxpart,4),omz,z,fac,ptDpg,pbDpg,ptDpb,pwsq,xr,
     & y,ymax
      integer:: j,k

      ndmax=1
      
      do j=-nf,nf
      do k=-nf,nf
        msqc(1,j,k)=0._dp
        incldip(1)=.true.
      enddo
      enddo

c--- special dipole for radiation in top decay
      call wtransform_generic(p,7,8,6,9,q,pbDpg,ptDpg,ptDpb)

      pwsq=2._dp*(q(7,4)*q(8,4)-q(7,1)*q(8,1)-q(7,2)*q(8,2)-q(7,3)*q(8,3))

c--- form of subtraction depends on whether b-quark in decay is massless or not
      if (mb < 1.e-6_dp) then
c----- massless case      
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)
        z=1._dp-omz
        xr=sqrt(pwsq/mt**2)
        ymax=(1._dp+xr)**2*z*omz/(z+xr**2*omz)
        y=2._dp*pbDpg/mt**2/(1._dp-xr)**2
        if ((z < 1._dp-aff) .and. (y > aff*ymax)) then
          incldip(1)=.false.
          return
        endif
        fac=gsq*cf*(1._dp/pbDpg*(2._dp/omz-1._dp-z)-(mt/ptDpg)**2)
      else
c----- massive case    
c-----  (no alpha-dependence at present)  
        fac=gsq*cf*((mt**2+mb**2-pwsq)/(ptDpg*pbDpg)
     &             -(mt/ptDpg)**2-(mb/pbDpg)**2)
      endif      
      
      call qqb_QQbdk(q,msq) 

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo
            
      return
      end

