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
 
      function C0fb2m(u,msq)
      implicit none
      include 'types.f'
      complex(dp):: C0fb2m
C     C0(Pc,Pg,0,msq,msq)=
C     C0(msq,0,u,0,msq,msq) (LT notation) 
C  C0fb2m(u,msq)=(Pi^2/6-Li2[u/m2])/(u-msq) with u=(Pc+Pg)^2;
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: u,msq,ubar,r,omr,ddilog
      complex(dp):: lnrat,wlogu,dilogu
      ubar=u-msq
      r=-ubar/msq
      omr=one-r
      if (omr > one) then 
         wlogu=lnrat(-ubar,msq)
         dilogu=cplx1(pisqo6-ddilog(r))-wlogu*cplx1(log(omr))
      else
         dilogu=cplx1(ddilog(omr))
      endif
      C0fb2m=(cplx1(pisqo6)-dilogu)/ubar
      return
      end
