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
 
      function C0fa2m(t,qsq,msq)
      implicit none
      include 'types.f'
      complex(dp):: C0fa2m
C     C0(Pc,Pg,0,msq,msq)=
C     C0(tsq,0,qsq,0,msq,msq) (LT notation) 
C     result for qsq<0,t<0 is 
C     C0fa2m(t,qsq,msq)=(Li2(qsq/msq)-Li2(t/msq))/(t-qsq)

       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: t,qsq,msq,r,omr,ddilog
      complex(dp):: lnrat,wlog,dilogt,dilogq

      r=one-qsq/msq
      omr=qsq/msq
      wlog=lnrat(msq-qsq,msq)
      if (omr > one) then 
         dilogq=cplx1(pisqo6-ddilog(r))-wlog*cplx1(log(omr))
      else
         dilogq=cplx1(ddilog(omr))
      endif

      r=one-t/msq
      omr=t/msq
      wlog=lnrat(msq-t,msq)
      if (omr > one) then 
         dilogt=cplx1(pisqo6-ddilog(r))-wlog*cplx1(log(omr))
      else
         dilogt=cplx1(ddilog(omr))
      endif
      C0fa2m=(dilogq-dilogt)/(t-qsq)
      return
      end
