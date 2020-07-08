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
 
      function ttbqqbg_sq(i1,i2,i9,i5,i3,i6,i8)
      implicit none
      include 'types.f'
      real(dp):: ttbqqbg_sq
c--- returns the summed squared helicity amplitudes for the
c--- ttbqqbg amplitudes
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i1,i2,i9,i5,i3,i6,i8
      complex(dp):: sq,tq,qq,rq,
     & ttbqqbsqpp,ttbqqbsqpm,ttbqqbsqmp,ttbqqbsqmm,
     & ttbqqbtqpp,ttbqqbtqpm,ttbqqbtqmp,ttbqqbtqmm,
     & ttbqqbqqpp,ttbqqbqqpm,ttbqqbqqmp,ttbqqbqqmm,
     & ttbqqbrqpp,ttbqqbrqpm,ttbqqbrqmp,ttbqqbrqmm
      real(dp):: appsq,apmsq,ampsq,ammsq
      
c-- for q-qb , there are four colour amplitudes:
c---    sq        proportional to Ta(it,i1)*delta(i2,ib)
c---    tq        proportional to Ta(i2,ib)*delta(it,i1)
c---    qq        proportional to Ta(i2,i1)*delta(it,ib)
c---    rq        proportional to Ta(it,ib)*delta(i2,i1)
      sq=ttbqqbsqpp(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqpp(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqpp(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqpp(i1,i2,i9,i5,i3,i6,i8)
      appsq=xn*cf*(
     & +xn*abs(sq)**2+xn*abs(tq)**2
     & +xn*abs(qq)**2+xn*abs(rq)**2
     & +2._dp*Dble((sq+tq)*conjg(rq+qq)))
      
      sq=ttbqqbsqpm(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqpm(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqpm(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqpm(i1,i2,i9,i5,i3,i6,i8)
      apmsq=xn*cf*(
     & +xn*abs(sq)**2+xn*abs(tq)**2
     & +xn*abs(qq)**2+xn*abs(rq)**2
     & +2._dp*Dble((sq+tq)*conjg(rq+qq)))
      
      sq=ttbqqbsqmp(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqmp(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqmp(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqmp(i1,i2,i9,i5,i3,i6,i8)
      ampsq=xn*cf*(
     & +xn*abs(sq)**2+xn*abs(tq)**2
     & +xn*abs(qq)**2+xn*abs(rq)**2
     & +2._dp*Dble((sq+tq)*conjg(rq+qq)))
      
      sq=ttbqqbsqmm(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqmm(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqmm(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqmm(i1,i2,i9,i5,i3,i6,i8)
      ammsq=xn*cf*(
     & +xn*abs(sq)**2+xn*abs(tq)**2
     & +xn*abs(qq)**2+xn*abs(rq)**2
     & +2._dp*Dble((sq+tq)*conjg(rq+qq)))
      
      ttbqqbg_sq=appsq+apmsq+ampsq+ammsq

      return
      end
      
