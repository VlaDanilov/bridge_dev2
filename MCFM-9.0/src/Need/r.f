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
 
      function rpure(pi,pj)
          use types
          implicit none
          include 'userap.f'
          
          real(dp) :: rpure
          real(dp), intent(in) :: pi(4),pj(4)

          real(dp) :: pti2, ptj2, ei, ej, r1, r2, biti, bitj, delphi, dely
          real(dp), parameter :: tiny = 1.e-6_dp

          pti2=pi(1)**2+pi(2)**2
          ptj2=pj(1)**2+pj(2)**2

          if (userap) then
!--- use rapidities (not pseudorapidities)
            ei=pi(4)
            ej=pj(4)
          else
            ei=sqrt(pti2+pi(3)**2)
            ej=sqrt(ptj2+pj(3)**2)
          endif

          biti=pi(3)/ei
          bitj=pj(3)/ej

          if  ((abs(1d0+biti) < tiny) .or. (abs(1d0-biti) < tiny)
     &     .or.(abs(1d0+bitj) < tiny) .or. (abs(1d0-bitj) < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
              dely=100d0
          else
              r1=(1._dp+biti)*(1._dp-bitj)/((1._dp+bitj)*(1._dp-biti))
              dely=0.5_dp*log(r1)
          endif

          r2= (pi(1)*pj(1)+pi(2)*pj(2))/sqrt(pti2*ptj2)
          if (r2 > +0.9999999_dp) r2=+1._dp
          if (r2 < -0.9999999_dp) r2=-1._dp
          delphi=acos(r2)

          rpure = sqrt(dely**2+delphi**2)
      end function
      
      function r(p,i,j)
      implicit none
      include 'types.f'
      real(dp):: r
c----calculate the jets separation between p(i) and p(j)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'userap.f'
      real(dp):: p(mxpart,4),r1,r2,dely,delphi,ei,ej,pti2,ptj2,
     & biti,bitj
      integer:: i,j
      real(dp), parameter:: tiny=1.e-6

      pti2=p(i,1)**2+p(i,2)**2
      ptj2=p(j,1)**2+p(j,2)**2

      if (userap) then
!--- use rapidities (not pseudorapidities)
        ei=p(i,4)
        ej=p(j,4)
      else
        ei=sqrt(pti2+p(i,3)**2)
        ej=sqrt(ptj2+p(j,3)**2)
      endif

! sanity check in case this is called when a momentum is not filled
      if ((ei < 1.e-13_dp) .or. (ej < 1.e-13_dp)) then
        r=0._dp
        return
      endif

c      r1= (ei+p(i,3))*(ej-p(j,3))/
c     &   ((ej+p(j,3))*(ei-p(i,3)))
      biti=p(i,3)/ei
      bitj=p(j,3)/ej
      if  ((abs(1d0+biti) < tiny) .or. (abs(1d0-biti) < tiny)
     & .or.(abs(1d0+bitj) < tiny) .or. (abs(1d0-bitj) < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
        dely=100d0
      else
        r1=(one+biti)*(one-bitj)/((one+bitj)*(one-biti))
        dely=0.5_dp*dlog(r1)
      endif

      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))/sqrt(pti2*ptj2)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      r=sqrt(dely**2+delphi**2)
      
      return
      end
      
      function delphi(pi,pj)
        include 'types.f'

        real(dp) :: delphi
        real(dp), intent(in) :: pi(4), pj(4)

        real(dp) :: pti2, ptj2

        pti2 = pi(1)**2+pi(2)**2
        ptj2 = pj(1)**2+pj(2)**2
        delphi = acos((pi(1)*pj(1)+pi(2)*pj(2))/sqrt(pti2*ptj2))

      end function

