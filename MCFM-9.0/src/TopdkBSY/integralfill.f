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
 
      subroutine integralfill(p)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
C-----Authors: John Campbell and Keith Ellis, November 2011
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'scale.f'
      include 'massiveintegrals.f'
      real(dp):: p(mxpart,4)
      integer:: nu,p2,p3,j
      real(dp):: mt2,xbeta2,s12,s13,s23,twop1Dp2,twop1Dp3,gram

      do j=1,2
      if (j==1) then 
      p2=2
      p3=3
      elseif (j==2) then
      p2=3
      p3=2
      endif

      s12=(p(1,4)+p(p2,4))**2
      s13=(p(1,4)+p(p3,4))**2
      s23=(p(p2,4)+p(p3,4))**2
      do nu=1,3
      s12=s12-(p(1,nu)+p(p2,nu))**2
      s13=s13-(p(1,nu)+p(p3,nu))**2
      s23=s23-(p(p2,nu)+p(p3,nu))**2
      enddo 
      mt2=mt**2
      twop1Dp2=s12-mt2
      twop1Dp3=s13-mt2
      gram=twop1Dp2*twop1Dp3-mt2*s23
      xbeta2=1d0-4d0*mt2/s23

      if (j == 1) then
      I41x2x4x3=qlI4(mt2,zip,mt2,zip,s12,s13,zip,mt2,mt2,zip,musq,0)
      I31x23x4=qlI3(s23,mt2,mt2,zip,zip,mt2,musq,0)
      I3m1x23x4=qlI3(s23,mt2,mt2,mt2,mt2,zip,musq,0)

      I32x3x41=qlI3(s23,zip,zip,zip,zip,zip,musq,0)
      I3m2x3x41=qlI3(s23,zip,zip,mt2,mt2,mt2,musq,0)

      I2m=qlI2(mt2,zip,mt2,musq,0)
      Im2m=qlI2(zip,mt2,mt2,musq,0)
      I23=qlI2(s23,zip,zip,musq,0)
      F2m23=qlI2(s23,mt2,mt2,musq,0)-Im2m
      I2h23=qlI2(s23,zip,zip,musq,0)-I2m+ctwo
      endif
      
      I41x2x3x4(j)=qlI4(mt2,zip,zip,mt2,s12,s23,mt2,zip,zip,zip,musq,0)
      I4m1x2x3x4(j)=qlI4(mt2,zip,zip,mt2,s12,s23,zip,mt2,mt2,mt2,musq,0)

      I312x3x4(j)=qlI3(s12,zip,mt2,mt2,zip,zip,musq,0)
      I3m12x3x4(j)=qlI3(s12,zip,mt2,zip,mt2,mt2,musq,0)

      I313x2x4(j)=qlI3(s13,zip,mt2,mt2,zip,zip,musq,0)
      I3m13x2x4(j)=qlI3(s13,zip,mt2,zip,mt2,mt2,musq,0)

      F41x2x3x4(j)=I41x2x3x4(j)-I31x23x4/twop1Dp2
      F4m1x2x3x4(j)=I4m1x2x3x4(j)-I3m1x23x4/twop1Dp2

      F212(j)=qlI2(s12,zip,mt2,musq,0)-I2m

      I461x2x3x4(j)=0.5d0/gram*(
     &   -twop1Dp2**2*s23*I41x2x3x4(j)
     &   +(twop1Dp2+2d0*mt2)*s23*I31x23x4
     &   +2d0*twop1Dp2**2*I312x3x4(j)
     &   +twop1Dp2*s23*I32x3x41)

      I46m1x2x3x4(j)=0.5d0/gram*(
     %  -twop1Dp2**2*s23*xbeta2*I4m1x2x3x4(j)
     &   +twop1Dp2*s23*xbeta2*I3m1x23x4
     &   +2d0*twop1Dp2*(twop1Dp2+2d0*mt2)*I3m12x3x4(j)
     &   +(twop1Dp2+2d0*mt2)*s23*I3m2x3x41)
      enddo
      return
 
      end
