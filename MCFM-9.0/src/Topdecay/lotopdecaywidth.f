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
 
      function lotopdecaywidth(mt,mb,mw,gamw)
        use singletop2_decaywidth_m
      implicit none
      include 'types.f'
      real(dp):: lotopdecaywidth
************************************************************************
*     Authors: R.K. Ellis and J. Campbell, February 2012               *
*                                                                      *
*     LO width of the top quark, including the effect of the           *
*     bottom quark mass and off shellness                              *
*                                                                      *
*     Formula is taken from Eq. (27) of                                *
*                                                                      *
*     \bibitem{Czarnecki:1990kv}                                       *
*     A.~Czarnecki,                                                    *
*     ``QCD corrections to the decay t ---> W b                        *
*       in dimensional regularization,''                               *
*     Phys.\ Lett.\ B {\bf 252}, 467 (1990).                           *
*     %%CITATION = PHLTA,B252,467;%%                                   *
*                                                                      *
*                                                                      *
************************************************************************
      
      include 'zerowidth.f'
      real(dp):: mb,mt,mt1,mw,om,omsq,be,besq,dgauss,Gamma0,
     & xlo,xhi,xi,gamw,ga,Gamma0int
      real(dp):: cachemass,cachewidth,tiny
      data cachemass,cachewidth,tiny/0._dp,0._dp,1.e-8_dp/
      common/transfer/mt1,besq,xi,ga
      save cachemass,cachewidth,tiny
!$omp threadprivate(cachemass,cachewidth,tiny,/transfer/)
      external Gamma0int


c--- check to see if result has already been computed
c     if (abs(mt*mw-mb-cachemass) < tiny) then
c       lotopdecaywidth=cachewidth
c       return
c     endif
       
      mt1=mt
      om=mw/mt
      be=mb/mt
      besq=be**2

      if (zerowidth) then
      omsq=om**2
      lotopdecaywidth=Gamma0(mt,besq,omsq)
      else
      xlo=0._dp
      xhi=(1._dp-be)**2
      ga=gamw/mw
      xi=(mt/mw)**2
      lotopdecaywidth=dgauss(Gamma0int,xlo,xhi,tiny)
      endif

c--- set-up caching variables
      cachemass=mt*mw-mb
      cachewidth=lotopdecaywidth

      return
      end



