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
 
      function higgsprop(s)
      implicit none
      include 'types.f'
      complex(dp):: higgsprop
c--- computes Higgs propagator for Higgs boson four-momentum squared s
c--- if CPscheme = .true. then it is computed in the complex pole
c---   scheme (Goria, Passarino, Rosco, arXiv:1112.5517,
c---           and c.f. Eq. (2.11) of arXiv:1206.4803)
c--- otherwise it takes the usual Breit-Wigner form
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'cpscheme.f'
      include 'first.f'
      real(dp):: s,mhbarsq,mhbar,gammahbar
      complex(dp):: cfac
      save mhbarsq,cfac

      if (CPscheme) then
c--- complex pole scheme propagator      
        if (first) then
          mhbarsq=hmass**2+hwidth**2
          mhbar=sqrt(mhbarsq)
          gammahbar=mhbar/hmass*hwidth
          cfac=cplx2(one,gammahbar/mhbar)
          first=.false.
        write(6,*)
        write(6,*)'****************************************************'
        write(6,*)'*  Using complex pole scheme for Higgs propagator  *'
        write(6,99) mhbar,gammahbar
        write(6,*)'****************************************************'
        write(6,*)
        endif
        higgsprop=cfac/(s*cfac-cplx2(mhbarsq,zip))
      else
c--- Breit Wigner propagator      
        higgsprop=one/cplx2(s-hmass**2,hmass*hwidth)
      endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end
      
      
      
