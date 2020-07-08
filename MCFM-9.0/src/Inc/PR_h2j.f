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
 
c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (H1(a,b,c,is) and leg 2 (H2(a,b,c,is)
c--- In each case the parton labelling uses the normal QM notation of
c--- putting everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

c--- NOTE: this is just an extension to PR_new which enables the
c---  division of matrix elements into colour structures in the
c---  same way as 'msq_struc.f'

c      integer:: igg_ab,igg_ba,igg_sym,iqq_a,iqq_b,iqq_i,
c     & igggg_a,igggg_b,igggg_c,iqr
c      parameter(igg_ab=4,igg_ba=5,igg_sym=6)
c--- Note that the 4-quark and 4-gluon pieces are never needed simultaneously,
c--- so we can reuse the same indices to save memory
c      parameter(iqq_a=1,iqq_b=2,iqq_i=3)
c      parameter(igggg_a=1,igggg_b=2,igggg_c=3)
c--- One extra parameter needed for non-identical quark pieces
c      parameter(iqr=7)

      real(dp):: 
     & H1(-1:1,-1:1,-1:1,8,3),H2(-1:1,-1:1,-1:1,8,3)
      common/RP_h2j/H1,H2
!$omp threadprivate(/RP_h2j/)
