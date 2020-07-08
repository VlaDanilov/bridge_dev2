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
c--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
c--- In each case the parton labelling is Using the normal QM notation of putting 
c--- everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

c--- Additional label (1-8) in this array, to represent the partons
c--- in the final state, according to the look-up parameters
c--- that are also defined here

c--- code is as follows:
c---  g=0, q=1, a=-1, r=2, b=-2 --- "f" for final
      integer:: gf_gf,qf_af,qf_qf,qf_rf,bf_rf,rf_bf,qf_gf,af_gf
      parameter (gf_gf=1,qf_af=2,qf_qf=3,qf_rf=4,
     & bf_rf=5,rf_bf=6,qf_gf=7,af_gf=8)
      real(dp):: 
     & S1(-1:1,-1:1,-1:1,8,0:2,3),S2(-1:1,-1:1,-1:1,8,0:2,3)
      common/SP_twojet/S1,S2
!$omp threadprivate(/SP_twojet/)
