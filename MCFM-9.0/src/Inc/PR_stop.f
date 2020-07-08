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
 
!--- The variables R and P provide the Regular and Plus pieces associated
!--- with radiation from leg 1 (B1(a,b,c,is) and leg 2 (B2(a,b,c,is)
!--- In each case the parton labelling is Using the normal QM notation of putting 
!--- everything backward
!---       emitted line after emission =   a
!---       emitter before emission     =    b
!---       spectator                   =    c
!--- There is no label for he or she who is emitted.
!--- Note that in general each piece will be composed of many different
!--- dipole contributions

!--- NOTE: this is just an extension to PR_new which includes the
!---  possibility of distinguishing a heavy quark from a light quark
!---  (b=2, rather than the usual q=1)

      real(dp):: B1(-1:2,-1:2,-1:2,3),B2(-1:2,-1:2,-1:2,3)
      common/RP_stop/B1,B2
!$omp threadprivate(/RP_stop/)
