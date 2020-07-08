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
 
      include 'nplot.f'
      integer maxnbin
      parameter(maxnbin=200)
      real(dp)::HIST(nplot,maxnbin),XHIS(nplot,maxnbin),
     & HDEL(nplot),
     & HMIN(nplot),HMAX(nplot),HAVG(nplot),HINT(nplot),HSIG(nplot)
      COMMON/HISTOR/HIST,XHIS,HDEL,HMIN,HMAX,HAVG,HINT,HSIG

      CHARACTER TITLE*100,BOOK*3
      COMMON/HISTOC/BOOK(nplot),TITLE(nplot)
      INTEGER NBIN(nplot),NHIST
      integer(kind=8) IHIS(nplot,maxnbin),IUSCORE(nplot),
     & IOSCORE(nplot),IENT(nplot)
      COMMON/HISTOI/IHIS,IUSCORE,IOSCORE,IENT,NBIN,NHIST

      integer, parameter:: tagbook=1, tagplot=2

