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
 
      subroutine nplotter_ZZlept(p,wt,wt2,switch)
      implicit none
      include 'types.f'
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer:: equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'interference.f'
      real(dp):: p(mxpart,4),wt,wt2,m3456,pt34,pttwo
      integer:: switch,n,nplotmax
      integer tag
      include 'first.f'
      common/nplotmax/nplotmax

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************
      
      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        goto 99
      else
c--- Add event in histograms
        tag=tagplot
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

      m3456=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &     -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &     -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &     -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2
      m3456=sqrt(max(m3456,zip))
      
      pt34=pttwo(3,4,p)
       
************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- "n" will count the number of histograms
      n=nextnplot

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale

c--- Plots of m(3456) in specific regions
      call bookplot(n,tag,'10 < m(3456) < 2010',
     & m3456,wt,wt2,10._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < m(3456) < 2010',
     & m3456,wt,wt2,130._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'300 < m(3456) < 2020',
     & m3456,wt,wt2,300._dp,2020._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < m(3456) < 130',
     & m3456,wt,wt2,10._dp,130._dp,5._dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'pt(Z)',
     & pt34,wt,wt2,0._dp,2._dp,0.02_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'+INTEGRAL+ pt(Z)',
     & pt34,wt,wt2,0._dp,10._dp,0.1_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'50 < m(3456) < 250',
     & m3456,wt,wt2,50._dp,250._dp,2._dp,'log')
      n=n+1
      
c--- usual plots for 3+4
      call autoplot2(p,34,3,4,tag,wt,wt2,n)

c--- usual plots for 5+6
      call autoplot2(p,56,5,6,tag,wt,wt2,n)

      if (interference) then
c--- usual plots for 3+6
        call autoplot2(p,36,3,6,tag,wt,wt2,n)

c--- usual plots for 4+5
        call autoplot2(p,45,4,5,tag,wt,wt2,n)
      endif

c--- usual plots for 3+4+5+6
      call autoplot4(p,3456,3,4,5,6,tag,wt,wt2,n)

c--- additional plots that may be present at NLO       
      if (abs(p(7,4)) > 1.e-8_dp) then
        call autoplot1(p,7,tag,wt,wt2,n)
      else
        n=n+2
      endif

************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1


c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
      return
      end
