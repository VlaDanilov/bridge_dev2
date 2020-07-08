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
 
      subroutine nplotter_trigam(p,wt,wt2,switch)
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
      real(dp):: p(mxpart,4),wt,wt2,pt,pord(mxpart,4),pt3,pt4,pt5
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

c--- order photons by pt
      pt3=pt(3,p)
      pt4=pt(4,p)
      pt5=pt(5,p)
      
      pord(:,:)=p(:,:)
      
      if ((pt3 > pt4) .and. (pt3 > pt5)) then
        pord(3,:)=p(3,:)
        if (pt4 > pt5) then
        pord(4,:)=p(4,:)
        pord(5,:)=p(5,:)
        else
        pord(4,:)=p(5,:)
        pord(5,:)=p(4,:)
        endif
      endif  
      if ((pt4 > pt3) .and. (pt4 > pt5)) then
        pord(3,:)=p(4,:)
        if (pt3 > pt5) then
        pord(4,:)=p(3,:)
        pord(5,:)=p(5,:)
        else
        pord(4,:)=p(5,:)
        pord(5,:)=p(3,:)
        endif
      endif  
      if ((pt5 > pt3) .and. (pt5 > pt4)) then
        pord(3,:)=p(5,:)
        if (pt3 > pt4) then
        pord(4,:)=p(3,:)
        pord(5,:)=p(4,:)
        else
        pord(4,:)=p(4,:)
        pord(5,:)=p(3,:)
        endif
      endif  
             
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

c--- usual plots for highest pt photon
      call autoplot1(pord,3,tag,wt,wt2,n)

c--- usual plots for 2n.e-_dphighest pt photon
      call autoplot1(pord,4,tag,wt,wt2,n)

c--- usual plots for 3r.e-_dphighest pt photon
      call autoplot1(pord,5,tag,wt,wt2,n)

c--- usual plots for 3+4+5
      call autoplot3(p,345,3,4,5,tag,wt,wt2,n)

c--- additional plots that may be present at NLO       
      if (abs(p(6,4)) > 1.e-8_dp) then
        call autoplot1(p,6,tag,wt,wt2,n)
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
