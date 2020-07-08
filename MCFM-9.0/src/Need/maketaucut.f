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
 
      subroutine maketaucut(pparton,pjet,jets,isub,passed,nd)
          use SCET
c---- J. Campbell, April 2015
c----
c---- Adapted from code written by R. Boughezal et al. for
c---- 1-jettiness calculation in W+jet events
c----
c---- Implements cut on 1-jettiness for processes with nqcdjets=1
c---- Implements cut on 0-jettiness for processes with nqcdjets=0
c----
c---  Only keeps events with tau > taucut (passed via common block)
c---  tau < taucut will be treated separately using SCET
c----
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'first.f'
      include 'kprocess.f'
      include 'hdecaymode.f'
      logical passed
      integer i,j,isub,ihard,jets
      real(dp) pparton(mxpart,4),pjet(mxpart,4),pparton_boosted(mxpart,4),
     & ptjet(3),nn1(4),nn2(4),nn3(4),tau,taua,taub,tauj,lambda,Q(4),ptest(4),
     & getdynamictau,tauc
      integer, intent(in) :: nd
      integer, save :: ipp
!$omp threadprivate(ipp)
!      logical failedndp 
!      common/failedndp/failedndp 
!!$omp threadprivate(/failedndp/) 

      logical :: bin
      common/bin/bin

c--- determine beginning of parton entries in plabel
      if (first) then
        first=.false.
        ipp=3
        do while ((ipp < mxpart ) .and. (plabel(ipp) .ne. 'pp') .and. (plabel(ipp) .ne. ''))
          ipp=ipp+1
        enddo
        if (ipp == mxpart) then
          write(6,*) 'Could not identify partons in maketaucut.f'
          stop
        endif
c        write(6,*) 'found ipp=',ipp  
      endif

      scetreweight(:) = 0._dp
      
      passed=.false.
      includeTaucutgrid(nd) = .true.

c--- safety cut on pt(H)
c      if (ntau == 1) then
c        if (pttwo(3,4,pparton) < 1._dp) return
c      endif

!---  Demand at least (ntau) hard jets after implementing the jet
!---  algorithm, otherwise the event is removed
      if (jets .lt. ntau) return

      
!--- if there are not enough partons (i.e. we're doing virtual) then return
      if ((npart .eq. ipp-3+ntau) .and. (isub .eq. 0)) then
        passed=.true.
        return
      endif

      if ((kcase == kWHbbar .or. kcase == kZHbbar .or.
     &      kcase == kWH1jet .or. kcase == kZH1jet) .and.
     &      hdecaymode =='bqba') then
          continue
      else
!--- routine is not designed for more than 3 jets
          if (jets .gt. 3) then
            write(6,*) 'Error: >3 jets found in maketaucut.f'
            write(6,*) ' jets=',jets
            stop
          endif

!--- routine is not designed for nqcdjets>2
          if (nqcdjets .gt. 2) then
            write(6,*) 'Error: unimplemented nqcdjets in maketaucut.f'
            write(6,*) ' nqcdjets=',nqcdjets
            stop
          endif
      endif

!--- reference directions for beam axes
      nn1(:)=(/0._dp,0._dp, 1._dp,1._dp/)
      nn2(:)=(/0._dp,0._dp,-1._dp,1._dp/)

!---  compute hard jet reference direction for 1-jettiness
      if (ntau > 0) then
        ptjet(:)=0._dp
        do i=0,jets-1
          ptjet(i+1)=sqrt(pjet(ipp+i,1)**2+pjet(ipp+i,2)**2)
        enddo
        ihard=ipp
        if (jets > 1) then
          if (ptjet(2) .gt. max(ptjet(1),ptjet(3))) ihard=ipp+1
        endif
        if (jets > 2) then
          if (ptjet(3) .gt. max(ptjet(1),ptjet(2))) ihard=ipp+2
        endif
        nn3(1:3)=pjet(ihard,1:3)
     &    /sqrt(pjet(ihard,1)**2+pjet(ihard,2)**2+pjet(ihard,3)**2)
        nn3(4)=1.0_dp
      endif

      if (tauboost) then 
        call getboostedvectors(pparton,pjet,isub,ipp,ihard,
     &                         nn1,nn2,nn3,pparton_boosted)
      endif

      tau=0._dp
!--- compute N-jettiness tau, by dotting the final state QCD parton
!--- 4-momentum into ni, for i=1,2 and 3 (ntau>0 only)
      do j = ipp,npart+2-isub     
        if (tauboost) then
          ptest(:)=pparton_boosted(j,:)
        else
          ptest(:)=pparton(j,:)
        endif
        taua=ptest(4)*nn1(4)-ptest(1)*nn1(1)
     &      -ptest(2)*nn1(2)-ptest(3)*nn1(3)
        taub=ptest(4)*nn2(4)-ptest(1)*nn2(1)
     &      -ptest(2)*nn2(2)-ptest(3)*nn2(3)
        if (ntau .eq. 0) then
          tau=tau+min(abs(taua),abs(taub))
        else
          tauj=ptest(4)*nn3(4)-ptest(1)*nn3(1)
     &        -ptest(2)*nn3(2)-ptest(3)*nn3(3)
          tau=tau+min(abs(taua),abs(taub),abs(tauj))
        endif
      enddo

!--- check to make sure no NaN
      if (tau .ne. tau) then
        call writeout(pparton)
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif
!      write(6,*) 'tau',tau

      if (dynamictau) then
        tauc=getdynamictau(pjet)
      else
        tauc=taucut
      endif
      
      if (bin .and. doMultitaucut) then
          ! no sampled value for taus will pass cuts, if smaller than the smallest taucut
          if (tau < smallestTaucut*(tauc/taucut)) then
              scetreweight(:) = 0._dp
              includeTaucutgrid(nd) = .false.
              return
          endif

          ! otherwise compute "weights" for other taucuts
          do j=1,size(tcutarray)
              if (tau < tcutarray(j)*(tauc/taucut)) then
                  scetreweight(j) = 0._dp
              else
                  scetreweight(j) = 1._dp
              endif
          enddo

          ! and postpone this cut for later in the *int routines
          ! and in the plotting
          if (tau < tauc) then
              includeTaucutgrid(nd) = .false.
          endif
      else
          if (tau < tauc) return
      endif

c--- at this point cut has been passed
      passed=.true.
c      if (failedndp) then 
c         write(6,*) 'WARNING, phase space point has passed, '
c         write(6,*) 'both tau cut AND has failed ndp check'
c      endif

      return
      end
      

      subroutine getboostedvectors(pparton,pjet,isub,ipp,ihard,
     & nn1,nn2,nn3,pparton_boosted)
c--- compute all four-vector quantities in Born frame where
c--- entire Born system is at rest
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'kprocess.f'
      integer j,isub,ipp,ihard
      real(dp) pparton(mxpart,4),pjet(mxpart,4),Q(4),ppin(4),ppout(4),
     & Qrest(4),Qsq,nn1(4),nn2(4),nn3(4),pparton_boosted(mxpart,4)

c--- create 4-vector Q representing Born system
      Q(:)=pparton(3,:)
      if (ipp > 4) Q(:)=Q(:)+pparton(4,:)
      if (ipp > 5) Q(:)=Q(:)+pparton(5,:)
      if (ipp > 6) Q(:)=Q(:)+pparton(6,:)
      if (ntau > 0) Q(:)=Q(:)+pjet(ihard,:)

c--- create 4-vector for Q in its own rest frame
      Qsq=Q(4)**2-Q(1)**2-Q(2)**2-Q(3)**2
      Qrest(:)=zip
      Qrest(4)=sqrt(Qsq)

c--- boost beam 1 into the same frame and compute corresponding direction, nn1
      ppin(:)=-pparton(1,:)
      call boostx(ppin,Q,Qrest,ppout)
      nn1(1:3)=ppout(1:3)/sqrt(ppout(1)**2+ppout(2)**2+ppout(3)**2)
      nn1(4)=1.0_dp

c--- boost beam 2 into the same frame and compute corresponding direction, nn2
      ppin(:)=-pparton(2,:)
      call boostx(ppin,Q,Qrest,ppout)
      nn2(1:3)=ppout(1:3)/sqrt(ppout(1)**2+ppout(2)**2+ppout(3)**2)
      nn2(4)=1.0_dp

c--- boost jet direction (if required) into the same frame and compute
c--- corresponding direction, nn3
      if (ntau > 0) then
        ppin(:)=pjet(ihard,:)
        call boostx(ppin,Q,Qrest,ppout)
        nn3(1:3)=ppout(1:3)/sqrt(ppout(1)**2+ppout(2)**2+ppout(3)**2)
        nn3(4)=1.0_dp
      endif

c--- boost partons into the same frame, filling pparton_boosted
      do j = ipp,npart+2-isub     
        ppin(:)=pparton(j,:)
        call boostx(ppin,Q,Qrest,ppout)
        pparton_boosted(j,:)=ppout(:)
      enddo

      return
      end
