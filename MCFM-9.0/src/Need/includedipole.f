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
 
      !----------------------------------------------------------------------
      !      This replaces the original includedipole
      !      It calls the original and then the user one
      logical function includedipole(nd,ptrans)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      real(dp) ptrans(mxpart,4)
      integer nd
      logical :: mcfm_includedipole
      
      ! it looks like all (nd .ne. 0) automatically have their momenta
      ! stored in the ptilde common; do this also for nd=0
      if (nd .eq. 0) call storeptilde(nd,ptrans)

      ! first call the original MCFM includedipole
      includedipole = mcfm_includedipole(nd,ptrans)

      end 

      logical function mcfm_includedipole(nd,ptrans) 
     &                 result(mcfmincdipole)
       use types
       use m_gencuts
c--- This function returns TRUE if the specified point ptrans,
c--- corresponding to dipole nd (nd=0 => real radiation),
c--- should be included 
      implicit none
      include 'mxpart.f'
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      include 'kprocess.f'
      include 'frag.f'
      include 'phot_dip.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      include 'notag.f'
      include 'taucut.f'
      include 'hdecaymode.f'
      include 'ewcorr.f'
      include 'energy.f'
      include 'runstring.f'
      include 'first.f'

      real(dp) ptrans(mxpart,4),pjet(mxpart,4),rcut,pt,pttwo,dot
      integer j,nd,isub
      logical failedgencuts,photoncuts,makecuts,filterWbbmas,
     &     photonfailed,filterW_bjet,is_photon
      integer count_photo,nphotons
      logical passed_frix,iso, passed_taucut

      common/rcut/rcut
      common/makecuts/makecuts

c--- default: include this contribution
      mcfmincdipole=.true.
      
c--- isub=1 for dipole subtractions, isub=0 for real radiation
      if (nd .gt. 0) then
        isub=1
      else
        isub=0
      endif

      nphotons=count_photo() 
      if (nphotons .gt. 0) then 
c--- Photons: Frixione isolation cuts if no fragmentation included 
         if (frag .eqv. .false.) then
            do j=3,mxpart 
               if(is_photon(j)) then 
                  call frix(ptrans,passed_frix,j,isub)
                  if(passed_frix.eqv..false.) then 
                     mcfmincdipole=.false.
                     return 
                  endif
               endif
            enddo
            call genclustphotons(ptrans,rcut,pjet,isub)
         else 
c--- Photons: not Frixione, need fragmentation and isolation
!---- do not want to cluster partons inside of jet cone, allow 
!---- isolation to describe these regions, therefore use 
!---- genclustphotons here 
            call genclustphotons(ptrans,rcut,pjet,isub)
!            call genclust2(ptrans,rcut,pjet,isub)
c---  Isolate photon
            do j=3,mxpart
            if (is_photon(j)) then 
               if (iso(ptrans,j,isub,nd) .eqv. .false.)then
                  mcfmincdipole=.false.
                  return 
               endif
            endif
            enddo
         endif
c--- check the photon cuts 
         photonfailed=photoncuts(pjet)
         if (photonfailed) then
            mcfmincdipole=.false.
            return
         endif 
      else
c--- No photons: the usual case
         call genclust2(ptrans,rcut,pjet,isub)
      endif


c--- perform mass cuts
      call masscuts(pjet,*999)
          
c--- fill ptilde array as persistent storage for the jet momenta
      ! GPS written in compact f90 form - in case we wish to copy it elsewhere
      ! (e.g. earlier, so that jets are defined even when includedipole is false)
      ptildejet(nd,1:npart+2,1:4)=pjet(1:npart+2,1:4)

c--- for the Wbb process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing 
c--- NOTE: only for process numbers > 400 (20 and 25 should be handled normally) 
      if ((kcase==kWbbmas) .and. (nproc .gt. 400)) then
        mcfmincdipole=filterWbbmas()
      if (mcfmincdipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) mcfmincdipole=.false.
        endif
      goto 99
      endif
      

c--- for the Wb+X process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing   
      if (kcase==kW_bjet) then
        mcfmincdipole=filterW_bjet()
      if (mcfmincdipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) mcfmincdipole=.false.
        endif
      goto 99
      endif

     
      if (usescet) then
c---  for SCET calculation do not check jets, make tau cut instead
        call maketaucut(ptrans,pjet,jets,isub,passed_taucut,nd)
        mcfmincdipole=passed_taucut
c           write(6,*) 'includedipole: nd,mcfmincdipole',nd,mcfmincdipole
c           if (passed_taucut .eqv. .false.) write(6,*) 'tau failed: ',nd
        if (mcfmincdipole .eqv. .false.) return
      else
c--- for a normal calculation,      
c--- if the number of jets is not correct, then do not include dipole
        if ((clustering .and. (jets .ne. nqcdjets-notag)
     &         .and. (inclusive .eqv. .false.)) .or.
     &      (clustering .and. (jets .lt. nqcdjets-notag)
     &         .and. (inclusive .eqv. .true.))) then
            mcfmincdipole=.false.
            return
        endif
      endif

c--- check the lepton cuts, if necessary
      if (makecuts) then
        failedgencuts=gencuts(pjet,jets)
        if (failedgencuts) then
          mcfmincdipole=.false.
          return
        endif
      endif
      
   99 continue
      
      return

  999 continue
      mcfmincdipole=.false.
      return    
      
      end

