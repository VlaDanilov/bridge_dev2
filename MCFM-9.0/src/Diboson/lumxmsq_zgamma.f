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
 
      subroutine lumxmsq_zgamma(p,xx,z1,z2,QB,order,xmsq,central)
          use SCET
          implicit none
          include 'types.f'
          include 'nf.f'
          include 'constants.f'
          include 'mxpart.f'
          include 'facscale.f'
          include 'qcdcouple.f'
          include 'tiny.f'
          include 'taucut.f'

          real(dp), intent(in) :: p(mxpart,4), xx(2), z1, z2, QB(2)
          integer, intent(in) :: order
          real(dp), intent(out) :: xmsq
          logical, intent(in) :: central

          real(dp) :: soft1(-1:1), soft2(-1:3)
          real(dp) :: beama0(-5:5), beamb0(-5:5),
     &                beama1(-5:5, -1:1), beamb1(-5:5, -1:1),
     &                beama2(-5:5, -1:3), beamb2(-5:5, -1:3)
          real(dp) :: hard(2), bit, assemble
          real(dp) :: assemble_zgamma

          real(dp) :: msq(-nf:nf, -nf:nf, 0:2), ggcontrib
          real(dp) :: Q, dot, msqpow(-5:5,-5:5)

          integer :: ih1, ih2
          common/density/ih1,ih2

          real(dp) :: origtaucut

          integer j,k,m

          call softqqbis(order,soft1,soft2)

          if (order >= 0) then
              call fdist(ih1,xx(1),facscale,beama0)
              call fdist(ih2,xx(2),facscale,beamb0)
          endif
          if (order >= 1) then
              call xbeam1bis(ih1,z1,xx(1),QB(1),beama1)
              call xbeam1bis(ih1,z2,xx(2),QB(2),beamb1)
          endif
          if (order >= 2) then
              call xbeam2bis(ih1,z1,xx(1),QB(1),beama2)
              call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2)
          endif

          xmsq = 0

          call zgam_mat(p,msq)

          xmsq = assemble_zgamma(p,xx,order,soft1,soft2,
     &              beama0,beamb0,beama1,beamb1,beama2,beamb2,msq)

          if ((order >= 2) .and. (onlypowcorr .eqv. .false.)) then
            call gg_zgam(p,ggcontrib)
            xmsq = xmsq + ggcontrib*beama0(0)*beamb0(0)
          endif

          ! only fill scetreweight for the central value binning
          ! not for pdf uncertainties and scale variation
          if (central .and. doMultitaucut) then
              scetreweight(:) = 0._dp
              ! is there any reasonable way for the other taucut values to give a non-zero result
              ! when xmsq is zero? I don't think so..
              if (xmsq /= 0._dp) then
                  origtaucut = taucut
                  do m=1,size(tcutarray)
                    taucut = tcutarray(m)
                    scetreweight(m) = assemble_zgamma(p,xx,order,soft1,soft2,
     &                      beama0,beamb0,beama1,beamb1,beama2,beamb2,msq)
                    if ((order >= 2) .and. (onlypowcorr .eqv. .false.)) then
                        scetreweight(m) = scetreweight(m) + ggcontrib*beama0(0)*beamb0(0)
                    endif
                  enddo
                  taucut = origtaucut
                  scetreweight(:) = scetreweight(:) / xmsq
              endif
          endif

      end subroutine

      function assemble_zgamma(p,xx,order,soft1,soft2,beama0,beamb0,beama1,beamb1,beama2,beamb2,msq)
          implicit none
          include 'types.f'
          include 'nf.f'
          include 'constants.f'
          include 'mxpart.f'
          include 'facscale.f'
          include 'qcdcouple.f'
          include 'tiny.f'
          include 'taucut.f'

          real(dp) :: assemble_zgamma
          real(dp), intent(in) :: p(mxpart,4), xx(2)
          integer, intent(in) :: order
          real(dp), intent(in) :: soft1(-1:1), soft2(-1:3)
          real(dp), intent(in) :: beama0(-5:5), beamb0(-5:5),
     &                beama1(-5:5, -1:1), beamb1(-5:5, -1:1),
     &                beama2(-5:5, -1:3), beamb2(-5:5, -1:3)
          real(dp), intent(in) :: msq(-nf:nf, -nf:nf, 0:2)

          real(dp) :: assemble, dot
          real(dp) :: Q, msqpow(-5:5,-5:5)
          real(dp) :: bit, hard(2), tauc, getdynamictau

          integer :: j,k

          if (dynamictau) then
            tauc=getdynamictau(p)
          else
            tauc=taucut
          endif

! compute power corrections if required
          if ((incpowcorr) .or. (onlypowcorr)) then
              Q=sqrt(two*dot(p,1,2))
              call powcorr_qa(order,tauc,xx(1),xx(2),Q,beama0,beamb0,msqpow)
          endif

          assemble_zgamma = 0._dp

          do j=-nf,nf
              k=-j
              if (j*k > 0 .or. j*k == 0) cycle
              if (msq(j,k,0) < tiny) cycle

              hard(1) = msq(j,k,1)/msq(j,k,0)/ason2pi
              hard(2) = msq(j,k,2)/msq(j,k,0)/ason2pi**2

              bit = assemble(order,tauc,beama0(j),beamb0(k),
     &                         beama1(j,:),beamb1(k,:),
     &                         beama2(j,:),beamb2(k,:),soft1,soft2,hard)
     
              if (incpowcorr) then
                bit=bit+msqpow(j,k)+msqpow(j,0)+msqpow(0,k)
              endif
              if (onlypowcorr) then
                bit=msqpow(j,k)+msqpow(j,0)+msqpow(0,k)
              endif

              bit = bit * msq(j,k,0)

              assemble_zgamma = assemble_zgamma + bit
        
          enddo 
      end function 
