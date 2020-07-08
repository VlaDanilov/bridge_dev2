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
 
      subroutine mcfm_exit(xinteg,xinteg_err)
          use omp_lib
          use types
          use MCFMStorage
          use MCFMPrint
          use SCET, only : doMultitaucut
          use PDFerrors, only : doPDFerrors
          use Scalevar, only : doScalevar
      implicit none
      include 'mpicommon.f'

      real(dp), intent(in) :: xinteg, xinteg_err
     
      if (rank == 0) then
          call finalizeStorage

#ifdef HAVE_LHAPDF
          call printallcross()

          if (doPDFerrors) then
            call printPDFuncertainties()
          endif
#else
          call printcross(xinteg,xinteg_err, chisqMax())
#endif

#ifdef HAVE_LHAPDF
          if (doPDFerrors) then
            call printPDFuncertainties()
          endif
#endif

          if (doScalevar) then
              call printScaleuncertainties()
          endif

          if (doMultitaucut) then
              call printTaucuts()
          endif

          call writeAllHistograms()
          call writereference()
      endif

      end
