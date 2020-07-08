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
 
      function lowint(r,wgt)
        use ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
        use MCFMStorage
        use Scalevar
        use PDFerrors
#ifdef HAVE_LHAPDF
        use LHAPDF, only: getalphas
#endif
        use SCET
        use singletop2_m
        use singletop2_scale_m
        use VVconfig_m
        use bbfrac_m, only : bbfrac
        use m_gencuts, only : reweight_user, enable_reweight_user
      implicit none
      real(dp):: lowint
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'noglue.f'
      include 'kprocess.f'
      include 'maxwt.f'
      include 'phasemin.f'
      include 'wts_bypart.f'
      include 'stopscales.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'outputoptions.f'
      include 'runstring.f'
      include 'energy.f'
      include 'VVstrong.f'
      include 'dm_params.f'
      include 'initialscales.f'
      include 'toploopgaga.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'badpoint.f'
      include 'ewcorr.f'
      include 'scalevar.f'
      include 'nflav.f'
      include 'x1x2.f'
      include 'bypart.f'
      include 'taucut.f' ! for usescet
c--- APPLgrid - to use grids
      include 'ptilde.f'
      include 'APPLinclude.f'
      double precision psCR
c--- APPLgrid - end
      integer:: pflav,pbarflav
c--- To use VEGAS random number sequence :
      real(dp):: ran2
      integer:: ih1,ih2,j,k,nvec,sgnj,sgnk,ii,i1,i2,i3,i4
      integer:: i,t
      integer:: itrial
      real(dp):: alphas,msqtrial,xmsqvar(2),
     & fx1up(-nf:nf),fx2up(-nf:nf),fx1dn(-nf:nf),fx2dn(-nf:nf)
      real(dp):: r(mxdim),W,xmsq,val,val2,ptmp,
     & fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     & pswt,rscalestart,fscalestart,
     & fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     & fxb1(-nf:nf),fxb2(-nf:nf)
      real(dp):: wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
      real(dp):: xmsqjk_noew,msq_noew(-nf:nf,-nf:nf),xmsq_noew,lowint_noew
      real(dp):: msq1(-nf:nf,-nf:nf),
     & msq4(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: flux,vol,vol_mass,vol3_mass,vol_wt,BrnRat
      logical:: bin,includedipole,checkpiDpjk
      external qg_tbq,BSYqqb_QQbdk_gvec,qqb_QQbdk,qg_tbqdk,qg_tbqdk_gvec,
     & qqb_Waa,qqb_Waa_mad
     & qqb_Zbbmas,qqb_Zbbmas,qqb_totttZ,qqb_totttZ_mad
      common/density/ih1,ih2
      common/bin/bin
      common/BrnRat/BrnRat
      include 'bqscale.f'
      external qq_tchan_ztq,qq_tchan_ztq_mad
      external qq_tchan_htq,qq_tchan_htq_mad,qq_tchan_htq_amp
      external qqb_gamgam_g,qqb_gmgmjt_gvec,gg_hzgamg,gg_hg_zgam_gvec

      lowint=0._dp
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0._dp
      
      W=sqrts**2
      p(:,:)=0._dp
      pjet(:,:)=0._dp

      currentNd = 0
      currentPDF=0

      if (maxPDFsets > 0) pdfreweight(:) = 0._dp

      call gen_lops(r,p,pswt,*999)

      if (all(.not. ieee_is_nan(p(1:npart,:))) .eqv. .false.) then
          if (debug) then
              write(6,*) 'Discarding NaN or infinite phase space point'
          endif
          goto 999
      endif

      nvec=npart+2
      call dotem(nvec,p,s)

c--- (moved to includedipole) impose cuts on final state
c      if (kcase.ne.kvlchk6 .and. kcase.ne.ktautau) then
c        call masscuts(p,*999)
c      endif

! small safety cuts
      call smallnew(p,npart,*999)

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      includeTaucutgrid(0) = .true.
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif

  771 continue    
      if (dynamicscale) then
          call scaleset(initscale,initfacscale,p)
      else
          call usescales(initscale,initfacscale)
      endif

      if (doPDFAlphas) then
          call updateAlphas(scale)
      endif
      
      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts
      if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp)
     &.or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
         goto 999
      endif

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)

      if ((doScalevar .and. currentPDF == 0 .and. bin) .and. (foundpow .eqv. .false.)) then
        itrial=1
      endif
   66 continue

c--- Calculate the required matrix elements      
      if     (kcase==kW_only) then
        call qqb_w(p,msq)
      elseif (kcase==kW_1jet) then
        call qqb_w_g(p,msq)
c        call qqb_w_gbis(p,msq1)
c        do j=-nf,nf
c        do k=-nf,nf
c        if (msq(j,k) .ne. 0._dp) write(6,*) msq(j,k)/msq1(j,k)
c        enddo
c        enddo
c        stop
      elseif (kcase==kWgamma) then
        call qqb_wgam(p,msq)
      elseif (kcase==kWgajet) then
         call qqb_wgam_g(p,msq)
      elseif (kcase==kWga2jt) then
         call qqb_waj_g(p,msq)
      elseif (kcase==kWbfrmc) then
        call qqb_wbfromc(p,msq)
      elseif (kcase==kW_cjet) then
        call qqb_w_cjet(p,msq)
      elseif (kcase==kWcjet0) then
        call qqb_w_cjet_massless(p,msq)
      elseif (kcase==kWbbmas) then
        call qqb_wbbm(p,msq)
      elseif (kcase==kWbbjem) then
        call qqb_wbbm_g(p,msq)
      elseif (kcase==kWttmas) then
        call qqb_wbbm(p,msq)
      elseif (kcase==kWbbbar) then
        call qqb_wbb(p,msq)
      elseif (kcase==kW_2jet) then
        call qqb_w2jet(p,msq)
      elseif (kcase==kW_3jet) then
        call qqb_w2jet_g(p,msq)
      elseif (kcase==kWbbjet) then
        call qqb_wbb_g(p,msq)
      elseif (kcase==kZ_only) then
        call qqb_z(p,msq)
        if (kewcorr /= knone) then
          msq_noew=msq
          if     (kewcorr == ksudakov) then
            call qqb_z_ew_sudakov(p,msq)
          elseif (kewcorr == kexact) then
            call qqb_z_ew_exact(p,msq)
          endif
        endif
      elseif (kcase==kgg2lep) then
        call ggdilep(p,msq)
      elseif (kcase==kZ_1jet) then
        call qqb_z1jet(p,msq)
c        call qqb_z1jetbis(p,msq1)
c        do j=-nf,nf
c        do k=-nf,nf
c        if (msq(j,k) .ne. 0._dp) write(6,*) msq(j,k)/msq1(j,k)
c        enddo
c        enddo
c        pause
      elseif (kcase==kZ_2jet) then
        call qqb_z2jet(p,msq)
      elseif (kcase==kZ_3jet) then
        call qqb_z2jet_g(p,msq)
      elseif (kcase==kZgamma) then
        call set_anomcoup(p)
        ! distinction between decays is inside the subroutine
        call qqb_zgam_new(p,msq)
      elseif (kcase==kZ_2gam) then
        call qqb_zaa(p,msq)
      elseif (kcase==kW_2gam) then
c        if (checkpiDpjk(p)) goto 999
        call qqb_Waa(p,msq)
      elseif (kcase==kZgajet) then
        if (decayChannel() == decayQuarks) then
          call qqb_zaj_vdecay(p,msq)
        else
          call set_anomcoup(p)
          call qqb_zaj(p,msq)
        endif
      elseif (kcase==kZ2gajt) then
        call qqb_zaa_g(p,msq)
      elseif (kcase==kZga2jt) then
        call set_anomcoup(p)
        call qqb_zaj_g(p,msq)
!        call qqb_zaj_g_mad(p,msq1)
!        call qqb_zaj_g(p,msq)
!        do j=-4,4
!        do k=-4,4
!        if (msq(j,k) .ne. zip) write(6,*) j,k,msq(j,k),msq(j,k)/msq1(j,k)
!        enddo
!        enddo
!        stop
      elseif (kcase==kZbbmas) then
        call qqb_zbbm(p,msq)
      elseif (kcase==kZbbbar) then
        call qqb_zbb(p,msq)
      elseif (kcase==kZbbjet) then
        call qqb_zbb_g(p,msq)
      elseif (kcase==kWWqqbr) then
        call qqb_ww(p,msq)
      elseif (kcase==kWWnpol) then
        call qqb_ww_unpol(p,msq)
      elseif (kcase==kWW_jet) then
        call qqb_ww_g(p,msq)
      elseif (kcase==kWpWp2j) then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpWp3j) then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpmZjj) then
        call qqb_WZjj(p,msq)
      elseif (kcase==kWpmZbj) then
        call qqb_WZbj(p,msq)
      elseif (kcase==kWpmZbb) then
        call qqb_WZbb(p,msq)
      elseif (kcase==kWZbbar) then
        call qqb_wz(p,msq)


      elseif (kcase==kWW_jet) then
        call qqb_ww_g(p,msq)
c--- Check of gvec routine
c       n(1)=1._dp
c       n(2)=0._dp
c       n(3)=0._dp
c       n(4)=0._dp
c       call qqb_ww_gvec(p,n,2,msqn)
c       n(1)=0._dp
c       n(2)=1._dp
c       n(3)=0._dp
c       n(4)=0._dp
c       call qqb_ww_gvec(p,n,2,msqmad)
c--- polarization vectors in general (for p7)
c       n(1)=p(7,2)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(2)=-p(7,1)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(3)=0._dp
c       n(4)=0._dp       
c       call qqb_ww_gvec(p,n,7,msqn)
c       n(1)=p(7,1)*p(7,3)/p(7,4)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(2)=p(7,2)*p(7,3)/p(7,4)/sqrt(p(7,1)**2+p(7,2)**2)
c       n(3)=-sqrt(p(7,1)**2+p(7,2)**2)/p(7,4)
c       n(4)=0._dp       
c       call qqb_ww_gvec(p,n,7,msqmad)
c       do j=-2,2
c       do k=-2,2
c       if (msq(j,k) .ne. 0._dp) write(6,'(2i3,2e18.8,f16.9)') 
c     &    j,k,msq(j,k),msqmad(j,k)+msqn(j,k),
c     &    msq(j,k)/(msqmad(j,k)+msqn(j,k))
c       enddo
c       enddo
c       pause
c--- Madgraph check
c        call qqb_ww_g_mad(p,msqmad)
c       do j=-4,4
c       do k=-4,4
c       if (msq(j,k) .ne. 0._dp)
c     &    write(6,'(2i3,2e18.8,f16.9)')
c     &    j,k,msq(j,k),msqmad(j,k),msq(j,k)/msqmad(j,k)
c       enddo
c       enddo
c       pause

c      elseif (kcase==kWW2jet) then
c        call qqb_wwg_g(p,msq)

c        call qqb_ww_gg_mad(p,msqmad)
c       do j=-4,4
c       do k=-4,4
c       if (msqmad(j,k) .ne. 0._dp) write(6,'(2i3,2e18.8,f16.9)') 
c     &    j,k,msq(j,k),msqmad(j,k),msq(j,k)/msqmad(j,k)
c       enddo
c       enddo
c       pause
      elseif (kcase==kWpWp2j) then
        call qqb_wpwp_qqb(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpWp3j) then
        call qqb_wpwp_qqb_g(p(1:12,:),msq(-5:5,-5:5))
      elseif (kcase==kWpmZjj) then
        call qqb_WZjj(p,msq)
      elseif (kcase==kWpmZbj) then
        call qqb_WZbj(p,msq)
      elseif (kcase==kWpmZbb) then
        call qqb_WZbb(p,msq)
      elseif (kcase==kWZbbar) then
        call qqb_wz(p,msq)
      elseif (kcase==kZZlept) then
        call qqb_zz(p,msq)
      elseif (kcase==kZZ_jet) then
        call qqb_zz_g(p,msq)
      elseif (kcase==kWHbbar) then
        call qqb_wh(p,msq)
      elseif (kcase==kWH1jet) then
        call qqb_WH1jet(p,msq)
      elseif (kcase==ktwojet) then
c         reweight=s(1,2)**2*reweight
        call qqb_twojet(p,msq)
        if (kewcorr /= knone) then
          msq_noew=msq
          if     (kewcorr == ksudakov) then
            call qqb_twojet_ew_sudakov(p,msq)
          elseif (kewcorr == kexact) then
            stop 'not implemented yet'
          endif
        endif
      elseif (kcase==ktwo_ew) then
c         reweight=s(1,2)**2*reweight
        call qqb_twojet_ew(p,msq)
        call qqb_twojet(p,msq_noew)
c      elseif (kcase==kthrjet) then
c        call qqb_3jet(p,msq)
      elseif (kcase==kdirgam) then
        call qqb_dirgam(p,msq)
      elseif (kcase==khflgam) then
        call qqb_hflgam(p,msq)
      elseif (kcase==kgamgam) then
        call qqb_gamgam(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        call qqb_gamgam_mad(p,msq)
c        write(6,*) 'gamgam',msq(1,-1),msq(2,-2),msq(-1,1),msq(-2,2) 
c        pause
      elseif (kcase==kgg2gam) then
         if(toploopgaga) then
            msq(:,:)=zip
!            call gg_2gam(p,msq)
            call gggaga_mt(p,msq(0,0))
         else
            call gg_2gam(p,msq)
         endif
      elseif (kcase==kgmgmjt) then
c      call checkgvec(+2, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 0,2,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0,-1,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 0, 2,1,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec( 1,-1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
c      call checkgvec(-1, 1,5,p,qqb_gamgam_g,qqb_gmgmjt_gvec)
!        call qqb_gamgam_g(p,msq)
         call qqb_gmgmjt(p,msq)
      elseif (kcase==ktrigam) then
        call qqb_trigam(p,msq)
      elseif (kcase==kfourga) then 
         call qqb_fourgam(p,msq)
      elseif (kcase==kgamjet) then
        call qqb_dirgam_g(p,msq)
      elseif (kcase==kWH__WW) then
        call qqb_wh_ww(p,msq)
      elseif (kcase==kWH__ZZ) then
        call qqb_wh_zz(p,msq)
      elseif (kcase==kWHgaga) then
        call qqb_wh_gaga(p,msq) 
      elseif (kcase==kZHbbar) then
        call qqb_zh(p,msq)
      elseif (kcase==kZH1jet) then
        call qqb_ZH1jet(p,msq)
      elseif (kcase==kZH__WW) then
        call qqb_zh_ww(p,msq)
      elseif (kcase==kZH__ZZ) then
        call qqb_zh_zz(p,msq)
      elseif (kcase==kZHgaga) then
        call qqb_zh_gaga(p,msq)
      elseif (kcase==kggfus0) then
        call gg_h(p,msq)
      elseif (kcase==kHigaga) then
        call gg_hgamgam(p,msq)
      elseif (kcase==kHi_Zga) then
        call gg_hzgam(p,msq)
      elseif (kcase==kHi_Zaj) then
!        call checkgvec(+2, 0,2,p,gg_hzgamg,gg_hg_zgam_gvec)
!        call checkgvec(-1, 0,2,p,gg_hzgamg,gg_hg_zgam_gvec)
!        call checkgvec( 0,-1,1,p,gg_hzgamg,gg_hg_zgam_gvec)
!        call checkgvec( 0, 2,1,p,gg_hzgamg,gg_hg_zgam_gvec)
!        call checkgvec( 1,-1,6,p,gg_hzgamg,gg_hg_zgam_gvec)
!        call checkgvec(-1, 1,6,p,gg_hzgamg,gg_hg_zgam_gvec)
        call gg_hzgamg(p,msq)
      elseif (kcase==kHWW_4l) then
        call qqb_hww(p,msq)
      elseif (kcase==kHWW2lq) then
        call qqb_hww(p,msq)
      elseif (kcase==kHWW_tb) then
        call qqb_hww_tb(p,msq)
      elseif ((kcase==kHWWint) .or. (kcase==kHWWHpi)
     &   .or. (kcase==kggWW4l)) then
        call gg_ww_int(p,msq)
      elseif (kcase==kggWWbx) then
        msq(:,:)=zip
        call gg_WW(p,msq(0,0))
      elseif (kcase==kHZZ_4l) then
        call qqb_hzz(p,msq)
      elseif (kcase==kHZZ_tb) then
        call gg_hzz_tb(p,msq)
      elseif (kcase==kHVV_tb) then
        call gg_hvv_tb(p,msq)
      elseif (kcase==kggVV4l) then
        call gg_VV_all(p,msq)
      elseif (kcase==kggVVbx) then
        msq(:,:)=0._dp
        call gg_VV(p,msq(0,0))
      elseif (kcase==kHZZint) then
        call gg_zz_int(p,msq)
      elseif (kcase==kHZZHpi) then
        call gg_zz_Hpi(p,msq)
      elseif (kcase==kggZZ4l) then
        call gg_zz_all(p,msq)
      elseif (kcase==kggZZbx) then
        msq(:,:)=0._dp
        call gg_ZZ(p,msq(0,0))
      elseif (kcase==kHZZqgI) then 
         call qg_Hint_ZZ(p,msq)
      elseif (kcase==kH_1jet) then
        call qqb_hg(p,msq)
      elseif (kcase==kttZbbl) then
        call qqbZtt(p,msq)
      elseif ((kcase==ktt_bbl) .or. (kcase==ktt_bbh)) then
        call qqb_QQbdk(p,msq)
      elseif (kcase==ktt_bbu) then
        call qqb_QQbdku(p,msq)
      elseif (kcase==kqq_ttg) then
       call qqb_QQbdk_g(p,msq)
      elseif (kcase==ktt_tot) then
        call qqb_QQb(p,msq)
        if (kewcorr /= knone) then
          msq_noew=msq
          if     (kewcorr == ksudakov) then
            call qqb_QQb_ew_sudakov(p,msq)
          elseif (kewcorr == kexact) then
            stop 'not implemented yet'
          endif
        endif
      elseif (kcase==kbb_tot) then
        call qqb_QQb(p,msq)
      elseif (kcase==kcc_tot) then
        call qqb_QQb(p,msq)
      elseif (kcase==ktt_glu) then
         call qqb_QQb_g(p,msq)
      elseif (kcase==ktopanom) then
        bbfrac = 0._dp
        call singletop2_tree(p,msq)
      elseif (kcase==kbq_tpq) then
        call bq_tpq(p,msq)
      elseif (kcase==kttdkay) then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (kcase==kt_bbar) then
      call qqb_tbbdk(p,msq)
      elseif (kcase==ktdecay) then
        write(6,*) 'This process is not a leading order contribution'
        stop
      elseif (kcase==kW_tndk) then
        call qqb_w_tndk(p,msq)
      elseif (kcase==kW_twdk) then
        call qqb_w_twdk(p,msq)
      elseif (kcase==kWtbwdk) then
        call qqb_wtbwdk(p,msq)
      elseif (kcase==kWtbndk) then
        call qqb_wtbndk(p,msq)
      elseif (kcase==ktottth) then
        call qqb_tottth(p,msq)
      elseif (kcase==kqq_tth) then
        call qqb_tth(p,msq)
      elseif (kcase==ktth_ww) then
        call qqb_tth(p,msq)
      elseif (kcase==kqq_ttz) then
        call qqb_ttz(p,msq)
      elseif (kcase==kqqtthz) then
        call qqb_ttz(p,msq)
      elseif (kcase==kqq_ttw) then
        call qqb_ttw(p,msq)
      elseif (kcase==khttjet) then
        call qqb_higgs(p,msq)
      elseif (kcase==kggfus1) then
        call gg_hg(p,msq)
      elseif (kcase==khjetma) then
        call hjetmass(p,msq)
      elseif (kcase==kHgagaj) then
        call gg_hgagag(p,msq)
      elseif (kcase==kHWWjet) then
        call gg_hWWg(p,msq)
      elseif (kcase==kHZZjet) then
        call gg_hZZg(p,msq)
      elseif (kcase==kHWW2jt) then
        call gg_hWWgg(p,msq)
      elseif (kcase==kHZZ2jt) then
        call gg_hZZgg(p,msq)
      elseif (kcase==kHWW3jt) then
        call gg_hWWggg(p,msq)
      elseif (kcase==kHZZ3jt) then
        call gg_hZZggg(p,msq)
      elseif (kcase==kattjet) then
        call qqb_higgs_odd(p,msq)
      elseif (kcase==kqq_Hqq) then
        call VV_hqq(p,msq)
      elseif (kcase==kqq_Hgg) then
        call VV_Hgaga(p,msq)
      elseif (kcase==kqqHqqg) then
        call VV_hqq_g(p,msq)
      elseif (kcase==kqq_HWW) then
        call VV_HWW(p,msq)
      elseif (kcase==kqq_HZZ) then
        call VV_HZZ(p,msq)
      elseif (kcase==ktautau) then
        call qqb_tautau(p,msq)
      elseif (kcase==kqg_tbq) then
        call qg_tbq(p,msq)
c--- Check of gvec routines
c      call checkgvec(+2,0,2,p,qg_tbq,qg_tbq_gvec)
c      call checkgvec(-1,0,2,p,qg_tbq,qg_tbq_gvec)
      elseif (kcase==kqqZZqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_ZZqqstrong(p,msq)
        else
          call qq_ZZqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqWWqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WWqqstrong(p,msq)
        else
          call qq_WWqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqVVqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
c          call qq_VVqqstrong(p,msq)
          write(6,*) 'Not yet implemented'
          stop
        else
          call qq_VVqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqWWss) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WWssstrong(p,msq)
        else
          call qq_WWss(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqqWZqq) then
c        call getvbfpoint(p)
        if (VVstrong) then
          call qq_WZqqstrong(p,msq)
        else
          call qq_WZqq(p,msq)
        endif
c        call comparevbf(msq)
      elseif (kcase==kqgtbqq) then
        call qg_tbq_g(p,msq)
      elseif (kcase==k4ftwdk) then
        call qg_tbqdk(p,msq)
      elseif (kcase==k4ftjet) then
        call qg_tbqdk_g(p,msq)
      elseif (kcase==kqq_tbg) then
        call qq_tbg(p,msq)
c--- Check of gvec routines
c      call checkgvec(2,-1,5,p,qq_tbg,qq_tbg_gvec)
c      call checkgvec(-1,2,5,p,qq_tbg,qq_tbg_gvec)
      elseif (kcase==kqqtbgg) then
        call qq_tbg_g(p,msq)
      elseif (kcase==kepem3j) then
        call epem3j(p,msq)
      elseif (kcase==kgQ__ZQ) then
        call gQ_zQ(p,msq)
      elseif (kcase==kZccmas) then
        call qqb_zccm(p,msq)
      elseif (kcase==kggfus2) then
        call gg_hgg(p,msq)
      elseif (kcase==kgagajj) then
        call gg_hgg(p,msq)
      elseif (kcase==kh2jmas) then
        badpoint = .false.
        call hjetmass_r(p,msq)
        if (badpoint) then
          msq(:,:)=zip
        endif
      elseif (kcase==kggfus3) then
        call gg_hggg(p,msq)
      elseif (kcase==kW_bjet) then
        call qqb_wbjet(p,msq)
      elseif (kcase==kWcjetg) then
        call qqb_w_cjet_massless_g(p,msq)
      elseif (kcase==kZ_bjet) then
        call qqb_zbjet(p,msq)
      elseif (kcase==kZbjetg) then
        call qqb_zbjet_g(p,msq)
      elseif (kcase==kH_tjet) then
        call qq_tchan_htq(p,msq)
      elseif (kcase==kH_tdkj) then
        call qq_tchan_htq_dk(p,msq)
      elseif (kcase==kZ_tjet) then
        call qq_tchan_ztq(p,msq)
      elseif (kcase==kZ_tdkj) then
         call qq_tchan_ztq_dk(p,msq)
      elseif (kcase==kZtdk2j) then
         call qq_tchan_ztqg_dk(p,msq)
      elseif (kcase==kZt2jet) then
        call qq_tchan_ztqg(p,msq)
      elseif (kcase==kHHpair) then
        call gg_HH(p,msq)
!        call gg_HH_phase(p,msq1)
!        write(6,*) 'msq(0,0),msq1(0,0),msq1(0,0)/msq(0,0)',msq(0,0),msq1(0,0),msq1(0,0)/msq(0,0)
      elseif (kcase==kdm_jet) then 
         call qqb_dm_monojet(p,msq)
      elseif (kcase==kdm_gam) then 
         call qqb_dm_monophot(p,msq)
      elseif ( kcase==kdm2jet) then 
         call qqb_dm_monojet_g(p,msq) 
      elseif ( kcase==kdm_gaj) then 
         call qqb_dm_monophot_g(p,msq)
      elseif (kcase==kvlchk2) then
        call qqb_vol(p,msq)
        flux=one/vol(W,2)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=1._dp
        fx2(-1)=1._dp
      elseif (kcase==kvlchk3) then
        call qqb_vol(p,msq)
        flux=one/vol(W,3)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=2._dp
        fx2(-1)=2._dp
      elseif (kcase==kvlchk4) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol(W,4)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp*xx(1)
        fx2(-1)=2._dp/xx(2)
      elseif (kcase==kvlchk5) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol(W,5)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp
        fx2(-1)=4._dp
      elseif (kcase==kvlchk6) then
        call qqb_vol(p,msq)
        flux=one/vol(W,6)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=6._dp*xx(1)
        fx2(-1)=4._dp/xx(2)
      elseif (kcase==kvlchk8) then
        call qqb_vol(p,msq)
        flux=one/vol(W,8)
        bbsqmax=W
        bbsqmin=0._dp
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=6._dp/xx(1)
        fx2(-1)=6._dp/xx(2)
      elseif (kcase==kvlchkm) then
        taumin=0.0001_dp
        bbsqmax=W
        bbsqmin=0._dp
        call qqb_vol(p,msq)
        flux=one/vol_mass(mb,W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=4._dp*xx(1)
        fx2(-1)=2._dp/xx(2)      
      elseif (kcase==kvlchm3) then
        taumin=(2._dp*mt/sqrts)**2
        call qqb_vol(p,msq)
        flux=one/vol3_mass(mt,W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=2._dp
        fx2(-1)=2._dp
      elseif ((kcase==kvlchwt) .or. (kcase==kvlchwn)
     &   .or. (kcase==kvlchwg) .or. (kcase==kvlchwh)) then
        taumin=0.0001_dp
        call qqb_vol(p,msq)
        flux=one/vol_wt(W)
        do j=-nf,nf
        fx1(j)=0._dp
        fx2(j)=0._dp
        enddo
        fx1(2)=1._dp
        fx2(-1)=1._dp
      else
        write(6,*) 'Unimplemented process in lowint : kcase=',kcase
        stop 
      endif
! code to find power of alpha-s for scale variation
      if ((doScalevar .and. currentPDF == 0 .and. bin) .and. (foundpow .eqv. .false.)) then
        if (itrial == 1) then
          msqtrial=maxval(msq)
          if (msqtrial == 0) goto 999
          as=as*two
          ason2pi=ason2pi*two
          ason4pi=ason4pi*two
          gsq=gsq*two
          itrial=itrial+1
          goto 66
        endif
        msqtrial=maxval(msq)/msqtrial
        alphaspow=-1
        if (abs(msqtrial-one) < 1.e-8) alphaspow=0
        if (abs(msqtrial-two) < 1.e-8) alphaspow=1
        if (abs(msqtrial-four) < 1.e-8) alphaspow=2
        if (abs(msqtrial-eight) < 1.e-8) alphaspow=3
        if (abs(msqtrial-16._dp) < 1.e-8) alphaspow=4
        if (abs(msqtrial-32._dp) < 1.e-8) alphaspow=5
        if (alphaspow == -1) then
          write(6,*) 'Unable to determine power of alpha-s for scale variation'
          stop
        endif
        as=as/two
        ason2pi=ason2pi/two
        ason4pi=ason4pi/two
        gsq=gsq/two
        foundpow=.true.
!        write(6,*) 'Found alpha-s power: ',alphaspow
        goto 66
      endif

c--- APPLgrid - initialize array
      if (creategrid.and.bin) then
         do j=-nf,nf
            do k=-nf,nf
               weightb(j,k) = 0d0
            enddo
         enddo
         weightfactor = 1d0
      endif
c--- APPLgrid - end     
c--- do not calculate the flux if we're only checking the volume      
c      if (case(1:4) .ne. 'vlch') then      
      flux=fbGeV2/(2._dp*xx(1)*xx(2)*W)
c      endif
      
c--- initialize a PDF set here, if calculating errors
  777 continue    
      xmsq=0._dp
      xmsq_noew=0._dp

c--- calculate PDF's  
      if (((kcase==kqg_tbq) .or. (kcase==k4ftwdk))
     &     .and. (dynamicscale .eqv. .false.)) then
c--- for single top + b, make sure to use two different scales
         if (doPDFerrors .and. doScalevar) then
             error stop "simultaneous PDF and scale variation not supported for this process"
         endif

         call fdist(ih1,xx(1),facscale_H,fx1_H)
         call fdist(ih2,xx(2),facscale_H,fx2_H)
         call fdist(ih1,xx(1),facscale_L,fx1_L)
         call fdist(ih2,xx(2),facscale_L,fx2_L)

         do j=-nf,nf
           if (j == 0) then   ! heavy quark line has gluon init. state
             fx1(j)=fx1_H(j)
             fx2(j)=fx2_H(j)
           else
             fx1(j)=fx1_L(j)
             fx2(j)=fx2_L(j)
           endif
         enddo
      elseif ((kcase == kqg_tbq .or. kcase==ktopanom) .and. use_DDIS) then
          if (doPDFerrors .and. doScalevar) then
              error stop "simultaneous PDF and scale variation not supported for DDIS scales"
          endif
          ! to support kqg_tbq here, we need to set fxi_H and fxi_L arrays

          call singletop2_scale_setup(p)
          ! onheavy or onlight does not matter for born kinematics
          call fdist(ih1,xx(1),facscale_beam1_isheavy_onheavy,fxb1)
          call fdist(ih2,xx(2),facscale_beam2_islight_onheavy,fx2)
          call fdist(ih1,xx(1),facscale_beam1_islight_onheavy,fx1)
          call fdist(ih2,xx(2),facscale_beam2_isheavy_onheavy,fxb2)
      else   
c--- usual case
            if ((maxPDFsets > 0) .and. bin) then
              call fdist(ih1,xx(1),facscale,fx1)
              call fdist(ih2,xx(2),facscale,fx2)
              ! this covers the case of PDF error AND scale variation for central PDF
              if (doScalevar .and. currentPDF == 0 .and. bin) then
                  call fdist(ih1,xx(1),facscale*two,fx1up)
                  call fdist(ih2,xx(2),facscale*two,fx2up)
                  call fdist(ih1,xx(1),facscale/two,fx1dn)
                  call fdist(ih2,xx(2),facscale/two,fx2dn)
                  xmsqvar(:)=zip
              endif
            else
              call fdist(ih1,xx(1),facscale,fx1)
              call fdist(ih2,xx(2),facscale,fx2)
              if (doScalevar .and. currentPDF == 0 .and. bin) then
                call fdist(ih1,xx(1),facscale*two,fx1up)
                call fdist(ih2,xx(2),facscale*two,fx2up)
                call fdist(ih1,xx(1),facscale/two,fx1dn)
                call fdist(ih2,xx(2),facscale/two,fx2dn)
                xmsqvar(:)=zip
            endif
        endif
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) cycle
      endif

      if (gqonly) then
      if (((j==0).and.(k==0)) .or. ((j.ne.0).and.(k.ne.0))) cycle
      endif
      
      if (noglue) then 
      if ((j==0) .or. (k==0)) cycle
      endif

      if (omitgg) then 
      if ((j==0) .and. (k==0)) cycle
      endif

      if     ((kcase==kbq_tpq .or. kcase==ktopanom) .and. use_DDIS) then
c--- special case for dynamic scale in t-channel single top
        if     (abs(j) == 5) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (abs(k) == 5) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0._dp
      endif
      elseif ((kcase==kqg_tbq) .and. use_DDIS) then
c--- special case for dynamic scale in t-channel single top
        if     (j == 0) then
          xmsqjk=fxb1(j)*fx2(k)*msq(j,k)
      elseif (k == 0) then
          xmsqjk=fx1(j)*fxb2(k)*msq(j,k)
      else
        xmsqjk=0._dp
      endif
      else
c--- DEFAULT
        xmsqjk=fx1(j)*fx2(k)*msq(j,k)
        if (doScalevar .and. currentPDF == 0 .and. bin) then
          xmsqvar(1)=xmsqvar(1)+fx1up(j)*fx2up(k)*msq(j,k)
          xmsqvar(2)=xmsqvar(2)+fx1dn(j)*fx2dn(k)*msq(j,k)
        endif
        if (kewcorr /= knone) xmsqjk_noew=fx1(j)*fx2(k)*msq_noew(j,k)
      endif

      xmsq=xmsq+xmsqjk
      if (kewcorr /= knone) xmsq_noew=xmsq_noew+xmsqjk_noew
      
      if     (j > 0) then
        sgnj=+1
      elseif (j < 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k > 0) then
        sgnk=+1
      elseif (k < 0) then
        sgnk=-1
      else
        sgnk=0
      endif

c--- APPLgrid - save weight
      if (currentPDF .eq. 0 .and. creategrid .and. bin ) then
c---- print*,j,k,msq(j,k)
           weightb(j,k) =  weightb(j,k) + msq(j,k)
      endif
c--- APPLgrid - end

      enddo
      enddo ! end loop over partons


      if (currentPDF == 0) then
        lowint=flux*pswt*xmsq/BrnRat
        if (kewcorr /= knone) lowint_noew=flux*pswt*xmsq_noew/BrnRat
      endif

! compute weights for scale variation
      if (doScalevar .and. currentPDF == 0 .and. bin) then
        if (abs(xmsq) > zip) then
          scalereweight(1)=(alphas(scale*two,amz,nlooprun)/as)**alphaspow
          scalereweight(2)=(alphas(scale/two,amz,nlooprun)/as)**alphaspow
          scalereweight(1)=scalereweight(1)*xmsqvar(1)/xmsq
          scalereweight(2)=scalereweight(2)*xmsqvar(2)/xmsq
          if (maxscalevar == 6) then
            scalereweight(3)=scalereweight(1)*xmsq/xmsqvar(1)
            scalereweight(4)=scalereweight(2)*xmsq/xmsqvar(2)
            scalereweight(5)=xmsqvar(1)/xmsq
            scalereweight(6)=xmsqvar(2)/xmsq
          endif
        else
          scalereweight(:)=zip
        endif
      endif
            
c--- loop over all PDF error sets, if necessary
#ifdef HAVE_LHAPDF
      if ((maxPDFsets > 0) .and. bin) then
        if (currentPDF > 0) then
            pdfreweight(currentPDF) = (lowint - flux*pswt*xmsq/BrnRat)*wgt
        endif

        currentPDF=currentPDF+1
        if (currentPDF <= maxPDFsets) then
            if (doPDFAlphas) then
                ! only if as(mz) has changed recompute matrix element
                if ( abs(getalphas(zmass) - amz) > 1d-5 ) then
                    goto 771
                endif
            endif
            goto 777
        endif
      endif
#endif

      call getptildejet(0,pjet)
      
      call dotem(nvec,pjet,s)

      val=lowint*wgt
      val2=val**2
      if (ieee_is_nan(val)) then
         write(6,*) 'lowint val = ',val
         write(6,*) 'Discarding point with random variables',r
         lowint=zip
         val=zip
         goto 999
      endif

      if ((abs(val) > wtmax)) then
        wtmax=abs(val)
      endif

      if (bin) then
c     APPLgrid - multiply by totalFactor
            if (creategrid) then ! P.S. scale with factor
               psCR = 1d0
               if ( (kcase==ktt_tot)
     &         .or. (kcase==kbb_tot)
     &         .or. (kcase==kcc_tot) 
     &         .or. (kcase==ktt_bbl)
     &         .or. (kcase==ktt_ldk)
     &         .or. (kcase==ktt_bbu)
     &         .or. (kcase==ktt_udk)
     &         .or. (kcase==ktt_bbh)
     &         .or. (kcase==ktt_hdk)
     &         .or. (kcase==ktthWdk)
     &         .or. (kcase==kqq_ttg) ) then
                  psCR = (1d0/ason2pi)**2
               elseif ( (kcase==kW_cjet)) then
                  psCR = (1d0/ason2pi)
               endif
               do j=-nflav,nflav
                  do k=-nflav,nflav
                    weightb(j,k)=weightb(j,k)*psCR
                 enddo
              enddo           
              contrib      = 100
              weightfactor = flux*pswt*wgt/BrnRat
              ag_xx1       = xx(1)
              ag_xx2       = xx(2)
              ag_scale     = facscale
              refwt        = val
              refwt2       = val2
C     print*,"  *******************************************"
C     print*, "meWeightFactor = ", weightfactor,
C     *             " me(2,-1) = " ,  weightb(2 ,-1) ," ", msq(2,-1),
C     *             " me(-1,2) = " ,  weightb(-1 ,2) ," ", msq(-1,2),
C     *             " me(1,-1) = " ,  weightb(1 ,-1) ," ", msq(1,-1),
C     *             " me(-2,2) = " ,  weightb(-2 ,2) ," ", msq(-2,2)
C     print*, " x1 = ",xx(1)," x2 = ",xx(2)," sca = ",facscale
C     print *, "rewt = ", refwt
C     print*,"  *********************************************"
C     flush(6)
              
           endif
c---  APPLgrid - end
c--- for EW corrections, make additional weight available inside common block
        if (kewcorr /= knone) then
          wt_noew=lowint_noew*wgt
        endif
        call nplotter(pjet,val,val2,0)
      endif

      if (includeTaucutgrid(0) .eqv. .false.) then
          lowint = 0._dp
      endif

      if (enable_reweight_user) then
          lowint = lowint * reweight_user(pjet)
      endif

      return

 999  continue
      lowint=0._dp
      
      return
      end


