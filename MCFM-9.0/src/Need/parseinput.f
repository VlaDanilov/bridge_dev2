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
 
      module parseinput
          use m_config
          use types
          implicit none

          type(CFG_T), public, save :: cfg

      public

      contains

      subroutine default_config()
          implicit none

          real(dp), parameter :: sqrts = 14000d0

          call cfg_add(cfg, "mcfm_version", "9.0", "MCFM version number")
          call cfg_add(cfg, "writerefs", .true., "write references")

          ! [histogram]
          call cfg_add(cfg, "histogram%writetop", .false., "write top-drawer histograms")
          call cfg_add(cfg, "histogram%writetxt", .false., "write raw table file for each histogram")

          ! [general]
          call cfg_add(cfg, "general%nproc", 1, "process number")
          call cfg_add(cfg, "general%part", "nlo", "part: lo, nlo, nlocoeff, nnlocoeff")
          call cfg_add(cfg, "general%runstring", "run", "string identifying the run")
          call cfg_add(cfg, "general%sqrts", sqrts, "center of mass energy")
          call cfg_add(cfg, "general%ih1", +1, "ih1: +1 for proton, -1 for antiproton")
          call cfg_add(cfg, "general%ih2", +1, "ih2: +1 for proton, -1 for antiproton")
          call cfg_add(cfg, "general%zerowidth", .false., "use zero width approx. (proc. dependent)")
          call cfg_add(cfg, "general%removebr", .false., "remove decay branching ratio (proc. dependent)")
          call cfg_add(cfg, "general%ewcorr", "none", "electroweak corrections: none, sudakov or exact")

          ! [nnlo]
          call cfg_add(cfg, "nnlo%dynamictau", .true., "event-by-event tau definition")
          call cfg_add(cfg, "nnlo%tcutarray", [real(dp) ::], "optional array "//
     &        "of taucut values that should be sampled on the fly in addition", dynamic_size=.true.)

          ! additional optional settings in [general]
          call cfg_add(cfg, "general%nevtrequested", 0, "")
          call cfg_add(cfg, "general%vdecayid", .false., "")
          call cfg_add(cfg, "general%v34id", "", "")
          call cfg_add(cfg, "general%v56id", "", "")

          ! [masses]
          call cfg_add(cfg, "masses%hmass", 125d0, "Higgs mass")
          call cfg_add(cfg, "masses%mt", 173d0, "Top-quark mass")
          call cfg_add(cfg, "masses%mb", 4.66d0, "Bottom-quark mass")
          call cfg_add(cfg, "masses%mc", 1.275d0, "Charm-quark mass")

          ! [scales]
          call cfg_add(cfg, "scales%renscale", 1d0, "Renormalization scale")
          call cfg_add(cfg, "scales%facscale", 1d0, "Factorization scale")
          call cfg_add(cfg, "scales%dynamicscale", "m(34)", "Dynamic scale")
          call cfg_add(cfg, "scales%doscalevar", .false., "perform scale variation?")
          call cfg_add(cfg, "scales%maxscalevar", 6, "can be 2 or 6 for 2 or 6-point variation")

          ! [integration]
          call cfg_add(cfg, "integration%usesobol", .true., "Low Sobol low discrepancy sequence")
          call cfg_add(cfg, "integration%seed", 0, "Seed to use for pseudo random number integration")
          call cfg_add(cfg, "integration%precisiongoal", 0.002d0, "Relative precision goal")
          call cfg_add(cfg, "integration%readin", .false., "Resume from previous integration snapshot")
          call cfg_add(cfg, "integration%writeintermediate", .true.,
     &              "write histograms and results after each vegas iteration")
          call cfg_add(cfg, "integration%warmupprecisiongoal", 0.25d0, 
     &              "Relative precision goal for each contribution during warmup")
          call cfg_add(cfg, "integration%warmupchisqgoal", 2.5d0,
     &              "Chisq/it goal for each contribution during warmup")

          ! [pdf]
          call cfg_add(cfg, "pdf%pdlabel", "CT14.NN", "PDF label for internal routines")

          ! [lhapdf]
          call cfg_add(cfg, "lhapdf%lhapdfset", ["CT14nnlo"], "LHAPDF PDF label", dynamic_size=.true.)
          call cfg_add(cfg, "lhapdf%lhapdfmember", [0], "LHAPDF PDF member number, -1 for PDF uncertainties",
     &                          dynamic_size=.true.)
          call cfg_add(cfg, "lhapdf%dopdferrors", .false., "calculate PDF uncertainties.")

          ! [basicjets]
          call cfg_add(cfg, "basicjets%inclusive", .true., "jet-inclusive cross-section")
          call cfg_add(cfg, "basicjets%algorithm", "ankt", "jet-algorithm: ankt, ktal, cone, hqrk, none")
          call cfg_add(cfg, "basicjets%ptjetmin", 30d0, "minimum jet pT")
          call cfg_add(cfg, "basicjets%etajetmax", 2.4d0, "maximum jet rapidity, absolute value")
          call cfg_add(cfg, "basicjets%Rcutjet", 0.5d0, "minimum jet separation in R")
          call cfg_add(cfg, "basicjets%userap", .true. , "use rapidity not pseudorapidity for jets")

          ! [masscuts]
          call cfg_add(cfg, "masscuts%m34min", 0d0, "minimum mass of 3-4 system")
          call cfg_add(cfg, "masscuts%m56min", 0d0, "minimum mass of 5-6 system")
          call cfg_add(cfg, "masscuts%m3456min", 0d0, "minimum mass of 3-4-5-6 system")
          ! max values are handled below

          ! [cuts]
          call cfg_add(cfg, "cuts%makecuts", .true., "enable these additional cuts")

          call cfg_add(cfg, "cuts%ptleptmin", 20d0, "minimum lepton pT")
          call cfg_add(cfg, "cuts%etaleptmax", 2.4d0, "maximum lepton rapidity, absolute value")
          call cfg_add(cfg, "cuts%etaleptveto", [0d0, 0d0], "lepton rapidity veto")
          call cfg_add(cfg, "cuts%ptminmiss", 30d0, "minimum missing pt")

          call cfg_add(cfg, "cuts%ptlept2min", 20d0, "additional leptons minimum pt")
          call cfg_add(cfg, "cuts%etalept2max", 2.4d0, "maximum rapidity of additional leptons, absolute value")
          call cfg_add(cfg, "cuts%etalept2veto", [0d0, 0d0], "additional lepton rapidity veto")

          call cfg_add(cfg, "cuts%m34transmin", 0d0, "minimum (3,4) transverse mass")
          call cfg_add(cfg, "cuts%Rjlmin", 0d0, "R(jet,lept)_min")
          call cfg_add(cfg, "cuts%Rllmin", 0d0, "R(lept,lept)_min")
          call cfg_add(cfg, "cuts%delyjjmin", 0d0, "Delta eta(jet,jet)_min")
          call cfg_add(cfg, "cuts%jetsopphem", .false., "force jets to be in opposite hemispheres")
          call cfg_add(cfg, "cuts%lbjscheme", 0, "lepbtwnjets scheme")
          call cfg_add(cfg, "cuts%ptbjetmin", 0d0, "b-jet minimum pT; can also add ptbjetmax")
          call cfg_add(cfg, "cuts%etabjetmax", 100d0, "b-jet maximum rapidity, absolute value; can also add ptbjetmin")

          ! [photon]
          call cfg_add(cfg, "photon%fragmentation", .false., "fragmentation included")
          call cfg_add(cfg, "photon%fragmentation_set", "GdRG__LO", "fragmentation set")
          call cfg_add(cfg, "photon%fragmentation_scale", 1d0, "fragmentation scale")

          call cfg_add(cfg, "photon%gammptmin", 40d0, "minimum photon pT; can also add gammptmax")
          call cfg_add(cfg, "photon%gammrapmax", 2.5d0, "maximum photon rapidity; can also add gammrapmin")
          call cfg_add(cfg, "photon%gammpt2", 25d0, "second photon minimum pT")
          call cfg_add(cfg, "photon%gammpt3", 25d0, "third photon minimum pT")
          call cfg_add(cfg, "photon%Rgalmin", 0d0, "R(photon,lepton)_min")
          call cfg_add(cfg, "photon%Rgagamin", 0.4d0, "R(photon,photon)_min")
          call cfg_add(cfg, "photon%Rgajetmin", 0d0, "R(photon,jet)_min")
          call cfg_add(cfg, "photon%cone_ang", 0.4d0, "cone size for isolation")
          call cfg_add(cfg, "photon%epsilon_h", 0.5d0, "epsilon_h, energy fraction for isolation")
          call cfg_add(cfg, "photon%n_pow", 1d0, "n_pow, exponent for smooth-cone isolation")

          ! [wz2jet]
          call cfg_add(cfg, "wz2jet%qflag", .true., "")
          call cfg_add(cfg, "wz2jet%gflag", .true., "")

          ! [hjetmass]
          call cfg_add(cfg, "hjetmass%mtex", 0, "controls approximation for 2-loop virtual corrections")

          ! [anom_higgs]
          call cfg_add(cfg, "anom_higgs%hwidth_ratio", 1d0, "Gamma_H / Gamma_H(SM)")
          call cfg_add(cfg, "anom_higgs%cttH", 1d0, "cttH")
          call cfg_add(cfg, "anom_higgs%cWWH", 1d0, "cWWH")

          ! [anom_wz]
          call cfg_add(cfg, "anom_wz%enable", .false., "enable anomalous W/Z couplings")
          call cfg_add(cfg, "anom_wz%delg1_z", 0d0, "Delta g1(Z)")
          call cfg_add(cfg, "anom_wz%delk_z", 0d0, "Delta K(Z)")
          call cfg_add(cfg, "anom_wz%delk_g", 0d0, "Delta K(gamma)")
          call cfg_add(cfg, "anom_wz%lambda_z", 0d0, "Lambda(Z)")
          call cfg_add(cfg, "anom_wz%lambda_g", 0d0, "Lambda(gamma)")
          call cfg_add(cfg, "anom_wz%h1Z", 0d0, "h1(Z)")
          call cfg_add(cfg, "anom_wz%h1gam", 0d0, "h1(gamma)")
          call cfg_add(cfg, "anom_wz%h2Z", 0d0, "h2(Z)")
          call cfg_add(cfg, "anom_wz%h2gam", 0d0, "h2(gamma)")
          call cfg_add(cfg, "anom_wz%h3Z", 0d0, "h3(Z)")
          call cfg_add(cfg, "anom_wz%h3gam", 0d0, "h3(gamma)")
          call cfg_add(cfg, "anom_wz%h4Z", 0d0, "h4(Z)")
          call cfg_add(cfg, "anom_wz%h4gam", 0d0, "h4(gamma)")
          call cfg_add(cfg, "anom_wz%tevscale", 2d0, "Form-factor scale, in TeV")

          ! [singletop]

          call cfg_add(cfg, "singletop%c_phiq", 0d0, "C_phiq (O1), real-valued")
          call cfg_add(cfg, "singletop%c_phiphi", [0d0,0d0], "C_phiphi (O2), real and imaginary part")
          call cfg_add(cfg, "singletop%c_tw", [0d0,0d0], "C_tW (O3), real and imaginary part")
          call cfg_add(cfg, "singletop%c_bw", [0d0,0d0], "C_bW (O4), real and imaginary part")
          call cfg_add(cfg, "singletop%c_tg", [0d0,0d0], "C_tG (O6), real and imaginary part")
          call cfg_add(cfg, "singletop%c_bg", [0d0,0d0], "C_bG (O7), real and imaginary part")
          call cfg_add(cfg, "singletop%lambda", 1000d0, "Lambda, scale of EFT breakdown in GeV")
          call cfg_add(cfg, "singletop%enable_lambda4", .false., "enable 1/Lambda^4 contributions")
          call cfg_add(cfg, "singletop%disable_sm", .false., "disable Standard Model contributions")
          call cfg_add(cfg, "singletop%mode_anomcoup", .false., "anomalous couplings mode (only LO)")

          ! [extra]
          call cfg_add(cfg, "extra%verbose", .false., "verbose output")
          call cfg_add(cfg, "extra%new_pspace", .true., "use new_pspace")
          call cfg_add(cfg, "extra%nohistograms", .false., "generate no additional histograms")
          call cfg_add(cfg, "extra%benchmark", "none", "benchmark to run")
      end subroutine

      subroutine read_config()
            use anomcoup_tbW
            use PDFerrors, only : doPDFerrors, PDFnames, PDFmembers, numPDFsets
            use Scalevar, only : doScalevar
            use m_gencuts, only : enable_reweight_user
          implicit none
          include 'cplx.h'
          include 'maxwt.f' ! for nevtrequested
          include 'outputoptions.f' ! writetop, etc.
          include 'vdecayid.f' ! vdecayid, v34id, v56id
          include 'nproc.f'
          include 'masses.f'
          include 'runstring.f'
          include 'energy.f' ! sqrts
          include 'zerowidth.f'
          include 'removebr.f'
          include 'pdlabel.f'
          include 'lhapdf.f'
          include 'clustering.f'
          include 'jetcuts.f'
          include 'Rcut.f'
          include 'makecuts.f'
          include 'lhcb.f'
          include 'limits.f'
          include 'leptcuts.f'
          include 'frag.f'
          include 'flags.f' ! qflag, gflag
          include 'noglue.f'
          include 'realwt.f'
          include 'lc.f'
          include 'verbose.f'
          include 'debug.f'
          include 'new_pspace.f'
          include 'asymptotic.f'
          include 'alfacut.f'
          include 'betacut.f'
          include 'anomHiggs.f'
          include 'anom_higgs.f'
          include 'anomcoup.f'
          include 'scalevar.f'
          include 'ewcorr.f'
          include 'cutoff.f'
          include 'rtsmin.f'
          include 'kpart.f'
          include 'kprocess.f'
          include 'taucut.f'
          include 'ewinput.f'
          character(len=cfg_string_len) :: mcfm_version, part, dynstring

          logical :: writerefs
          common/writerefs/writerefs

          integer :: ih1,ih2
          common/density/ih1,ih2

          logical :: spira
          common/spira/spira

          integer :: nmin, nmax
          common/nmin/nmin
          common/nmax/nmax

          real(dp) :: renscale_in, facscale_in
          real(dp) :: m34min, m34max, m56min, m56max
          real(dp) :: leptveto(2), lept2veto(2)

          real(dp) :: c_phiq_in, c_phiphi_in(2), c_tw_in(2), c_bg_in(2)
          real(dp) :: c_bw_in(2), c_tg_in(2), lambda_in

          call cfg_get(cfg, "mcfm_version", mcfm_version)
          if (mcfm_version(1:3) /= "9.0") then
              write (*,*) "Unsupported input file format version ", mcfm_version(1:3)
              stop
          endif

c [general]
          call cfg_get(cfg, "writerefs", writerefs)
          call cfg_get(cfg, "general%nevtrequested", nevtrequested)
          call cfg_get(cfg, "general%vdecayid", vdecayid)
          call cfg_get(cfg, "general%v34id", v34id)
          call cfg_get(cfg, "general%v56id", v56id)

          call cfg_get(cfg, "general%nproc", nproc)
          call cfg_get(cfg, "general%part", part)
          call parse_part(part)
          call cfg_get(cfg, "general%runstring", runstring)
          call cfg_get(cfg, "general%sqrts", sqrts)
          call cfg_get(cfg, "general%ih1", ih1)
          call cfg_get(cfg, "general%ih2", ih2)
          call cfg_get(cfg, "general%zerowidth", zerowidth)
          call cfg_get(cfg, "general%removebr", removebr)

c [nnlo]
          call cfg_get(cfg, "nnlo%dynamictau", dynamictau)
          call parse_taucut()
          call parse_ewcorr()

c [histogram]
          call cfg_get(cfg, "histogram%writetop", writetop)
          call cfg_get(cfg, "histogram%writetxt", writeroot)

c [masses]
          call cfg_get(cfg, "masses%hmass", hmass)
          call cfg_get_add(cfg, "masses%wmass", wmass_inp, wmass_inp, "W mass")
          call cfg_get_add(cfg, "masses%zmass", zmass_inp, zmass_inp, "Z mass")
          call cfg_get(cfg, "masses%mt", mt)
          call cfg_get(cfg, "masses%mb", mb)
          call cfg_get(cfg, "masses%mc", mc)

          if (abs(mb) > 1.e-8_dp) then
            mbsq=mb**2
          else
            mbsq=4.75_dp**2
          endif
          if (abs(mc) > 1.e-8_dp) then
            mcsq=mc**2
          else
            mcsq=1.5_dp**2
          endif


c [pdf]
          call cfg_get(cfg, "pdf%pdlabel", pdlabel)

c [lhapdf]
          if (cfg_var_size(cfg, "lhapdf%lhapdfset") > 1) then
              numpdfsets = cfg_var_size(cfg, "lhapdf%lhapdfset") 
              if (numpdfsets /= cfg_var_size(cfg, "lhapdf%lhapdfmember")) then
                  error stop "Please specify a pdf member for each pdf set"
              endif
              write (*,*) "Running with multiple PDF sets!"
          else
              numpdfsets = 1
          endif

          allocate(PDFnames(numpdfsets))
          allocate(PDFmembers(numpdfsets))
          call cfg_get(cfg, "lhapdf%lhapdfset", PDFnames)
          call cfg_get(cfg, "lhapdf%lhapdfmember", PDFmembers)
          call cfg_get(cfg, "lhapdf%dopdferrors", doPDFerrors)


#ifndef HAVE_LHAPDF
          if (doPDFerrors) then
              error stop "LHAPDF is required for PDF uncertainties!"
          endif
#endif

c [scales]
          ! note: setup ewcorr before
          call parse_scales()

c [basicjets]
          call cfg_get(cfg, "basicjets%inclusive", inclusive)
          call parse_jetalgo()
          call cfg_get(cfg, "basicjets%ptjetmin", ptjetmin)
          call cfg_get_add(cfg, "basicjets%ptjetmax", ptjetmax, sqrts, "")
          call cfg_get(cfg, "basicjets%etajetmax", etajetmax)
          call cfg_get_add(cfg, "basicjets%etajetmin", etajetmin, 0d0, "")
          call cfg_get(cfg, "basicjets%Rcutjet", Rcut)

         if ((etajetmin < 0._dp) .or. (etajetmax < 0._dp)) then
            write(6,*) 'etajetmin and etajetmax are absolute values,'
            write(6,*) ' please reset to a positive value.'
            stop
         endif

c [masscuts]
          call cfg_get(cfg, "masscuts%m34min", m34min)
          wsqmin=m34min**2
          call cfg_get_add(cfg, "masscuts%m34max", m34max, sqrts, "")
          if (m34max > sqrts*0.9999_dp) m34max = sqrts*0.9999_dp
          wsqmax=m34max**2

          call cfg_get(cfg, "masscuts%m56min", m56min)
          bbsqmin=m56min**2
          call cfg_get_add(cfg, "masscuts%m56max", m56max, sqrts, "")
          if (m56max > sqrts*0.9999_dp) m56max = sqrts*0.9999_dp
          bbsqmax=m56max**2

          call cfg_get(cfg, "masscuts%m3456min", m3456min)
          call cfg_get_add(cfg, "masscuts%m3456max", m3456max, sqrts, "")
          if (m3456max > sqrts) m3456max = sqrts

c [cuts]
          call cfg_get(cfg, "cuts%makecuts", makecuts)

          call cfg_get(cfg, "cuts%ptleptmin", leptptmin)
          call cfg_get_add(cfg, "cuts%ptleptmax", leptptmax, sqrts, "")
          call cfg_get(cfg, "cuts%etaleptmax", leptrapmax)
          call cfg_get_add(cfg, "cuts%etaleptmin", leptrapmin, 0d0, "")
          call cfg_get(cfg, "cuts%etaleptveto", leptveto)
          leptveto1min = leptveto(1)
          leptveto1max = leptveto(2)
          call cfg_get(cfg, "cuts%ptminmiss", misspt)

          call cfg_get(cfg, "cuts%ptlept2min", leptpt2min)
          call cfg_get_add(cfg, "cuts%ptlept2max", leptpt2max, sqrts, "")
          call cfg_get(cfg, "cuts%etalept2max", leptrap2max)
          call cfg_get_add(cfg, "cuts%etalept2min", leptrap2min, 0d0, "")

          call cfg_get(cfg, "cuts%etalept2veto", lept2veto)
          leptveto2min = lept2veto(1)
          leptveto2max = lept2veto(2)

          call cfg_get(cfg, "cuts%m34transmin", mtrans34cut)
          call cfg_get(cfg, "cuts%Rjlmin", Rjlmin)
          call cfg_get(cfg, "cuts%Rllmin", Rllmin)
          call cfg_get(cfg, "cuts%delyjjmin", delyjjmin)
          call cfg_get(cfg, "cuts%jetsopphem", jetsopphem)
          call cfg_get(cfg, "cuts%lbjscheme", lbjscheme)
          call cfg_get(cfg, "cuts%ptbjetmin", ptbjetmin)
          call cfg_get_add(cfg, "cuts%ptbjetmax", ptbjetmax, sqrts, "")
          call cfg_get(cfg, "cuts%etabjetmax", etabjetmax)
          call cfg_get_add(cfg, "cuts%etabjetmin", etabjetmin, 0d0, "")

c [photon]
          call cfg_get(cfg, "photon%fragmentation", frag)
          call cfg_get(cfg, "photon%fragmentation_set", fragset)
          call cfg_get(cfg, "photon%fragmentation_scale", frag_scale)
          frag_scalestart = frag_scale

          call cfg_get(cfg, "photon%gammptmin", gammptmin)
          call cfg_get_add(cfg, "photon%gammptmax", gammptmax, sqrts, "")

          call cfg_get(cfg, "photon%gammrapmax", gammrapmax)
          call cfg_get_add(cfg, "photon%gammrapmin", gammrapmin, 0d0, "")

          call cfg_get(cfg, "photon%gammpt2", gammpt2)
          call cfg_get(cfg, "photon%gammpt3", gammpt3)
          call cfg_get(cfg, "photon%Rgalmin", Rgalmin)
          call cfg_get(cfg, "photon%Rgagamin", Rgagamin)
          call cfg_get(cfg, "photon%Rgajetmin", Rgajetmin)
          call cfg_get(cfg, "photon%cone_ang", cone_ang)
          call cfg_get(cfg, "photon%epsilon_h", epsilon_h)
          call cfg_get(cfg, "photon%n_pow", n_pow)

c [lhcb]
          if (cfg_var_configadded(cfg, "lhcb%cut_mode")) then
              ! assume that block with lhcb config settings exists
              call cfg_get(cfg, "lhcb%cut_mode", cut_mode)
              call cfg_get(cfg, "lhcb%dir_mode", dir_mode)
              call cfg_get(cfg, "lhcb%nl_min", nl_min)
              call cfg_get(cfg, "lhcb%nj_min", nj_min)
              call cfg_get(cfg, "lhcb%nb_min", nb_min)

              if (cut_mode > 0) then
                  call lhcb_config()
              endif
          endif

c [integration]
          call parse_seed()

c [wz2jet]
          call cfg_get(cfg, "wz2jet%qflag", qflag)
          call cfg_get(cfg, "wz2jet%gflag", gflag)

c [hjetmass]
          call cfg_get(cfg, "hjetmass%mtex", mtex)

c [anom_higgs]
          call cfg_get(cfg, "anom_higgs%hwidth_ratio", hwidth_ratio)
          call cfg_get(cfg, "anom_higgs%cttH", cttH)
          call cfg_get(cfg, "anom_higgs%cWWH", cWWH)

c [anom_wz]
          call cfg_get(cfg, "anom_wz%enable", anomtgc)
          call cfg_get(cfg, "anom_wz%delg1_z", delg1_z)
          call cfg_get(cfg, "anom_wz%delk_z", delk_z)
          call cfg_get(cfg, "anom_wz%delk_g", delk_g)
          call cfg_get(cfg, "anom_wz%lambda_z", lambda_z)
          call cfg_get(cfg, "anom_wz%lambda_g", lambda_g)
          call cfg_get(cfg, "anom_wz%h1Z", h1Z)
          call cfg_get(cfg, "anom_wz%h1gam", h1gam)
          call cfg_get(cfg, "anom_wz%h2Z", h2Z)
          call cfg_get(cfg, "anom_wz%h2gam", h2gam)
          call cfg_get(cfg, "anom_wz%h3Z", h3Z)
          call cfg_get(cfg, "anom_wz%h3gam", h3gam)
          call cfg_get(cfg, "anom_wz%h4Z", h4Z)
          call cfg_get(cfg, "anom_wz%h4gam", h4gam)
          call cfg_get(cfg, "anom_wz%tevscale", tevscale)

c [singletop]

          call cfg_get(cfg, "singletop%c_phiq", c_phiq_in)
          call anomcoup_tbW_set_c1(cmplx(c_phiq_in,0._dp,dp))

          call cfg_get(cfg, "singletop%c_phiphi", c_phiphi_in)
          call anomcoup_tbW_set_c2(cmplx(c_phiphi_in(1), c_phiphi_in(2), dp))

          call cfg_get(cfg, "singletop%c_tw", c_tw_in)
          call anomcoup_tbW_set_c3(cmplx(c_tw_in(1), c_tw_in(2), dp))

          call cfg_get(cfg, "singletop%c_bw", c_bw_in)
          call anomcoup_tbW_set_c4(cmplx(c_bw_in(1), c_bw_in(2), dp))

          call cfg_get(cfg, "singletop%c_tg", c_tg_in)
          call anomcoup_tbW_set_c6(cmplx(c_tg_in(1), c_tg_in(2), dp))

          call cfg_get(cfg, "singletop%c_bg", c_bg_in)
          call anomcoup_tbW_set_c7(cmplx(c_bg_in(1), c_bg_in(2), dp))

          call cfg_get(cfg, "singletop%lambda", lambda_in)
          call anomcoup_tbW_set_lambda(lambda_in)

          call cfg_get(cfg, "singletop%enable_lambda4",  enable_lambda4)
          if (enable_lambda4) then
              enable_lambda2 = .true.
          endif

          call cfg_get(cfg, "singletop%disable_sm", disable_sm)
          call cfg_get(cfg, "singletop%mode_anomcoup", mode_anomcoup)

c [extra]
          ! these are flags that are usually set with the runstring
          ! or are optional technical parameters that should usually not be set

          call cfg_get_add(cfg, "extra%debug", debug, .false., "debug")
          call cfg_get(cfg, "extra%verbose", verbose)

          call cfg_get_add(cfg, "extra%toponly", toponly, .false., "")
          call cfg_get(cfg, "extra%new_pspace", new_pspace)
          call cfg_get_add(cfg, "extra%spira", spira, .true., "")

          call cfg_get_add(cfg, "extra%noglue", noglue, .false., "no gluons in initial state")
          call cfg_get_add(cfg, "extra%ggonly", ggonly, .false., "only gluon-gluon initial state")
          call cfg_get_add(cfg, "extra%gqonly", gqonly, .false., "only gluon-quark initial state")
          call cfg_get_add(cfg, "extra%omitgg", omitgg, .false., "omit gluon-gluon initial state")

          nmin = 1
          nmax = 2

          call cfg_get_add(cfg, "extra%clustering", clustering, .true., "")
          call cfg_get_add(cfg, "extra%realwt", realwt, .false., "")
          call cfg_get_add(cfg, "extra%colourchoice", colourchoice, 0, "")

          call cfg_get_add(cfg, "extra%cutoff", cutoff, 1d-9, "")
          call cfg_get_add(cfg, "extra%cutoff_s", cutoff_s, 1d-5, "")
          call cfg_get_add(cfg, "extra%rtsmin", rtsmin, 1d-8, "")

          call cfg_get_add(cfg, "extra%reweight", enable_reweight_user, .false.,
     &              "enable use of reweight_user routine")

c [dipoles]
          ! these are also optional
          call cfg_get_add(cfg, "dipoles%aii", aii, 1d0, "aii")
          call cfg_get_add(cfg, "dipoles%aif", aif, 1d0, "aif")
          call cfg_get_add(cfg, "dipoles%afi", afi, 1d0, "afi")
          call cfg_get_add(cfg, "dipoles%aff", aff, 1d0, "aff")
          call cfg_get_add(cfg, "dipoles%bfi", bfi, 1d0, "bfi")
          call cfg_get_add(cfg, "dipoles%bff", bff, 1d0, "bff")

          ! possibility for setting a common value
          if (cfg_var_configadded(cfg, "dipoles%alpha")) then
              call cfg_get_add(cfg, "dipoles%alpha", aii, 1d0, "common alpha parameter")
              aif = aii
              afi = aii
              aff = aii
              bfi = aii
              bff = aii
          endif

          ! if the user sets these to 1.0, we might end up with 1d0 + eps
          if (aii > 1d0) aii = 1d0
          if (aif > 1d0) aif = 1d0
          if (afi > 1d0) afi = 1d0
          if (aff > 1d0) aff = 1d0
          if (bfi > 1d0) bfi = 1d0
          if (bff > 1d0) bff = 1d0

ccc FINAL SANITY CHECKS AND WARNINGS (AND SETTINGS)
          if ((doscalevar) .and. (kewcorr /= knone)) then
            write(6,*) 'Cannot compute EW corrections and scale variation in a single run'
            stop
          endif

          if (noglue) then
            write(6,*) 'WARNING: no gluon contribution included in PDF'
          write(6,*)
          endif
          if (ggonly) then
            write(6,*) 'WARNING: only gluon-gluon flux included'
          write(6,*)
          endif
          if (gqonly) then
            write(6,*) 'WARNING: only gluon-quark flux included'
          write(6,*)
          endif
          if (omitgg) then
            write(6,*) 'WARNING: no gluon-gluon contribution included'
          write(6,*)
          endif

          if (kewcorr == kexact) then
              cutoff = max(cutoff, 1d-5)
          endif

ccc CHOOSER, set-up the variables for the process we wish to consider
          call chooser
ccc additional settings might depend on the chooser output below (n2,n3)

          ! E-M gauge invariance requires that delg1_g=0
          delg1_g=0._dp


c--- check that we have a valid value of 'part'
          if ( (kpart.ne.klord) .and. (kpart.ne.kreal) .and.
     &         (kpart.ne.kvirt) .and. (kpart.ne.ktota) ) then
            if    ( (kpart==ktodk) .and.
     &              ((kcase==kbq_tpq) .or. (kcase==kt_bbar)
     &          .or. (kcase==kW_twdk) .or. (kcase==ktt_bbl)
     &          .or. (kcase==ktt_bbh) .or. (kcase==k4ftwdk)
     &          .or. (kcase==kHWW2lq) .or. (kcase==kqq_ttw)
     &          .or. (kcase==kWWqqbr) .or. (kcase==ktt_bbu)
     &          .or. (kcase==kWHbbar) .or. (kcase==kZHbbar)) ) then
c--- this is an allowed combination
            elseif ( (kpart==kfrag) .and.
     &            ((kcase==kWgamma)
     &        .or. (kcase==kgamgam) .or. (kcase==kgg2gam)
     &        .or. (kcase==kdirgam) .or. (kcase==kdm_gam)
     &        .or. (kcase==kgmgmjt) .or .(kcase==ktrigam)
     &        .or. (kcase==kfourga)
     &        .or. (kcase==kZ_2gam) .or. (kcase==kZgajet)
     &        .or. (kcase==kW_2gam)) ) then
c--- this is an allowed combination
            elseif ((kpart==knnlo) .or. (kpart==ksnlo)) then 
c--- this is an allowed combination
            else 
             write(6,*) 'part=',part,' is not a valid option'
             write(6,*) 'for this process number.'
             stop     
            endif
          endif


c--- check that we are not trying to calculate radiation in decay at LO
          if    ( (kpart==klord) .and.
     &           ((kcase==kWWqqdk) .or. (kcase==kHWWdkW)
     &       .or. (kcase==ktt_ldk) .or. (kcase==ktt_udk)
     &       .or. (kcase==ktt_hdk) .or. (kcase==ktthWdk)
     &       .or. (kcase==kttdkay) .or. (kcase==kWtdkay)
     &       .or. (kcase==kdk_4ft) .or. (kcase==kttwldk)
     &       .or. (kcase==kWHbbdk) .or. (kcase==kZHbbdk)) ) then
              write(6,*) 'This process number cannot be used for'
              write(6,*) 'a LO calculation.'
              stop
          endif

c--- check that EW corrections are included for this process, if required
          if (kewcorr /= knone) then
            if ( (kcase == kZ_only) .or. (kcase == ktt_tot)
     &      .or. (kcase == ktt_mix) .or. (kcase == ktwojet)
     &      .or. (kcase == ktwo_ew) ) then
              continue
            else
              write(6,*) 'EW corrections not available for this process'
              stop
            endif
          endif

          ! can we remove this additional scale choice?
          ! rather let the user supply the scale?
          call negativescalehack()


      end subroutine

      subroutine parse_jetalgo()
          implicit none
          include 'clustering.f'

          call cfg_get(cfg, "basicjets%algorithm", algorithm)

          ! Assign choice of jet algorithm to integer variable
          if     (algorithm == 'ktal') then
            jetalgorithm=kt
          elseif (algorithm == 'ankt') then
            jetalgorithm=antikt
          elseif (algorithm == 'cone') then
            jetalgorithm=Rsepcone
          elseif (algorithm == 'hqrk') then
            jetalgorithm=hqrk
          elseif (algorithm == 'none') then
            jetalgorithm=noclustering
          else
            write(6,*) 'Invalid choice of jet algorithm: should be one of'
            write(6,*) 'ktal, ankt, cone, hqrk, none'
            stop
          endif

      end subroutine

      subroutine parse_ewcorr()
          implicit none
          include 'ewcorr.f'
          include 'cutoff.f'

          call cfg_get(cfg, "general%ewcorr", ewcorr)

          if     (ewcorr == 'none') then
            kewcorr=knone
          elseif (ewcorr == 'sudakov') then
            kewcorr=ksudakov
          elseif (ewcorr == 'exact') then
            kewcorr=kexact
          else
            write(6,*) 'Unexpected EW correction in input file: ',ewcorr
            stop
          endif

      end subroutine

      subroutine parse_seed()
          use cxx11random
          use omp_lib
          use iso_c_binding, only: c_loc
          implicit none
          include 'mpicommon.f'

          integer :: seed, origSeed
          common/seedBlock/seed
!$omp threadprivate(/seedBlock/)

          real(dp), save :: realSeed
!$omp threadprivate(realSeed)

          integer, target, dimension(:), allocatable :: seeds

          integer :: tid

          call cfg_get(cfg, "integration%seed", origSeed)

          allocate(seeds(omp_get_max_threads()))

!$omp parallel do
          do tid=0,omp_get_max_threads()-1
            !random seed if started with special seed value of 0
            if (origSeed == 0) then
                call random_seed()
                call random_number(realSeed)
                seed = -(nint(realSeed * huge(0)) + rank)
                seeds(tid+1) = -seed
            else
              seed = -(origSeed + omp_get_thread_num() + rank*100)
              seeds(tid+1) = -seed
            end if
          enddo
!$omp end parallel do
     
          call cxx11_init_random(c_loc(seeds))

      end subroutine

      subroutine parse_taucut()
          use SCET
          implicit none
          include 'taucut.f'
          include 'kpart.f'
          include 'nproc.f'

          if ((kpart==knnlo) .or. (kpart==ksnlo)) then
            call setuptau()
            call cfg_get_add(cfg, "nnlo%taucut", taucut, taucut, "taucut")

            usescet=.true.
            abovecut=.false.
            if (taucut < 0) then
              write(6,*) 'Must specify taucut > 0 for SCET calculation'
              stop
            endif
          else
            usescet=.false.
            abovecut=.false.
          endif

          call cfg_get_add(cfg, "nnlo%tauboost", tauboost, .true., "")
          call cfg_get_add(cfg, "nnlo%incpowcorr", incpowcorr, .true., "")
          call cfg_get_add(cfg, "nnlo%onlypowcorr", onlypowcorr, .false., "")


          if (cfg_var_configadded(cfg, "nnlo%tcutarray") .and. usescet) then
              allocate(tcutarray(cfg_var_size(cfg, "nnlo%tcutarray")))
              call cfg_get(cfg, "nnlo%tcutarray", tcutarray)
              doMultitaucut = .true.
          elseif (usescet) then
              allocate(tcutarray(5))
              tcutarray = [taucut*2, taucut*4, taucut*8, taucut*20, taucut*40]
              doMultitaucut = .true.
          else
              allocate(tcutarray(0))
              doMultitaucut = .false.
          endif

          if (doMultitaucut) then
              if ( (nproc /= 1) .and. (nproc /= 6) .and.
     &              (nproc /= 31) .and. (nproc /= 300) .and.
     &              (nproc /= 111) .and. (nproc /= 112) .and.
     &              (nproc /= 119) .and. (nproc /= 285) .and.
     &              (nproc < 91) .and. (nproc > 110) ) then
                  error stop "this process does not work with multitaucut yet"
              endif
          endif

!$omp parallel
          allocate(scetreweight(size(tcutarray)))
!$omp end parallel

          smallestTaucut = min(minval(tcutarray),taucut)

          ! this must be initialized to .true.
          ! only maketaucut modifies this, and for a non-scet run
          ! realint depends on this being true
!$omp parallel
          includeTaucutgrid(:) = .true.
!$omp end parallel

      end subroutine

      subroutine parse_scales()
          use singletop2_scale_m, only: singletop2_scale_init
          use Scalevar
          implicit none
          include 'nf.f'
          include 'nproc.f'
          include 'initialscales.f'
          include 'facscale.f'
          include 'scale.f'
          include 'stopscales.f'
          include 'dynamicscale.f'
          include 'scalevar.f'
          include 'ewcorr.f'

          real(dp) :: renscale_in, facscale_in

          call cfg_get(cfg, "scales%renscale", scale)
          call cfg_get(cfg, "scales%facscale", facscale)

          initscale=scale
          initfacscale=facscale

          ! sets use_DDIS flag if dynscale is 'DDIS'
          call singletop2_scale_init()

c  special case for stop+b process
          initrenscale_L=0._dp
          initfacscale_L=0._dp
          initrenscale_H=0._dp
          initfacscale_H=0._dp

          if (((nproc >= 231) .and. (nproc <= 240)) .and.      
     &         (scale == 0._dp) .and. (facscale == 0._dp)) then

              call cfg_get(cfg, "scales%renscale_L", initrenscale_L)
              call cfg_get(cfg, "scales%facscale_L", initfacscale_L)
              call cfg_get(cfg, "scales%renscale_H", initrenscale_H)
              call cfg_get(cfg, "scales%facscale_H", initfacscale_H)

              renscale_L=initrenscale_L
              facscale_L=initfacscale_L
              renscale_H=initrenscale_H
              facscale_H=initfacscale_H
              scale=initrenscale_H
              facscale=initfacscale_H        
          endif

c         dynamic scale

          call cfg_get(cfg, "scales%dynamicscale", dynstring)

          call cfg_get(cfg, "scales%doscalevar", doscalevar)
          call cfg_get(cfg, "scales%maxscalevar", maxscalevar)

          if (doscalevar .and. maxscalevar /= 2 .and. maxscalevar /= 6) then
              error stop "maxscalevar must be 2 or 6 when doing scale variation"
          endif
      
c---  create logical:: variable dynamicscale for use in other routines
          if (  (dynstring == 'no') .or. (dynstring == '.false.')
     &     .or. (dynstring == 'none') ) then 
             dynamicscale=.false. 
          else
             dynamicscale=.true. 
          endif




      end subroutine

      subroutine parse_part(part)
          implicit none
          include 'kpart.f'
          character(len=*), intent(in) :: part

          coeffonly=.false.
          kpart=0
          if     ((part == 'lo') .or. (part == 'lord')) then
            kpart=klord
          elseif (part == 'virt') then
            kpart=kvirt
          elseif (part == 'real') then
            kpart=kreal
          elseif ((part == 'nlo') .or. (part == 'tota')
     &       .or. (part == 'nlocoeff') .or. (part == 'totacoeff')) then
            kpart=ktota
          elseif (part == 'frag') then
            kpart=kfrag
          elseif ((part == 'todk') .or. (part == 'nlodk')) then
            kpart=ktodk
          elseif ((part == 'snlo') .or. (part == 'scetnlo')
     &       .or. (part == 'snlocoeff') .or. (part == 'scetnlocoeff')) then
            kpart=ksnlo
          elseif ((part == 'nnlo') .or. (part == 'nnlocoeff')
     &       .or. (part == 'nnloVV') .or. (part == 'nnloVVcoeff')
     &       .or. (part == 'nnloRV') .or. (part == 'nnloRVcoeff')
     &       .or. (part == 'nnloRR') .or. (part == 'nnloRRcoeff') ) then
            kpart=knnlo
            if     ((part == 'nnloVV') .or. (part == 'nnloVVcoeff')) then
              knnlopart=knnloVV
            elseif ((part == 'nnloRV') .or. (part == 'nnloRVcoeff')) then
              knnlopart=knnloRV
            elseif ((part == 'nnloRR') .or. (part == 'nnloRRcoeff')) then
              knnlopart=knnloRR
            else
              knnlopart=0
            endif
          endif
          if (index(part,'coeff') > 0) then
            coeffonly=.true.
          endif

          origKpart = kpart

          if (kpart == 0) then
            write(6,*) 'Invalid value of part = ',part
            stop
          endif


      end subroutine

      subroutine negativescalehack()
          implicit none
          include 'breit.f'
          include 'scale.f'
          include 'facscale.f'
          include 'kprocess.f'
          include 'masses.f'
          include 'qcdcouple.f'
          include 'couple.f'
          include 'constants.f'
          include 'nlooprun.f'

          real(dp) :: alphas, factor

c--- set up the default choices of static scale, if required
      if (scale < 0._dp) then
      if     (scale == -2._dp) then
        factor=0.25_dp
      elseif (scale == -3._dp) then
        factor=0.5_dp
      elseif (scale == -4._dp) then
        factor=0.75_dp
      elseif (scale == -5._dp) then
        factor=1._dp
      elseif (scale == -6._dp) then
        factor=2._dp
      elseif (scale == -7._dp) then
        factor=4._dp
        else
        factor=1._dp
      endif        
        if ((n2+n3 .ne. 0) .or. (kcase==ktt_tot)) then
c--- special case for t-tbar production
        if (kcase==ktt_tot) then
          scale=factor*mt
c--- special cases where Higgs mass is neither mass2 nor mass3
c        elseif ((case(1:1) == 'H') .or. (kcase==kWHbbar)
c     &     .or. (kcase==kZHbbar) .or. (kcase==kqq_Hqq)) then
c          scale=factor*hmass
        else     
          scale=factor*(n2*mass2+n3*mass3)/real(n2+n3,dp)
        endif
        as=alphas(scale,amz,nlooprun)
        ason2pi=as/twopi
        ason4pi=as/fourpi
        gsq=fourpi*as
        musq=scale**2
        write(6,*)
        write(6,*)'************* Strong coupling, alpha_s  ************'
        write(6,*)'*                                                  *'
        write(6,49)'alpha_s (scale)',gsq/fourpi
        write(6,49)'alpha_s (zmass)',amz
        write(6,50)' (using ',nlooprun,'-loop running of alpha_s)'  
        write(6,*)'****************************************************'
        write(6,*)
        write(6,*)'****************************************************'
        write(6,76) scale
        write(6,*)'****************************************************'
        else
        write(6,*) 'Invalid choice of renormalization scale!'
        stop
        endif
      endif

      if (facscale < 0._dp) then
      if     (facscale == -2._dp) then
        factor=0.25_dp
      elseif (facscale == -3._dp) then
        factor=0.5_dp
      elseif (facscale == -4._dp) then
        factor=0.75_dp
      elseif (facscale == -5._dp) then
        factor=1._dp
      elseif (facscale == -6._dp) then
        factor=2._dp
      elseif (facscale == -7._dp) then
        factor=4._dp
        else
        factor=1._dp
      endif        
        if ((n2+n3 .ne. 0) .or. (kcase==ktt_tot)) then
c--- special case for t-tbar production
        if (kcase==ktt_tot) then
          facscale=factor*mt
c--- special cases where Higgs mass is neither mass2 nor mass3
c        elseif ((case(1:1) == 'H') .or. (kcase==kWHbbar)
c     &     .or. (kcase==kZHbbar) .or. (kcase==kqq_Hqq)) then
c          facscale=factor*hmass
        else     
          facscale=factor*(n2*mass2+n3*mass3)/real(n2+n3,dp)
        endif
        write(6,*)
        write(6,*)'****************************************************'
        write(6,77) facscale
        write(6,*)'****************************************************'
       else
        write(6,*) 'Invalid choice of factorization scale!'
        stop
        endif
      endif

   49 format(' *  ',a20,f12.8,16x,'*')
   50 format(' *  ',6x,a8,i1,a25,8x,'*')
   76 format(' *      Renormalization scale =',f7.2,'              *')
   77 format(' *        Factorization scale =',f7.2,'              *')
   99 format(a90)


      
      end subroutine


      end module
