# Generated automatically by createMakefiles.pl
 
SRCS =	m_alloc_mods.f90 mo_advection.f90 mo_aero_gads.f90 mo_aero_tanre.f90 \
	mo_albedo.f90 mo_base_types.f90 mo_buffer_fft.f90 mo_call_trans.f90 \
	mo_clim.f90 mo_cloud.f90 mo_column.f90 mo_constants.f90 \
	mo_control.f90 mo_convect_tables.f90 mo_couple.f90 \
	mo_cumulus_flux.f90 mo_decomposition.f90 mo_diag_amip2.f90 \
	mo_diag_dynamics.f90 mo_diag_radiation.f90 mo_diag_tendency.f90 \
	mo_diff.f90 mo_doctor.f90 mo_essl_dft.f90 mo_exception.f90 \
	mo_fft991.f90 mo_fft992.f90 mo_field.f90 mo_filename.f90 \
	mo_forecast_switches.f90 mo_gaussgrid.f90 mo_geoloc.f90 mo_gl1.f90 \
	mo_global_op.f90 mo_greenhouse_gases.f90 mo_gwspectrum.f90 \
	mo_hd_highres_io.f90 mo_hdiff.f90 mo_hyb.f90 mo_hydrology.f90 \
	mo_interpo.f90 mo_io.f90 mo_kind.f90 mo_legendre.f90 \
	mo_linked_list.f90 mo_lw_clop.f90 mo_machine.f90 mo_memory_base.f90 \
	mo_memory_f.f90 mo_memory_g1a.f90 mo_memory_g1b.f90 mo_memory_g2a.f90 \
	mo_memory_g2b.f90 mo_memory_g3a.f90 mo_memory_g3b.f90 \
	mo_memory_gl.f90 mo_memory_ls.f90 mo_memory_sp.f90 \
	mo_memory_streams.f90 mo_midatm.f90 mo_mkl_dft.f90 mo_mpi.f90 \
	mo_namelist.f90 mo_netcdf.f90 mo_nmi.f90 mo_nudging.f90 \
	mo_nudging_buffer.f90 mo_nudging_constants.f90 mo_nudging_init.f90 \
	mo_nudging_io.f90 mo_nudging_pattern.f90 mo_nudging_sst.f90 \
	mo_nudging_utils.f90 mo_o3_lwb.f90 mo_o3clim.f90 mo_orbit.f90 \
	mo_output.f90 mo_param_switches.f90 mo_parameters.f90 mo_parrrtm.f90 \
	mo_physc1.f90 mo_physc2.f90 mo_port_test.f90 mo_post.f90 \
	mo_profile.f90 mo_radiation.f90 mo_real_timer.f90 mo_rrta_140gp.f90 \
	mo_rrtab.f90 mo_rrtbg2.f90 mo_rrtftr.f90 mo_rrtrf.f90 mo_rrtwn.f90 \
	mo_scan_buffer.f90 mo_semi_impl.f90 mo_semi_lagrangian.f90 mo_so4.f90 \
	mo_solmon.f90 mo_spectral.f90 mo_spitfire.f90 mo_ssodrag.f90 \
	mo_ssortns.f90 mo_sst.f90 mo_sub_echam.f90 mo_sub_nml.f90 mo_sw.f90 \
	mo_sw_clop.f90 mo_test_trans.f90 mo_time_base.f90 mo_time_control.f90 \
	mo_time_conversion.f90 mo_time_event.f90 mo_time_manager.f90 \
	mo_timer.f90 mo_timestamp.f90 mo_tmp_buffer.f90 mo_tpcore.f90 \
	mo_tracdef.f90 mo_tracer.f90 mo_transpose.f90 mo_tropopause.f90 \
	mo_truncation.f90 mo_upper_sponge.f90 mo_util_string.f90 \
	mo_vegetation.f90 mo_vsop87.f90

OBJS := $(SRCS:.f90=.o)

all: $(OBJS)

clean:
	rm -f $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

m_alloc_mods.o: mo_control.o mo_diff.o mo_hdiff.o mo_hyb.o mo_physc2.o \
	mo_semi_impl.o
mo_aero_gads.o: mo_constants.o mo_kind.o mo_memory_g1a.o mo_memory_g3b.o \
	mo_radiation.o
mo_aero_tanre.o: mo_decomposition.o mo_gaussgrid.o mo_geoloc.o mo_kind.o \
	mo_transpose.o
mo_albedo.o: mo_constants.o mo_control.o mo_exception.o mo_kind.o
mo_base_types.o: mo_kind.o
mo_buffer_fft.o: mo_decomposition.o mo_kind.o
mo_call_trans.o: mo_buffer_fft.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_linked_list.o mo_memory_f.o mo_memory_g1a.o mo_memory_ls.o \
	mo_memory_sp.o mo_mpi.o mo_scan_buffer.o mo_test_trans.o mo_timer.o \
	mo_transpose.o
mo_clim.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o mo_io.o \
	mo_kind.o mo_mpi.o mo_transpose.o
mo_cloud.o: mo_constants.o mo_control.o mo_doctor.o mo_exception.o mo_kind.o \
	mo_mpi.o
mo_column.o: mo_advection.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_filename.o mo_gaussgrid.o mo_kind.o mo_memory_g1a.o \
	mo_memory_g2a.o mo_memory_gl.o mo_mpi.o mo_namelist.o \
	mo_scan_buffer.o mo_time_control.o mo_transpose.o \
	$(INCLUDE)/columnctl.inc
mo_constants.o: mo_kind.o
mo_control.o: mo_kind.o
mo_convect_tables.o: mo_constants.o mo_exception.o mo_kind.o
mo_couple.o: mo_constants.o mo_control.o mo_decomposition.o mo_exception.o \
	mo_gaussgrid.o mo_kind.o mo_memory_g3b.o mo_mpi.o mo_time_control.o \
	mo_time_conversion.o mo_transpose.o
mo_cumulus_flux.o: mo_control.o mo_doctor.o mo_exception.o mo_kind.o mo_mpi.o
mo_decomposition.o: mo_doctor.o mo_exception.o mo_mpi.o mo_util_string.o
mo_diag_amip2.o: mo_constants.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_filename.o mo_gaussgrid.o mo_geoloc.o mo_io.o mo_kind.o mo_mpi.o \
	mo_netcdf.o mo_time_control.o mo_transpose.o
mo_diag_dynamics.o: mo_constants.o mo_control.o mo_decomposition.o \
	mo_doctor.o mo_forecast_switches.o mo_gaussgrid.o mo_kind.o \
	mo_memory_g3a.o mo_memory_g3b.o mo_memory_gl.o mo_mpi.o \
	mo_scan_buffer.o mo_semi_impl.o mo_time_control.o mo_transpose.o
mo_diag_radiation.o: mo_control.o mo_doctor.o mo_kind.o mo_mpi.o
mo_diag_tendency.o: mo_control.o mo_essl_dft.o mo_exception.o mo_fft992.o \
	mo_kind.o mo_legendre.o mo_linked_list.o mo_memory_base.o \
	mo_memory_sp.o mo_mkl_dft.o mo_mpi.o mo_output.o mo_spectral.o \
	mo_time_control.o mo_truncation.o
mo_diff.o: mo_parameters.o
mo_essl_dft.o: mo_kind.o
mo_exception.o: mo_doctor.o mo_mpi.o
mo_fft991.o: mo_doctor.o mo_kind.o mo_mpi.o
mo_fft992.o: mo_exception.o mo_kind.o
mo_field.o: mo_kind.o
mo_filename.o: mo_control.o mo_doctor.o mo_exception.o mo_mpi.o \
	mo_time_control.o
mo_gaussgrid.o: mo_constants.o mo_control.o mo_kind.o
mo_geoloc.o: mo_decomposition.o mo_gaussgrid.o mo_kind.o mo_memory_base.o \
	mo_radiation.o mo_transpose.o
mo_gl1.o: mo_kind.o
mo_global_op.o: mo_decomposition.o mo_exception.o mo_kind.o mo_mpi.o \
	mo_transpose.o
mo_greenhouse_gases.o: mo_constants.o mo_exception.o mo_io.o mo_kind.o \
	mo_mpi.o mo_netcdf.o mo_radiation.o mo_time_control.o \
	mo_time_conversion.o
mo_gwspectrum.o: mo_kind.o
mo_hd_highres_io.o: mo_exception.o mo_filename.o mo_kind.o mo_time_control.o
mo_hdiff.o: mo_kind.o
mo_hyb.o: mo_constants.o mo_control.o mo_kind.o mo_parameters.o
mo_hydrology.o: mo_constants.o mo_control.o mo_decomposition.o mo_exception.o \
	mo_filename.o mo_gaussgrid.o mo_hd_highres_io.o mo_io.o mo_kind.o \
	mo_memory_g3b.o mo_mpi.o mo_time_control.o mo_transpose.o
mo_interpo.o: mo_kind.o
mo_io.o: mo_control.o mo_couple.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_filename.o mo_gaussgrid.o mo_kind.o \
	mo_linked_list.o mo_memory_base.o mo_mpi.o mo_netcdf.o \
	mo_parameters.o mo_time_control.o mo_transpose.o mo_util_string.o
mo_legendre.o: mo_constants.o mo_control.o mo_decomposition.o mo_gaussgrid.o \
	mo_kind.o
mo_linked_list.o: mo_doctor.o mo_exception.o mo_filename.o mo_kind.o \
	mo_netcdf.o mo_util_string.o
mo_lw_clop.o: mo_kind.o mo_parrrtm.o
mo_machine.o: mo_doctor.o mo_exception.o mo_kind.o mo_mpi.o
mo_memory_base.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o \
	mo_kind.o mo_linked_list.o mo_netcdf.o mo_time_control.o \
	mo_time_event.o
mo_memory_f.o: mo_kind.o mo_linked_list.o mo_memory_base.o mo_netcdf.o
mo_memory_g1a.o: mo_kind.o mo_linked_list.o mo_memory_base.o mo_memory_gl.o \
	mo_netcdf.o mo_tracdef.o
mo_memory_g1b.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_g2a.o: mo_kind.o mo_linked_list.o mo_memory_base.o mo_netcdf.o
mo_memory_g2b.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_g3a.o: mo_kind.o mo_memory_g3b.o
mo_memory_g3b.o: mo_control.o mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_gl.o: mo_filename.o mo_kind.o mo_linked_list.o mo_memory_base.o \
	mo_netcdf.o mo_tracdef.o
mo_memory_ls.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_sp.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_streams.o: mo_buffer_fft.o mo_control.o mo_decomposition.o \
	mo_doctor.o mo_filename.o mo_memory_base.o mo_memory_f.o \
	mo_memory_g1a.o mo_memory_g1b.o mo_memory_g2a.o mo_memory_g2b.o \
	mo_memory_g3a.o mo_memory_g3b.o mo_memory_gl.o mo_memory_ls.o \
	mo_memory_sp.o mo_mpi.o mo_port_test.o mo_sub_nml.o mo_tracer.o
mo_midatm.o: mo_constants.o mo_control.o mo_doctor.o mo_exception.o \
	mo_gaussgrid.o mo_geoloc.o mo_gwspectrum.o mo_kind.o mo_memory_g1a.o \
	mo_memory_g2a.o
mo_mkl_dft.o: mo_kind.o
mo_mpi.o: mo_doctor.o mo_kind.o mo_util_string.o
mo_namelist.o: mo_exception.o mo_util_string.o
mo_netcdf.o: mo_control.o mo_doctor.o mo_exception.o mo_filename.o mo_kind.o \
	mo_parameters.o
mo_nmi.o: mo_constants.o mo_control.o mo_decomposition.o mo_exception.o \
	mo_hyb.o mo_kind.o mo_linked_list.o mo_memory_base.o mo_memory_sp.o \
	mo_mpi.o mo_namelist.o mo_nudging_buffer.o mo_output.o mo_semi_impl.o \
	mo_spectral.o mo_time_control.o mo_time_conversion.o mo_truncation.o
mo_nudging.o: mo_control.o mo_decomposition.o mo_exception.o mo_kind.o \
	mo_memory_sp.o mo_nmi.o mo_nudging_buffer.o mo_nudging_constants.o \
	mo_nudging_io.o mo_nudging_pattern.o mo_time_control.o
mo_nudging_buffer.o: mo_control.o mo_decomposition.o mo_exception.o mo_kind.o \
	mo_linked_list.o mo_memory_base.o mo_mpi.o mo_nudging_constants.o \
	mo_output.o mo_time_control.o
mo_nudging_constants.o: mo_kind.o mo_parameters.o
mo_nudging_init.o: mo_control.o mo_decomposition.o mo_exception.o \
	mo_filename.o mo_kind.o mo_memory_base.o mo_mpi.o mo_namelist.o \
	mo_nudging_buffer.o mo_nudging_constants.o mo_nudging_io.o \
	mo_nudging_sst.o mo_nudging_utils.o mo_physc2.o mo_time_control.o \
	mo_time_conversion.o mo_transpose.o mo_truncation.o
mo_nudging_io.o: mo_control.o mo_decomposition.o mo_exception.o mo_filename.o \
	mo_kind.o mo_mpi.o mo_nmi.o mo_nudging_buffer.o \
	mo_nudging_constants.o mo_nudging_pattern.o mo_nudging_utils.o \
	mo_time_control.o mo_transpose.o
mo_nudging_pattern.o: mo_control.o mo_exception.o mo_kind.o mo_memory_sp.o \
	mo_nudging_buffer.o mo_nudging_constants.o
mo_nudging_sst.o: mo_decomposition.o mo_exception.o mo_filename.o mo_kind.o \
	mo_memory_g3b.o mo_mpi.o mo_nudging_constants.o mo_nudging_utils.o \
	mo_physc2.o mo_time_control.o mo_transpose.o
mo_nudging_utils.o: mo_control.o mo_diag_tendency.o mo_exception.o mo_kind.o \
	mo_mpi.o mo_nudging_buffer.o mo_nudging_constants.o mo_post.o \
	mo_spectral.o mo_time_control.o mo_time_conversion.o
mo_o3_lwb.o: mo_decomposition.o mo_gaussgrid.o mo_geoloc.o mo_kind.o \
	mo_transpose.o
mo_o3clim.o: mo_constants.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_filename.o mo_interpo.o mo_io.o mo_kind.o mo_mpi.o \
	mo_netcdf.o mo_radiation.o mo_time_control.o mo_transpose.o
mo_orbit.o: mo_constants.o mo_exception.o mo_kind.o mo_radiation.o \
	mo_vsop87.o
mo_output.o: mo_constants.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_filename.o mo_gaussgrid.o mo_kind.o \
	mo_linked_list.o mo_memory_base.o mo_mpi.o mo_netcdf.o \
	mo_time_control.o mo_time_conversion.o mo_tracdef.o mo_transpose.o \
	mo_util_string.o $(INCLUDE)/cdi.inc
mo_physc1.o: mo_kind.o
mo_physc2.o: mo_kind.o mo_parameters.o
mo_port_test.o: mo_linked_list.o mo_memory_base.o mo_time_event.o
mo_radiation.o: mo_kind.o mo_time_control.o
mo_real_timer.o: mo_doctor.o mo_exception.o mo_kind.o mo_mpi.o \
	mo_util_string.o
mo_rrta_140gp.o: mo_kind.o
mo_rrtab.o: mo_kind.o
mo_rrtbg2.o: mo_kind.o
mo_rrtftr.o: mo_kind.o mo_parrrtm.o
mo_rrtrf.o: mo_kind.o
mo_rrtwn.o: mo_kind.o
mo_scan_buffer.o: mo_decomposition.o mo_kind.o mo_tracer.o
mo_semi_impl.o: mo_kind.o
mo_semi_lagrangian.o: mo_advection.o mo_constants.o mo_control.o \
	mo_decomposition.o mo_doctor.o mo_exception.o mo_gaussgrid.o \
	mo_global_op.o mo_hyb.o mo_kind.o mo_memory_base.o mo_memory_g1a.o \
	mo_memory_gl.o mo_mpi.o mo_parameters.o mo_scan_buffer.o mo_tracer.o \
	mo_transpose.o
mo_so4.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o mo_io.o \
	mo_kind.o mo_mpi.o mo_netcdf.o mo_time_control.o mo_transpose.o
mo_solmon.o: mo_control.o mo_doctor.o mo_exception.o mo_filename.o mo_kind.o \
	mo_mpi.o mo_radiation.o mo_time_control.o mo_time_conversion.o
mo_spectral.o: mo_exception.o mo_kind.o
mo_spitfire.o: mo_advection.o mo_constants.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_gaussgrid.o mo_hyb.o mo_kind.o mo_memory_g1a.o \
	mo_mpi.o mo_scan_buffer.o mo_time_control.o mo_tracer.o \
	mo_transpose.o
mo_ssodrag.o: mo_control.o mo_exception.o mo_kind.o
mo_ssortns.o: mo_constants.o mo_geoloc.o mo_kind.o mo_ssodrag.o \
	mo_time_control.o
mo_sst.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o mo_io.o \
	mo_kind.o mo_mpi.o mo_netcdf.o mo_time_control.o mo_transpose.o
mo_sub_echam.o: mo_kind.o mo_time_control.o mo_tracer.o
mo_sub_nml.o: mo_exception.o mo_kind.o mo_memory_base.o mo_mpi.o \
	mo_namelist.o mo_time_control.o mo_time_event.o mo_tracdef.o \
	mo_tracer.o
mo_sw.o: mo_kind.o
mo_sw_clop.o: mo_kind.o
mo_test_trans.o: mo_decomposition.o mo_doctor.o mo_exception.o mo_kind.o \
	mo_linked_list.o mo_memory_base.o mo_mpi.o mo_transpose.o
mo_time_base.o: mo_exception.o mo_kind.o
mo_time_control.o: mo_constants.o mo_control.o mo_exception.o mo_interpo.o \
	mo_kind.o mo_machine.o mo_mpi.o mo_param_switches.o mo_time_base.o \
	mo_time_conversion.o mo_time_event.o mo_time_manager.o
mo_time_conversion.o: mo_exception.o mo_kind.o mo_time_base.o
mo_time_event.o: mo_exception.o mo_kind.o mo_time_base.o mo_time_conversion.o
mo_time_manager.o: mo_exception.o mo_kind.o mo_time_base.o \
	mo_time_conversion.o
mo_timer.o: mo_real_timer.o
mo_tmp_buffer.o: mo_kind.o
mo_tpcore.o: mo_advection.o mo_constants.o mo_control.o mo_decomposition.o \
	mo_exception.o mo_gaussgrid.o mo_kind.o mo_memory_g1a.o mo_mpi.o \
	mo_scan_buffer.o mo_time_control.o mo_tracer.o mo_transpose.o
mo_tracdef.o: mo_kind.o mo_linked_list.o mo_parameters.o mo_time_conversion.o
mo_tracer.o: mo_advection.o mo_control.o mo_doctor.o mo_exception.o mo_kind.o \
	mo_linked_list.o mo_memory_base.o mo_memory_g1a.o mo_memory_gl.o \
	mo_mpi.o mo_parameters.o mo_semi_impl.o mo_time_conversion.o \
	mo_tracdef.o mo_util_string.o
mo_transpose.o: mo_buffer_fft.o mo_decomposition.o mo_exception.o mo_kind.o \
	mo_mpi.o
mo_tropopause.o: mo_cloud.o mo_constants.o mo_kind.o mo_time_control.o
mo_truncation.o: mo_control.o mo_doctor.o mo_kind.o mo_mpi.o mo_parameters.o
mo_upper_sponge.o: mo_control.o mo_decomposition.o mo_kind.o mo_memory_sp.o \
	mo_time_control.o
mo_vegetation.o: mo_kind.o
mo_vsop87.o: mo_kind.o

ifeq ($(strip $(ARCH)), SX)
mo_spitfire.o: mo_spitfire.f90
	$(F90) $(F90FLAGS) -pi auto exp=minmod exp=medan exp=putyslice exp=cfdot1dp2 exp=cfint1x2 line=2000 -c mo_spitfire.f90
mo_tpcore.o: mo_tpcore.f90
#warning: Don't change line to 2000. This will give wrong code on the SX!!!
	$(F90) $(F90FLAGS) -pi auto exp=xmist,fxppm,kmppm,lmppm,xtp noexp=map1_ppm_gp,ppm2m,steepz nest=3 line=1000 -c mo_tpcore.f90
mo_transpose.o: mo_transpose.f90
	$(F90) $(F90FLAGS) -Npi -c mo_transpose.f90
endif


ifeq ($(strip $(ARCH)), ES)
mo_spitfire.o: mo_spitfire.f90
	$(F90) $(F90FLAGS) -pi auto exp=minmod exp=medan exp=putyslice exp=cfdot1dp2 exp=cfint1x2 line=2000 -c mo_spitfire.f90
mo_tpcore.o: mo_tpcore.f90
#warning: Don't change line to 2000. This will give wrong code on the SX!!!
	$(F90) $(F90FLAGS) -pi auto exp=xmist,fxppm,kmppm,lmppm,xtp noexp=map1_ppm_gp,ppm2m,steepz nest=3 line=1000 -c mo_tpcore.f90
mo_transpose.o: mo_transpose.f90
	$(F90) $(F90FLAGS) -Npi -c mo_transpose.f90
endif

ifeq ($(strip $(ARCH)), CRAY_PVP)
mo_buffer_fft.o: mo_buffer_fft.f90
	$(F90) $(F90FLAGS) -Ovector1 -c mo_buffer_fft.f90
mo_grib.o: mo_grib.f90
	$(F90) $(F90FLAGS) -Ovector1 -c mo_grib.f90
endif

ifeq ($(strip $(ARCH)), CRAY_X1)
	$(F90) $(F90FLAGS) -Ovector1 -Oscalar1 -Ostream1 -c mo_io.f90
endif

ifeq ($(strip $(ARCH)), CRAY_XT3)
mo_transpose.o: mo_transpose.f90
	( F90FLAGS=" " ; $(F90) -O2  -c mo_transpose.f90 )
endif

ifeq ($(strip $(ARCH)), CRAY_XD1)
mo_transpose.o: mo_transpose.f90
	( F90FLAGS=" " ; $(F90) -O2  -c mo_transpose.f90 )
endif
