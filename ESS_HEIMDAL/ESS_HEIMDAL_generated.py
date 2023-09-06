#!/usr/bin/env python3
# Automatically generated file. 
# Format:    Python script code
# McStas <http://www.mcstas.org>
# Instrument: ESS_HEIMDAL.instr (ESS_HEIMDAL)
# Date:       Wed Sep  6 07:19:16 2023
# File:       ESS_HEIMDAL_generated.py

import mcstasscript as ms

# Python McStas instrument description
instr = ms.McStas_instr("ESS_HEIMDAL_generated", author = "McCode Py-Generator", origin = "ESS DMSC")

# Add collected DEPENDENCY strings
instr.set_dependency('')

# *****************************************************************************
# * Start of instrument 'ESS_HEIMDAL' generated code
# *****************************************************************************
# MCSTAS system dir is "/Users/pkwi/McStas/mcstas/3.x-dev/"


# *****************************************************************************
# * instrument 'ESS_HEIMDAL' and components DECLARE
# *****************************************************************************

# Instrument parameters:

HighRes = instr.add_parameter('double', 'HighRes', value=0, comment='Parameter type (double) added by McCode py-generator')
T = instr.add_parameter('double', 'T', value=25, comment='Parameter type (double) added by McCode py-generator')
tau = instr.add_parameter('double', 'tau', value=0, comment='Parameter type (double) added by McCode py-generator')
dmin = instr.add_parameter('double', 'dmin', value=0.5, comment='Parameter type (double) added by McCode py-generator')
dmax = instr.add_parameter('double', 'dmax', value=5.5, comment='Parameter type (double) added by McCode py-generator')

component_definition_metadata = {
}
instr.append_declare(r'''
  double tof,tofmin,tofmax;
  double det_rad,det_Linst,det_t0,det_th1,det_th2;
  double hm = 2*PI*K2V; // h/m_n
  double samplechoice = 0;
  double Tc = 525;
  double lam0;
  double dlam;
  double dPulse;
  double lmin = 0.65;
  double lmax = 2.35;
  double use_CoFe2O4=1, use_CoFe=1, use_magn=1, f1=0.5, f2=1, f3=1;
  double dfactor = 1;
  int dbins;
  
  char monopts[256];
  char monopts2[256];
  
  struct pieceswisepars {
    double min;
    double max;
    double* uppers;
    double* slopes;
    double* consts;
  };
  double Itime_CoFe2O4; // reaction time dependent intensity factor for CoFe2O4
  double Itime_CoFe; // reaction time dependent intensity factor for CoFe
  struct pieceswisepars pars_CoFe2O4;
  struct pieceswisepars pars_CoFe;
  int piecewiselin(double x, double* val, struct pieceswisepars* fct, int sz);

  //
  // determine Itime_CoFe2O4 and Itime_CoFe
  int piecewiselin(double x, double* val, struct pieceswisepars* fct, int sz) {
    // NOTE: the user is responsible
    double* u = fct->uppers;
    double* s = fct->slopes;
    double* c = fct->consts;
    double min = fct->min;

    // use first found segment
    int i = 0;
    double low, high;
    for (i=0;i<sz;i++) {
      // determine segment lower limit
      if (i==0) low = min;
      else low = u[i-1];

      // determine segment upper limit
      high = u[i];
      if (low<=x && x<=high) {
        fprintf(stdout, "\nl:%1f h:%1f s:%1f c:%1f x:%1f", low, high, s[i], c[i], x);
        *val = s[i]*x + c[i];
        return s[i]*x + c[i];
        //return 0;
      }
    }
    return 1;
  }

''')


instr.append_initialize(r'''


// angular range
det_th1 = 10;
det_th2 = 170;
det_rad = 1.5;
det_Linst = 150+det_rad;
tof = det_Linst/hm*lam0*1e6; // tof in [us]
tofmin = det_Linst/hm*lmin*1e6; // tof in [us]
tofmax = det_Linst/hm*lmax*1e6; // tof in [us]
lam0 = (lmax+lmin)/2.0;
dlam = (lmax-lmin)/2.0;

dbins = ceil((dmax-dmin)/0.002);
 
//
// determine temperature-dependent d-spacing stretch
dfactor = 0.02*T/Tc + 1; // latices are scaled by 1% at Tc...



// NOTE:  To enable simple particle weight multiplication, the chosen
//        piecewiselin functions have a max value of 1 or less.

double min = 0;
double u[2] = {50, 600};
double s[2] = {0, 1/250.0};
double c[2] = {0, -50*1/250.0};
pars_CoFe.min = min;
pars_CoFe.uppers = u;
pars_CoFe.slopes = s;
pars_CoFe.consts = c;

piecewiselin(tau, &Itime_CoFe, &pars_CoFe, 2);

double min2 = 0;
double u2[1] = {300};
double s2[1] = {-0.4/300};
double c2[1] = {0.4};
pars_CoFe2O4.min = min2;
pars_CoFe2O4.uppers = u2;
pars_CoFe2O4.slopes = s2;
pars_CoFe2O4.consts = c2;

piecewiselin(tau, &Itime_CoFe2O4, &pars_CoFe2O4, 1);

// DEBUG OUTPUT:
//printf("\nCHECK: CoFe=%1f CoFe2O4=%1f\n", Itime_CoFe, Itime_CoFe2O4);
//exit(1);

 if (HighRes==0) {
   dPulse=3.86e-3*0.05;
 } else {
   dPulse=3.86e-3*0.0017;
 }

 sprintf(monopts,"banana theta limits=[%g %g] bins=600, tof limits=[%g %g] bins=800",det_th1, det_th2, 0.95*tofmin/1e6, tofmax*1.05/1e6);
 sprintf(monopts2,"banana theta limits=[%g %g] bins=600, tof limits=[%g %g] bins=800",-170.0, -55.0, 0.95*tofmin/1e6, tofmax*1.05/1e6);
''')


# *****************************************************************************
# * instrument 'ESS_HEIMDAL' TRACE
# *****************************************************************************

# Comp instance origin, placement and parameters
origin = instr.add_component('origin','Progress_bar')

origin.profile = '"NULL"'
origin.percent = '10'
origin.flag_save = '0'
origin.minutes = '0'

# Comp instance source, placement and parameters
source = instr.add_component('source','Source_simple', AT=['0', '0', '0'], AT_RELATIVE='origin', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='origin')
# EXTEND at source
source.append_EXTEND(r'''
  samplechoice=rand01();
  t=rand01()*dPulse;
  p*=dPulse/3.86e-3;
''')



source.radius = '0.1'
source.yheight = '0'
source.xwidth = '0'
source.dist = '0'
source.focus_xw = '0.01'
source.focus_yh = '0.07'
source.E0 = '0'
source.dE = '0'
source.lambda0 = 'lam0'
source.dlambda = 'dlam'
source.flux = '1e14'
source.gauss = '0'
source.target_index = '3'

# Comp instance SampleArm, placement and parameters
SampleArm = instr.add_component('SampleArm','Arm', AT=['0', '0', '150'], AT_RELATIVE='source', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='source')
# SPLIT 32 times at SampleArm
SampleArm.set_SPLIT('32')


# Comp instance sample_2, placement and parameters
sample_2 = instr.add_component('sample_2','PowderN_dspace_factor', AT=['0', '0', '150'], AT_RELATIVE='source', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='source')
# WHEN ( samplechoice >= 0.43 && samplechoice < 0.86 && use_CoFe == 1 ) at sample_2
sample_2.set_WHEN('( samplechoice >= 0.43 && samplechoice < 0.86 && use_CoFe == 1 )')
# EXTEND at sample_2
sample_2.append_EXTEND(r'''
    // peak intensity depends on reaction time
    p *= Itime_CoFe;
    p *= f1;
''')



sample_2.reflections = '"CoFe_56273.txt"'
sample_2.geometry = '"NULL"'
sample_2.format = '{ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }'
sample_2.radius = '0.005'
sample_2.yheight = '0.07'
sample_2.xwidth = '0'
sample_2.zdepth = '0'
sample_2.thickness = '0'
sample_2.pack = '1'
sample_2.Vc = '0'
sample_2.sigma_abs = '0'
sample_2.sigma_inc = '0'
sample_2.delta_d_d = '0'
sample_2.p_inc = '0.1'
sample_2.p_transmit = '0'
sample_2.DW = '0'
sample_2.nb_atoms = '1'
sample_2.d_omega = '0'
sample_2.d_phi = '10'
sample_2.tth_sign = '0'
sample_2.p_interact = '0.8'
sample_2.dspace_factor = 'dfactor'
sample_2.concentric = '0'
sample_2.density = '0'
sample_2.weight = '0'
sample_2.barns = '1'
sample_2.Strain = '0'
sample_2.focus_flip = '0'
sample_2.target_index = '0'

# Comp instance sample_1, placement and parameters
sample_1 = instr.add_component('sample_1','PowderN_dspace_factor', AT=['0', '0', '150'], AT_RELATIVE='source', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='source')
# WHEN ( samplechoice < 0.43 && use_CoFe2O4 == 1 ) at sample_1
sample_1.set_WHEN('( samplechoice < 0.43 && use_CoFe2O4 == 1 )')
# EXTEND at sample_1
sample_1.append_EXTEND(r'''
    // peak intensity varies as a function of reaction time
    p *= Itime_CoFe2O4;
    p *= f2;
''')



sample_1.reflections = '"CoFe2O4_109044.txt"'
sample_1.geometry = '"NULL"'
sample_1.format = '{ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }'
sample_1.radius = '0.005'
sample_1.yheight = '0.07'
sample_1.xwidth = '0'
sample_1.zdepth = '0'
sample_1.thickness = '0'
sample_1.pack = '1'
sample_1.Vc = '0'
sample_1.sigma_abs = '0'
sample_1.sigma_inc = '0'
sample_1.delta_d_d = '0'
sample_1.p_inc = '0.1'
sample_1.p_transmit = '0'
sample_1.DW = '0'
sample_1.nb_atoms = '1'
sample_1.d_omega = '0'
sample_1.d_phi = '10'
sample_1.tth_sign = '0'
sample_1.p_interact = '0.8'
sample_1.dspace_factor = 'dfactor'
sample_1.concentric = '0'
sample_1.density = '0'
sample_1.weight = '0'
sample_1.barns = '1'
sample_1.Strain = '0'
sample_1.focus_flip = '0'
sample_1.target_index = '0'

# Comp instance Magnet, placement and parameters
Magnet = instr.add_component('Magnet','PowderN_dspace_factor', AT=['0', '0', '150'], AT_RELATIVE='source', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='source')
# WHEN ( samplechoice >= 0.86 ) at Magnet
Magnet.set_WHEN('( samplechoice >= 0.86 )')
# EXTEND at Magnet
Magnet.append_EXTEND(r'''
    // peak intensity diminishes as a function of temperature
    if (INSTRUMENT_GETPAR(T)<Tc && INSTRUMENT_GETPAR(T)>=0) {
      p *= (Tc-INSTRUMENT_GETPAR(T))/(Tc);
      p *= f3;
      p *= Itime_CoFe2O4;
    }
    else
        ABSORB;
''')



Magnet.reflections = '"magnetic_peak.txt"'
Magnet.geometry = '"NULL"'
Magnet.format = '{ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }'
Magnet.radius = '0.001'
Magnet.yheight = '0.07'
Magnet.xwidth = '0'
Magnet.zdepth = '0'
Magnet.thickness = '0'
Magnet.pack = '1'
Magnet.Vc = '0'
Magnet.sigma_abs = '0'
Magnet.sigma_inc = '0'
Magnet.delta_d_d = '0'
Magnet.p_inc = '0'
Magnet.p_transmit = '0'
Magnet.DW = '0'
Magnet.nb_atoms = '1'
Magnet.d_omega = '0'
Magnet.d_phi = '10'
Magnet.tth_sign = '0'
Magnet.p_interact = '0.8'
Magnet.dspace_factor = 'dfactor'
Magnet.concentric = '0'
Magnet.density = '0'
Magnet.weight = '0'
Magnet.barns = '1'
Magnet.Strain = '0'
Magnet.focus_flip = '0'
Magnet.target_index = '0'

# Comp instance psdtof, placement and parameters
psdtof = instr.add_component('psdtof','NPI_tof_theta_monitor', AT=['0', '0', '0'], AT_RELATIVE='sample_1', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='sample_1')

psdtof.filename = '"tof_theta.dat"'
psdtof.radius = '1.5'
psdtof.yheight = '1'
psdtof.tmin = '0.95 * tofmin'
psdtof.tmax = 'tofmax * 1.05'
psdtof.amin = 'det_th1'
psdtof.amax = 'det_th2'
psdtof.restore_neutron = '1'
psdtof.verbose = '0'
psdtof.nt = '800'
psdtof.na = '600'

# Comp instance npi_tof_dhkl_detector, placement and parameters
npi_tof_dhkl_detector = instr.add_component('npi_tof_dhkl_detector','NPI_tof_dhkl_detector', AT=['0', '0', '0'], AT_RELATIVE='sample_1', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='sample_1')

npi_tof_dhkl_detector.filename = '"dhkl.dat"'
npi_tof_dhkl_detector.radius = 'det_rad'
npi_tof_dhkl_detector.yheight = '1.0'
npi_tof_dhkl_detector.zdepth = '0.01'
npi_tof_dhkl_detector.amin = 'det_th1'
npi_tof_dhkl_detector.amax = 'det_th2'
npi_tof_dhkl_detector.nd = 'dbins'
npi_tof_dhkl_detector.na = '800'
npi_tof_dhkl_detector.nt = '800'
npi_tof_dhkl_detector.nev = '1'
npi_tof_dhkl_detector.d_min = 'dmin'
npi_tof_dhkl_detector.d_max = 'dmax'
npi_tof_dhkl_detector.time0 = '0'
npi_tof_dhkl_detector.Linst = '152'
npi_tof_dhkl_detector.Lc = '0'
npi_tof_dhkl_detector.res_x = '0.002'
npi_tof_dhkl_detector.res_y = '0.005'
npi_tof_dhkl_detector.res_t = '1e-6'
npi_tof_dhkl_detector.mu = '1.0'
npi_tof_dhkl_detector.mod_shift = '0'
npi_tof_dhkl_detector.modulation = '0'
npi_tof_dhkl_detector.mod_dt = '0'
npi_tof_dhkl_detector.mod_twidth = 'dPulse'
npi_tof_dhkl_detector.mod_d0_table = '"dhkl.tab"'
npi_tof_dhkl_detector.verbose = '0'
npi_tof_dhkl_detector.restore_neutron = '1'

# Instruct McStasscript not to 'check everythng'
instr.settings(checks=False)
# (also, this is where one can add e.g. seed=1000, ncount=1e7, mpi=8, openacc=True, force_compile=False etc.)


# Visualise with default parameters.
instr.show_instrument()


# Generate a dataset with default parameters.
data = instr.backengine()

# Overview plot:
ms.make_sub_plot(data)


# Other useful commands follow...

# One plot pr. window
#ms.make_plot(data)

# Load another dataset
#data2 = ms.load_data('some_other_folder')

# Adjusting a specific plot
#ms.name_plot_options("PSD_4PI", data, log=1, colormap="hot", orders_of_mag=5)


# Bring up the 'interface' - only relevant in Jupyter
#%matplotlib widget
#import mcstasscript.jb_interface as ms_widget
#ms_widget.show(data)


# Bring up the simulation 'interface' - only relevant in Jupyter
#%matplotlib widget
#import mcstasscript.jb_interface as ms_widget
#sim_widget = ms_widget.SimInterface(instr)
#sim_widget.show_interface()
#data = sim_widget.get_data()


# end of generated Python code ESS_HEIMDAL_generated.py 
