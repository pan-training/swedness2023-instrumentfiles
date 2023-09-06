#!/usr/bin/env python3
# Automatically generated file. 
# Format:    Python script code
# McStas <http://www.mcstas.org>
# Instrument: Nanodisc_Loki.instr (Loki)
# Date:       Wed Sep  6 07:19:52 2023
# File:       Nanodisc_Loki_generated.py

import mcstasscript as ms

# Python McStas instrument description
instr = ms.McStas_instr("Loki_generated", author = "McCode Py-Generator", origin = "ESS DMSC")

# Add collected DEPENDENCY strings
instr.set_dependency(' -Wl,-rpath,CMD(mcpl-config --show libdir) -LCMD(mcpl-config --show libdir) -lmcpl -ICMD(mcpl-config --show includedir)')

# *****************************************************************************
# * Start of instrument 'Loki' generated code
# *****************************************************************************
# MCSTAS system dir is "/Users/pkwi/McStas/mcstas/3.x-dev/"


# *****************************************************************************
# * instrument 'Loki' and components DECLARE
# *****************************************************************************

# Instrument parameters:

Instrument_Mode = instr.add_parameter('int', 'Instrument_Mode', value=0, comment='Parameter type (int) added by McCode py-generator')
Detector_Distance = instr.add_parameter('double', 'Detector_Distance', value=2.0, comment='Parameter type (double) added by McCode py-generator')
SAMPLE = instr.add_parameter('int', 'SAMPLE', value=0, comment='Parameter type (int) added by McCode py-generator')
Solvent_Fraction = instr.add_parameter('double', 'Solvent_Fraction', value=1.0, comment='Parameter type (double) added by McCode py-generator')
Mix_Fraction = instr.add_parameter('double', 'Mix_Fraction', value=0.0, comment='Parameter type (double) added by McCode py-generator')
Incoherent_Solvent = instr.add_parameter('int', 'Incoherent_Solvent', value=0, comment='Parameter type (int) added by McCode py-generator')

component_definition_metadata = {
}
instr.append_declare(r'''



double DetectorRadius;
double DetectorSize;
double Mix_Chance;

  //Strings
  char virtualoutput_filename[256];
  char calibration_filename[256];
  char mcploutput_filename[256];
  char Sample_String[256];

  // For sample controle:


int Sam;
int SAMPLE_HOLD;
int SAMPLE_HOLD2;
int SAMPLE_Select;
int Repeat_Number;
int Plot_Frac;

//int VI;
int Beamstop;
double Rmin;
double T_min;
double T_max;
double P_Interact;

double sigma_Head_Deutorated;
double sigma_CH2_Deutorated;
double sigma_CH3_Deutorated;

double sigma_Head_Hydrated;
double sigma_CH2_Hydrated;
double sigma_CH3_Hydrated;

double sigma_Head_Mixed;
double sigma_CH2_Mixed;
double sigma_CH3_Mixed;

double RhoH2O;
double RhoD2O;
double RhoSolv;

double AxisRatio;
double N_Lipids;
double AreaPrHead;
double H_1MSP;
double V_1MSP;
double V_1MSPE3;
double V_1Headgruope;
double V_1CH2Tail;
double V_1CH3Tail;
double SLD_1MSP;
double SLD_1MSPE3;
double Rough;
double Conc;
double AbsCrossection;

double sigma_abs_H2O;
double sigma_inc_H2O;
double sigma_abs_D2O;
double sigma_inc_D2O;

double sigma_abs_solvent;
double sigma_inc_solvent;
double RhoSolv_Inc;
''')


instr.append_initialize(r'''

//All of these could be input perameters
Repeat_Number=5;
DetectorSize = 1.0;
Beamstop = 0;
P_Interact = 0.0;
Sam = 999;
Rmin = sqrt(0.07*0.07+0.07*0.07);
  // Calculation of sample statistics:

if(Instrument_Mode < 0){Instrument_Mode = 0;}
if(Instrument_Mode > 5){Instrument_Mode = 0;}
if(Instrument_Mode > 2){Repeat_Number = 2;}

if(Detector_Distance == 2.0){
    if(Instrument_Mode == 0){
        sprintf(virtualoutput_filename,"inputs/C3_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal0_2.dat");}
    if(Instrument_Mode == 1){
        sprintf(virtualoutput_filename,"inputs/C5_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal1_2.dat");}
    if(Instrument_Mode == 2){
        sprintf(virtualoutput_filename,"inputs/C8_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal2_2.dat");}
    if(Instrument_Mode == 3){
        sprintf(virtualoutput_filename,"inputs/C3_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal3_2.dat");}
    if(Instrument_Mode == 4){
        sprintf(virtualoutput_filename,"inputs/C5_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal4_2.dat");}
    if(Instrument_Mode == 5){
        sprintf(virtualoutput_filename,"inputs/C8_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal5_2.dat");}
}

else if(Detector_Distance == 10.0){
    if(Instrument_Mode == 0){
        sprintf(virtualoutput_filename,"inputs/C3_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal0_10.dat");}
    if(Instrument_Mode == 1){
        sprintf(virtualoutput_filename,"inputs/C5_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal1_10.dat");}
    if(Instrument_Mode == 2){
        sprintf(virtualoutput_filename,"inputs/C8_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal2_10.dat");}
    if(Instrument_Mode == 3){
        sprintf(virtualoutput_filename,"inputs/C3_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal3_10.dat");}
    if(Instrument_Mode == 4){
        sprintf(virtualoutput_filename,"inputs/C5_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal4_10.dat");}
    if(Instrument_Mode == 5){
        sprintf(virtualoutput_filename,"inputs/C8_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal5_10.dat");}
}
else{
    if(Instrument_Mode == 0){
        sprintf(virtualoutput_filename,"inputs/C3_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal0_5.dat");}
    if(Instrument_Mode == 1){
        sprintf(virtualoutput_filename,"inputs/C5_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal1_5.dat");}
    if(Instrument_Mode == 2){
        sprintf(virtualoutput_filename,"inputs/C8_L115.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal2_5.dat");}
    if(Instrument_Mode == 3){
        sprintf(virtualoutput_filename,"inputs/C3_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal3_5.dat");}
    if(Instrument_Mode == 4){
        sprintf(virtualoutput_filename,"inputs/C5_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal4_5.dat");}
    if(Instrument_Mode == 5){
        sprintf(virtualoutput_filename,"inputs/C8_L214.mcpl.gz");
        sprintf(calibration_filename,"calibrations/cal5_5.dat");}
}

T_min=22000.0;
T_max=170000.0;
DetectorRadius = sqrt((DetectorSize*DetectorSize) + (DetectorSize*DetectorSize));


SAMPLE_HOLD = SAMPLE;
SAMPLE_HOLD2 = SAMPLE;

AxisRatio=1.2;
AreaPrHead=60.13;
H_1MSP=24.0;
N_Lipids=120.0;
V_1MSP=27146.4;
SLD_1MSP=9.23E-10;
V_1Headgruope=319.0;
V_1CH2Tail=754.2;
V_1CH3Tail=108.6;
Rough=2.0;
Conc=50.0;
AbsCrossection=0.0;

sigma_Head_Deutorated = 247.55E-13*(1-0.5 * Mix_Fraction) + 60.06E-13 * 0.5 * Mix_Fraction;
sigma_CH2_Deutorated = 479.91E-13*(1-0.5 * Mix_Fraction) -20.05E-13 * 0.5 * Mix_Fraction;
sigma_CH3_Deutorated = 53.34E-13*(1-0.5 * Mix_Fraction) -9.16E-13 * 0.5 * Mix_Fraction;

sigma_Head_Hydrated = 60.06E-13*(1-0.5 * Mix_Fraction) + 247.55E-13 * 0.5 * Mix_Fraction;
sigma_CH2_Hydrated = -20.05E-13*(1-0.5 * Mix_Fraction) + 479.91E-13 * 0.5 * Mix_Fraction;
sigma_CH3_Hydrated = -9.16E-13*(1-0.5 * Mix_Fraction) + 53.34E-13 * 0.5 * Mix_Fraction;

sigma_Head_Mixed = 153.80E-13;
sigma_CH2_Mixed = 229.93E-13;
sigma_CH3_Mixed = 22.09E-13;

RhoH2O = -1.678E-13/30.;
RhoD2O = 19.15E-13/30.;
RhoSolv = RhoH2O*(1.0-Solvent_Fraction) + RhoD2O*Solvent_Fraction;

sigma_abs_H2O = 2*0.3326+0.00019;
sigma_inc_H2O = 2*80.27+0.0;
sigma_abs_D2O = 2*0.000519+0.00019;
sigma_inc_D2O = 2*2.05+0.0;

sigma_abs_solvent = sigma_abs_H2O*(1.0-Solvent_Fraction) + sigma_abs_D2O*Solvent_Fraction;
sigma_inc_solvent = sigma_inc_H2O*(1.0-Solvent_Fraction) + sigma_inc_D2O*Solvent_Fraction;


RhoSolv_Inc = sigma_inc_solvent/30.;


if(Incoherent_Solvent == 1){
  Conc = 0.0;
}
if(Incoherent_Solvent == 2){
  RhoSolv_Inc = 0.0;
}

''')


# *****************************************************************************
# * instrument 'Loki' TRACE
# *****************************************************************************

# Comp instance Origin, placement and parameters
Origin = instr.add_component('Origin','Progress_bar')

Origin.profile = '"NULL"'
Origin.percent = '10'
Origin.flag_save = '0'
Origin.minutes = '0'

# Comp instance virtual_input, placement and parameters
virtual_input = instr.add_component('virtual_input','MCPL_input', AT=['0', '0', '0'], AT_RELATIVE='Origin', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Origin')
# EXTEND at virtual_input
virtual_input.append_EXTEND(r'''
if(SAMPLE_HOLD2 == 2){

SAMPLE_Select = rand() % 2;
SAMPLE_HOLD = SAMPLE_Select;

}
''')



virtual_input.filename = 'virtualoutput_filename'
virtual_input.polarisationuse = '1'
virtual_input.verbose = '1'
virtual_input.Emin = '0'
virtual_input.Emax = 'FLT_MAX'
virtual_input.repeat_count = 'Repeat_Number'
virtual_input.E_smear = '0.01'
virtual_input.pos_smear = '0.0001'
virtual_input.dir_smear = '0.0001'
virtual_input.preload = '0'

# Comp instance l_monitor_before_sample, placement and parameters
l_monitor_before_sample = instr.add_component('l_monitor_before_sample','L_monitor', AT=['0', '0', '0.004'], AT_RELATIVE='virtual_input', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='virtual_input')

l_monitor_before_sample.nL = '100'
l_monitor_before_sample.filename = '"Lambda_before_sample"'
l_monitor_before_sample.nowritefile = '0'
l_monitor_before_sample.xmin = '-0.05'
l_monitor_before_sample.xmax = '0.05'
l_monitor_before_sample.ymin = '-0.05'
l_monitor_before_sample.ymax = '0.05'
l_monitor_before_sample.xwidth = 'DetectorRadius'
l_monitor_before_sample.yheight = 'DetectorRadius'
l_monitor_before_sample.Lmin = '0.0'
l_monitor_before_sample.Lmax = '23.0'
l_monitor_before_sample.restore_neutron = '1'

# Comp instance Hydrated_Disc, placement and parameters
Hydrated_Disc = instr.add_component('Hydrated_Disc','SANSNanodiscsFast_Loki', AT=['0', '0', '0.005'], AT_RELATIVE='virtual_input', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='virtual_input')
# WHEN ( SAMPLE_HOLD == 0 ) at Hydrated_Disc
Hydrated_Disc.set_WHEN('( SAMPLE_HOLD == 0 )')

Hydrated_Disc.NumberOfQBins = '2000'
Hydrated_Disc.AxisRatio = 'AxisRatio'
Hydrated_Disc.NumberOfLipids = 'N_Lipids'
Hydrated_Disc.AreaPerLipidHeadgroup = 'AreaPrHead'
Hydrated_Disc.HeightOfMSP = 'H_1MSP'
Hydrated_Disc.VolumeOfOneMSP = 'V_1MSP'
Hydrated_Disc.VolumeOfHeadgroup = 'V_1Headgruope'
Hydrated_Disc.VolumeOfCH2Tail = 'V_1CH2Tail'
Hydrated_Disc.VolumeOfCH3Tail = 'V_1CH3Tail'
Hydrated_Disc.ScatteringLengthOfOneMSP = 'SLD_1MSP'
Hydrated_Disc.ScatteringLengthOfHeadgroup = 'sigma_Head_Hydrated'
Hydrated_Disc.ScatteringLengthOfCH2Tail = 'sigma_CH2_Hydrated'
Hydrated_Disc.ScatteringLengthOfCH3Tail = 'sigma_CH3_Hydrated'
Hydrated_Disc.Roughness = 'Rough'
Hydrated_Disc.Concentration = 'Conc'
Hydrated_Disc.RhoSolvent = 'RhoSolv'
Hydrated_Disc.RhoSolventIncoherent = 'RhoSolv_Inc'
Hydrated_Disc.AbsorptionCrosssection = 'AbsCrossection'
Hydrated_Disc.AbsorptionCrosssectionSolvent = 'sigma_abs_solvent'
Hydrated_Disc.xwidth = '0.1'
Hydrated_Disc.yheight = '0.1'
Hydrated_Disc.zdepth = '0.005'
Hydrated_Disc.SampleToDetectorDistance = 'Detector_Distance'
Hydrated_Disc.DetectorRadius = 'DetectorRadius'
Hydrated_Disc.qMin = '0.0'
Hydrated_Disc.qMax = '10.0'
Hydrated_Disc.beamwidth_x = '0.01'
Hydrated_Disc.beamwidth_y = '0.01'

# Comp instance Deutorated_Disc, placement and parameters
Deutorated_Disc = instr.add_component('Deutorated_Disc','SANSNanodiscsFast_Loki', AT=['0', '0', '0.005'], AT_RELATIVE='virtual_input', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='virtual_input')
# WHEN ( SAMPLE_HOLD == 1 ) at Deutorated_Disc
Deutorated_Disc.set_WHEN('( SAMPLE_HOLD == 1 )')

Deutorated_Disc.NumberOfQBins = '2000'
Deutorated_Disc.AxisRatio = 'AxisRatio'
Deutorated_Disc.NumberOfLipids = 'N_Lipids'
Deutorated_Disc.AreaPerLipidHeadgroup = 'AreaPrHead'
Deutorated_Disc.HeightOfMSP = 'H_1MSP'
Deutorated_Disc.VolumeOfOneMSP = 'V_1MSP'
Deutorated_Disc.VolumeOfHeadgroup = 'V_1Headgruope'
Deutorated_Disc.VolumeOfCH2Tail = 'V_1CH2Tail'
Deutorated_Disc.VolumeOfCH3Tail = 'V_1CH3Tail'
Deutorated_Disc.ScatteringLengthOfOneMSP = 'SLD_1MSP'
Deutorated_Disc.ScatteringLengthOfHeadgroup = 'sigma_Head_Deutorated'
Deutorated_Disc.ScatteringLengthOfCH2Tail = 'sigma_CH2_Deutorated'
Deutorated_Disc.ScatteringLengthOfCH3Tail = 'sigma_CH3_Deutorated'
Deutorated_Disc.Roughness = 'Rough'
Deutorated_Disc.Concentration = 'Conc'
Deutorated_Disc.RhoSolvent = 'RhoSolv'
Deutorated_Disc.RhoSolventIncoherent = 'RhoSolv_Inc'
Deutorated_Disc.AbsorptionCrosssection = 'AbsCrossection'
Deutorated_Disc.AbsorptionCrosssectionSolvent = 'sigma_abs_solvent'
Deutorated_Disc.xwidth = '0.1'
Deutorated_Disc.yheight = '0.1'
Deutorated_Disc.zdepth = '0.005'
Deutorated_Disc.SampleToDetectorDistance = 'Detector_Distance'
Deutorated_Disc.DetectorRadius = 'DetectorRadius'
Deutorated_Disc.qMin = '0.0'
Deutorated_Disc.qMax = '10.0'
Deutorated_Disc.beamwidth_x = '0.01'
Deutorated_Disc.beamwidth_y = '0.01'

# Comp instance beamstop, placement and parameters
beamstop = instr.add_component('beamstop','Beamstop', AT=['0', '0', '1'], AT_RELATIVE='virtual_input', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='virtual_input')
# WHEN ( Beamstop == 1 ) at beamstop
beamstop.set_WHEN('( Beamstop == 1 )')

beamstop.xmin = '-0.05'
beamstop.xmax = '0.05'
beamstop.ymin = '-0.05'
beamstop.ymax = '0.05'
beamstop.xwidth = '0.012'
beamstop.yheight = '0.012'
beamstop.radius = '0'

# Comp instance tof_sans_monitor, placement and parameters
tof_sans_monitor = instr.add_component('tof_sans_monitor','TOF_SANS_spacial_monitor', AT=['0', '0', 'Detector_Distance + 0.005'], AT_RELATIVE='virtual_input', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='virtual_input')

tof_sans_monitor.nr = '200'
tof_sans_monitor.nt = '1000'
tof_sans_monitor.nq = '500'
tof_sans_monitor.qFilename = '"QDetector"'
tof_sans_monitor.filename = '"TOF_SANS"'
tof_sans_monitor.calibration_file = 'calibration_filename'
tof_sans_monitor.rmax = 'DetectorRadius'
tof_sans_monitor.rmin = 'Rmin'
tof_sans_monitor.tmin = 'T_min'
tof_sans_monitor.tmax = 'T_max'
tof_sans_monitor.SDD = 'Detector_Distance + 0.005'
tof_sans_monitor.restore_neutron = '1'
tof_sans_monitor.instrumentlength = 'Detector_Distance + 0.005 + 23.5'

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


# end of generated Python code Nanodisc_Loki_generated.py 
