#!/usr/bin/env python3
# Automatically generated file. 
# Format:    Python script code
# McStas <http://www.mcstas.org>
# Instrument: reflectometer.instr (reflectometer)
# Date:       Wed Sep  6 07:17:09 2023
# File:       reflectometer_generated.py

import mcstasscript as ms

# Python McStas instrument description
instr = ms.McStas_instr("reflectometer_generated", author = "McCode Py-Generator", origin = "ESS DMSC")

# Add collected DEPENDENCY strings
instr.set_dependency('')

# *****************************************************************************
# * Start of instrument 'reflectometer' generated code
# *****************************************************************************
# MCSTAS system dir is "/Users/pkwi/McStas/mcstas/3.x-dev/"


# *****************************************************************************
# * instrument 'reflectometer' and components DECLARE
# *****************************************************************************

# Instrument parameters:

lambda_min = instr.add_parameter('double', 'lambda_min', value=5.3, comment='Parameter type (double) added by McCode py-generator')
lambda_max = instr.add_parameter('double', 'lambda_max', value=5.45, comment='Parameter type (double) added by McCode py-generator')
slittranslation = instr.add_parameter('double', 'slittranslation', value=0, comment='Parameter type (double) added by McCode py-generator')
sampletranslation = instr.add_parameter('double', 'sampletranslation', value=0, comment='Parameter type (double) added by McCode py-generator')
slitwidth = instr.add_parameter('double', 'slitwidth', value=0.001, comment='Parameter type (double) added by McCode py-generator')
slitheight = instr.add_parameter('double', 'slitheight', value=0.002, comment='Parameter type (double) added by McCode py-generator')
dist_source2slit = instr.add_parameter('double', 'dist_source2slit', value=1, comment='Parameter type (double) added by McCode py-generator')
dist_slit2slit = instr.add_parameter('double', 'dist_slit2slit', value=3.2, comment='Parameter type (double) added by McCode py-generator')
dist_slit2sample = instr.add_parameter('double', 'dist_slit2sample', value=0.18, comment='Parameter type (double) added by McCode py-generator')
dist_sample2detector = instr.add_parameter('double', 'dist_sample2detector', value=2, comment='Parameter type (double) added by McCode py-generator')
sampletype = instr.add_parameter('double', 'sampletype', value=1, comment='Parameter type (double) added by McCode py-generator')
samplesize = instr.add_parameter('double', 'samplesize', value=0.15, comment='Parameter type (double) added by McCode py-generator')
substratethickness = instr.add_parameter('double', 'substratethickness', value=0.003, comment='Parameter type (double) added by McCode py-generator')
MR_Qc = instr.add_parameter('double', 'MR_Qc', value=0.15, comment='Parameter type (double) added by McCode py-generator')
sampleangle = instr.add_parameter('double', 'sampleangle', value=2.5, comment='Parameter type (double) added by McCode py-generator')
detectorangle = instr.add_parameter('double', 'detectorangle', value=5, comment='Parameter type (double) added by McCode py-generator')

component_definition_metadata = {
}
instr.append_declare(r'''
double blocktranslation;
''')


instr.append_initialize(r'''
blocktranslation = -slittranslation;
''')


# *****************************************************************************
# * instrument 'reflectometer' TRACE
# *****************************************************************************

# Comp instance Origin, placement and parameters
Origin = instr.add_component('Origin','Progress_bar')

Origin.profile = '"NULL"'
Origin.percent = '10'
Origin.flag_save = '0'
Origin.minutes = '0'

# Comp instance Source, placement and parameters
Source = instr.add_component('Source','Source_Maxwell_3', AT=['0', '0', '0'], AT_RELATIVE='Origin', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Origin')

Source.size = '0.12'
Source.yheight = '0'
Source.xwidth = '0'
Source.Lmin = 'lambda_min'
Source.Lmax = 'lambda_max'
Source.dist = 'dist_source2slit + dist_slit2slit'
Source.focus_xw = 'slitwidth'
Source.focus_yh = 'slitheight'
Source.T1 = '150.42'
Source.T2 = '38.72'
Source.T3 = '14.84'
Source.I1 = '3.67E11'
Source.I2 = '3.64E11'
Source.I3 = '0.95E11'
Source.target_index = '+ 1'
Source.lambda0 = '0'
Source.dlambda = '0'

# Comp instance Slit1, placement and parameters
Slit1 = instr.add_component('Slit1','Slit', AT=['0', '0', 'dist_source2slit'], AT_RELATIVE='Source', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Source')

Slit1.xmin = 'UNSET'
Slit1.xmax = 'UNSET'
Slit1.ymin = 'UNSET'
Slit1.ymax = 'UNSET'
Slit1.radius = 'UNSET'
Slit1.xwidth = 'slitwidth'
Slit1.yheight = 'slitheight'

# Comp instance Slit2, placement and parameters
Slit2 = instr.add_component('Slit2','Slit', AT=['0', '0', 'dist_slit2slit'], AT_RELATIVE='Slit1', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Slit1')

Slit2.xmin = 'UNSET'
Slit2.xmax = 'UNSET'
Slit2.ymin = 'UNSET'
Slit2.ymax = 'UNSET'
Slit2.radius = 'UNSET'
Slit2.xwidth = 'slitwidth'
Slit2.yheight = 'slitheight'

# Comp instance Arm_sampleNOROTNOTRANS, placement and parameters
Arm_sampleNOROTNOTRANS = instr.add_component('Arm_sampleNOROTNOTRANS','Arm', AT=['blocktranslation', '0', 'dist_slit2sample'], AT_RELATIVE='Slit2', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Slit2')


# Comp instance Arm_sampleNOROT, placement and parameters
Arm_sampleNOROT = instr.add_component('Arm_sampleNOROT','Arm', AT=['sampletranslation', '0', '0'], AT_RELATIVE='Arm_sampleNOROTNOTRANS', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Arm_sampleNOROTNOTRANS')


# Comp instance Arm_sample, placement and parameters
Arm_sample = instr.add_component('Arm_sample','Arm', AT=['0', '0', '0'], AT_RELATIVE='Arm_sampleNOROT', ROTATED=['0', 'sampleangle', '0'], ROTATED_RELATIVE='Arm_sampleNOROT')


# Comp instance Sample_Mirror, placement and parameters
Sample_Mirror = instr.add_component('Sample_Mirror','Mirror', AT=['0', '0', '0'], AT_RELATIVE='Arm_sample', ROTATED=['0', '90', '0'], ROTATED_RELATIVE='Arm_sample')
# WHEN ( sampletype == 1 ) at Sample_Mirror
Sample_Mirror.set_WHEN('( sampletype == 1 )')

Sample_Mirror.reflect = '0'
Sample_Mirror.xwidth = 'samplesize'
Sample_Mirror.yheight = 'samplesize'
Sample_Mirror.R0 = '0.99'
Sample_Mirror.Qc = 'MR_Qc'
Sample_Mirror.alpha = '6.07'
Sample_Mirror.m = '1'
Sample_Mirror.W = '0.003'
Sample_Mirror.center = '1'
Sample_Mirror.transmit = '0'

# Comp instance Sample_Mirror_backside, placement and parameters
Sample_Mirror_backside = instr.add_component('Sample_Mirror_backside','Isotropic_Sqw', AT=['0', '0', '- substratethickness / 2 -1e-6'], AT_RELATIVE='Sample_Mirror', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Sample_Mirror')
# WHEN ( sampletype == 1 ) at Sample_Mirror_backside
Sample_Mirror_backside.set_WHEN('( sampletype == 1 )')

Sample_Mirror_backside.powder_format = '{ 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }'
Sample_Mirror_backside.Sqw_coh = '0'
Sample_Mirror_backside.Sqw_inc = '0'
Sample_Mirror_backside.geometry = '0'
Sample_Mirror_backside.radius = '0'
Sample_Mirror_backside.thickness = '0'
Sample_Mirror_backside.xwidth = 'samplesize'
Sample_Mirror_backside.yheight = 'samplesize'
Sample_Mirror_backside.zdepth = 'substratethickness'
Sample_Mirror_backside.threshold = '1e-20'
Sample_Mirror_backside.order = '0'
Sample_Mirror_backside.T = '0'
Sample_Mirror_backside.verbose = '1'
Sample_Mirror_backside.d_phi = '0'
Sample_Mirror_backside.concentric = '0'
Sample_Mirror_backside.rho = '1 / 13.827'
Sample_Mirror_backside.sigma_abs = '500.08'
Sample_Mirror_backside.sigma_coh = '0'
Sample_Mirror_backside.sigma_inc = '4.935'
Sample_Mirror_backside.classical = '-1'
Sample_Mirror_backside.powder_Dd = '0'
Sample_Mirror_backside.powder_DW = '0'
Sample_Mirror_backside.powder_Vc = '0'
Sample_Mirror_backside.density = '0'
Sample_Mirror_backside.weight = '0'
Sample_Mirror_backside.p_interact = '-1'
Sample_Mirror_backside.norm = '-1'
Sample_Mirror_backside.powder_barns = '1'
Sample_Mirror_backside.quantum_correction = '"Frommhold"'

# Comp instance Sample_Multilayer1, placement and parameters
Sample_Multilayer1 = instr.add_component('Sample_Multilayer1','Mirror', AT=['0', '0', '0'], AT_RELATIVE='Arm_sample', ROTATED=['0', '90', '0'], ROTATED_RELATIVE='Arm_sample')
# WHEN ( sampletype == 2 ) at Sample_Multilayer1
Sample_Multilayer1.set_WHEN('( sampletype == 2 )')

Sample_Multilayer1.reflect = '"d54DMPC-D2O.dat"'
Sample_Multilayer1.xwidth = 'samplesize'
Sample_Multilayer1.yheight = 'samplesize'
Sample_Multilayer1.R0 = '0.99'
Sample_Multilayer1.Qc = '0.021'
Sample_Multilayer1.alpha = '6.07'
Sample_Multilayer1.m = '2'
Sample_Multilayer1.W = '0.003'
Sample_Multilayer1.center = '1'
Sample_Multilayer1.transmit = '0'

# Comp instance Sample_Multilayer2, placement and parameters
Sample_Multilayer2 = instr.add_component('Sample_Multilayer2','Mirror', AT=['0', '0', '0'], AT_RELATIVE='Arm_sample', ROTATED=['0', '90', '0'], ROTATED_RELATIVE='Arm_sample')
# WHEN ( sampletype == 3 ) at Sample_Multilayer2
Sample_Multilayer2.set_WHEN('( sampletype == 3 )')

Sample_Multilayer2.reflect = '"d54DMPC-H2O.dat"'
Sample_Multilayer2.xwidth = 'samplesize'
Sample_Multilayer2.yheight = 'samplesize'
Sample_Multilayer2.R0 = '0.99'
Sample_Multilayer2.Qc = '0.021'
Sample_Multilayer2.alpha = '6.07'
Sample_Multilayer2.m = '2'
Sample_Multilayer2.W = '0.003'
Sample_Multilayer2.center = '1'
Sample_Multilayer2.transmit = '0'

# Comp instance Sample_Multilayer3, placement and parameters
Sample_Multilayer3 = instr.add_component('Sample_Multilayer3','Mirror', AT=['0', '0', '0'], AT_RELATIVE='Arm_sample', ROTATED=['0', '90', '0'], ROTATED_RELATIVE='Arm_sample')
# WHEN ( sampletype == 5 ) at Sample_Multilayer3
Sample_Multilayer3.set_WHEN('( sampletype == 5 )')

Sample_Multilayer3.reflect = '"hDMPC-D2O.dat"'
Sample_Multilayer3.xwidth = 'samplesize'
Sample_Multilayer3.yheight = 'samplesize'
Sample_Multilayer3.R0 = '0.99'
Sample_Multilayer3.Qc = '0.021'
Sample_Multilayer3.alpha = '6.07'
Sample_Multilayer3.m = '2'
Sample_Multilayer3.W = '0.003'
Sample_Multilayer3.center = '1'
Sample_Multilayer3.transmit = '0'

# Comp instance Sample_Multilayer4, placement and parameters
Sample_Multilayer4 = instr.add_component('Sample_Multilayer4','Mirror', AT=['0', '0', '0'], AT_RELATIVE='Arm_sample', ROTATED=['0', '90', '0'], ROTATED_RELATIVE='Arm_sample')
# WHEN ( sampletype == 6 ) at Sample_Multilayer4
Sample_Multilayer4.set_WHEN('( sampletype == 6 )')

Sample_Multilayer4.reflect = '"hDMPC-H2O.dat"'
Sample_Multilayer4.xwidth = 'samplesize'
Sample_Multilayer4.yheight = 'samplesize'
Sample_Multilayer4.R0 = '0.99'
Sample_Multilayer4.Qc = '0.021'
Sample_Multilayer4.alpha = '6.07'
Sample_Multilayer4.m = '2'
Sample_Multilayer4.W = '0.003'
Sample_Multilayer4.center = '1'
Sample_Multilayer4.transmit = '0'

# Comp instance Sample_Multilayer5, placement and parameters
Sample_Multilayer5 = instr.add_component('Sample_Multilayer5','Mirror', AT=['0', '0', '0'], AT_RELATIVE='Arm_sample', ROTATED=['0', '90', '0'], ROTATED_RELATIVE='Arm_sample')
# WHEN ( sampletype == 6 ) at Sample_Multilayer5
Sample_Multilayer5.set_WHEN('( sampletype == 6 )')

Sample_Multilayer5.reflect = '"silicon-D2O.dat"'
Sample_Multilayer5.xwidth = 'samplesize'
Sample_Multilayer5.yheight = 'samplesize'
Sample_Multilayer5.R0 = '0.99'
Sample_Multilayer5.Qc = '0.021'
Sample_Multilayer5.alpha = '6.07'
Sample_Multilayer5.m = '2'
Sample_Multilayer5.W = '0.003'
Sample_Multilayer5.center = '1'
Sample_Multilayer5.transmit = '0'

# Comp instance Sample_Multilayer6, placement and parameters
Sample_Multilayer6 = instr.add_component('Sample_Multilayer6','Mirror', AT=['0', '0', '0'], AT_RELATIVE='Arm_sample', ROTATED=['0', '90', '0'], ROTATED_RELATIVE='Arm_sample')
# WHEN ( sampletype == 7 ) at Sample_Multilayer6
Sample_Multilayer6.set_WHEN('( sampletype == 7 )')

Sample_Multilayer6.reflect = '"silicon-H2O.dat"'
Sample_Multilayer6.xwidth = 'samplesize'
Sample_Multilayer6.yheight = 'samplesize'
Sample_Multilayer6.R0 = '0.99'
Sample_Multilayer6.Qc = '0.021'
Sample_Multilayer6.alpha = '6.07'
Sample_Multilayer6.m = '2'
Sample_Multilayer6.W = '0.003'
Sample_Multilayer6.center = '1'
Sample_Multilayer6.transmit = '0'

# Comp instance Arm_detectorONLYROT, placement and parameters
Arm_detectorONLYROT = instr.add_component('Arm_detectorONLYROT','Arm', AT=['0', '0', '0'], AT_RELATIVE='Arm_sampleNOROTNOTRANS', ROTATED=['0', 'detectorangle', '0'], ROTATED_RELATIVE='Source')


# Comp instance Arm_detector, placement and parameters
Arm_detector = instr.add_component('Arm_detector','Arm', AT=['0', '0', 'dist_sample2detector'], AT_RELATIVE='Arm_detectorONLYROT', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Arm_detectorONLYROT')


# Comp instance Detector, placement and parameters
Detector = instr.add_component('Detector','PSD_monitor', AT=['0', '0', '0'], AT_RELATIVE='Arm_detector', ROTATED=['0.0', '0.0', '0.0'], ROTATED_RELATIVE='Arm_detector')

Detector.nx = '200'
Detector.ny = '200'
Detector.filename = '"mon_detector"'
Detector.xmin = '-0.05'
Detector.xmax = '0.05'
Detector.ymin = '-0.05'
Detector.ymax = '0.05'
Detector.xwidth = '0.025'
Detector.yheight = '0.05'
Detector.restore_neutron = '1'
Detector.nowritefile = '0'

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


# end of generated Python code reflectometer_generated.py 
