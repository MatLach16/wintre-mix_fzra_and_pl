import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime
import os
import subprocess #install to run bash lines from python
import metpy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import matplotlib
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from scipy.optimize import bisect
import microphy_p3_qv_sat

# '''
# Functions written by Mathieu Lachapelle in 2024 to run the Kinematic models form a python interactive interface
# '''


def prep_files_for_wrapper_run(path_canam,path_model,with_vd=False):
    
    paths_to_add_to_kin1d_model = []
    paths_to_add_to_kin1d_model.append( path_canam + 'external/P3/src/microphy_p3.f90')
    paths_to_add_to_kin1d_model.append( path_canam + 'external/P3/interfaces/canam/mp_p3_wrapper_canam.f90')
    paths_to_add_to_kin1d_model.append( path_canam + 'external/P3/interfaces/canam/statcld5_wrapper_p3.f90')
    paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/cond11.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/statcld5.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/zcr.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/thcal4.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'glue/phys_consts.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'glue/phys_parm_defs.F90'   )
    paths_to_add_to_kin1d_model.append( path_canam + 'glue/tracers_info_mod.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'external/PAM/gen/wrn.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'glue/psizes_19.F90')
    paths_to_add_to_kin1d_model.append( path_canam + 'glue/msizes.F90')    
    if path_model!='src-canam':
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/diagnostics/phys_functions.F90')
    if with_vd==True:
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/vrtdftke.F90')
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/vrtdf22.F90')
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/nlclmxtke.F90')
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/nlclmx8.F90')        
        
        paths_to_add_to_kin1d_model.append( path_canam + 'external/P3/interfaces/canam/vrtdf_in_p3_wrapper.f90')
        # paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/tgcal2.F90')
        paths_to_add_to_kin1d_model.append( path_canam + 'glue/times_mod.F90')
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/implvd7.F90')
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/abcvdq6.F90')  
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/abcvdm6.F90')  
        paths_to_add_to_kin1d_model.append( path_canam + 'physics/atmosphere/phys/vrossr.F90')
    for file in paths_to_add_to_kin1d_model:
        subprocess.run(['cp', file, path_model])
        
    '''Modify files in src_canam'''
    code_to_add_to_p3 = path_model+'code_to_add_to_p3.txt'
    insert_code_after_line_from_file(path_model + 'microphy_p3.f90', code_to_add_to_p3, match_line='    endif compute_procs')
    insert_code_after_line(path_model + 'microphy_p3.f90', 'use input_kinematic_model\n', match_line='use statcld5_wrapper_p3, only : statcld5_call_in_p3')

    insert_code_after_line(path_model + 'cond11.F90', '  integer, parameter :: iso4 = 3\n', match_line='  real, parameter :: one = 1.')
    insert_code_after_line(path_model + 'cond11.F90', '  integer, parameter :: iiwc = 2\n', match_line='  real, parameter :: one = 1.')
    insert_code_after_line(path_model + 'cond11.F90', '  integer, parameter :: ilwc = 1\n', match_line='  real, parameter :: one = 1.')

    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: irc = 4\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iso4 = 3\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iiwc = 2\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: ilwc = 1\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: icn = 5\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: irn = 6\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iiwc1 = 7\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iirc1 = 8\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iin1 = 9\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iirv1 = 10\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iiz1 = 11\n', match_line='    implicit none')
    insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    integer, parameter :: iilc1 = 12\n', match_line='    implicit none')

    # remove_lines_with_value(path_model + 'mp_p3_wrapper_canam.f90','    use phys_parm_defs, only : ap_cdnc_fac, ap_cdnc_exp,ap_facacc,ap_facaut,ap_facacc')
    # insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    real, parameter :: ap_cdnc_fac = 60.\n',   match_line=' ! Interactive aerosols')
    # insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    real, parameter :: ap_cdnc_exp = 0.2\n',   match_line=' ! Interactive aerosols')
    # insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    real, parameter :: ap_facaut = 1.\n',   match_line=' ! Interactive aerosols')
    # insert_code_after_line(path_model + 'mp_p3_wrapper_canam.f90', '    real, parameter :: ap_facacc = 1.\n',   match_line=' ! Interactive aerosols')    
    if with_vd==True:
        remove_lines_with_value(path_model + 'vrtdftke.F90', '  use tracers_info_mod, only: modl_tracers, iiwc, ilwc')
        insert_code_after_line(path_model +  'vrtdftke.F90', '   integer, parameter :: iiwc = 2\n', match_line='  implicit none')
        insert_code_after_line(path_model +  'vrtdftke.F90', '    integer, parameter :: ilwc = 1\n', match_line='  implicit none')
        insert_code_after_line(path_model +  'vrtdftke.F90', '      if (n==ilwc .or. n==iiwc) then\n', match_line='      if (modl_tracers(n)%phs/=0 .or. n==ilwc .or. n==iiwc) then')
        remove_lines_with_value(path_model + 'vrtdftke.F90', '      if (modl_tracers(n)%phs/=0 .or. n==ilwc .or. n==iiwc) then')
        remove_lines_with_value(path_model + 'nlclmxtke.F90', '  use tracers_info_mod, only: modl_tracers, iiwc, ilwc')
        insert_code_after_line(path_model +  'nlclmxtke.F90', '   integer, parameter :: iiwc = 2\n', match_line='  implicit none')
        insert_code_after_line(path_model +  'nlclmxtke.F90', '    integer, parameter :: ilwc = 1\n', match_line='  implicit none')
        insert_code_after_line(path_model +  'nlclmxtke.F90', '      if (n==ilwc .or. n==iiwc) then\n', match_line='    if (modl_tracers(n)%phs /=0 .or. n==ilwc .or. n==iiwc) then')
        remove_lines_with_value(path_model + 'nlclmxtke.F90', '    if (modl_tracers(n)%phs /=0 .or. n==ilwc .or. n==iiwc) then')
        remove_lines_with_value(path_model + 'vrtdf_in_p3_wrapper.f90', 'xit(')
        remove_lines_with_value(path_model + 'vrtdftke.F90', 'xit(')
        remove_lines_with_value(path_model + 'vrtdf22.F90', 'xit(')
        # Add tracers to psizes_19
        replace_line(path_model +  'msizes.F90', '  integer,parameter :: ntrac = 4\n', match_line='  integer :: ntrac')
        replace_line(path_model +  'msizes.F90', '  integer,parameter :: ilev=49\n', match_line='  integer :: ilev')
        replace_line(path_model +  'msizes.F90', '  integer,parameter :: im=1\n', match_line='  integer :: im')        
        

    remove_lines_with_value(path_model + 'psizes_19.F90', 'delz_20')
    remove_lines_with_value(path_model + 'psizes_19.F90', 'zbot_20')
    remove_lines_with_value(path_model + 'statcld5.F90', 'xit(')


    remove_lines_with_value(path_model + 'mp_p3_wrapper_canam.f90', 'phys_diag')
    remove_lines_with_value(path_model + 'cond11.F90', 'tracers_info_mod')
    insert_code_after_line(path_model + 'statcld5.F90', ' csigma = 0.2\n', match_line='  csigma = ap_csigma')
    insert_code_after_line(path_model + 'statcld5.F90', ' ap_scale_cvsg = 3.0E-3 \n', match_line='  csigma = ap_csigma')
    remove_lines_with_value(path_model + 'mp_p3_wrapper_canam.f90', 'csigma = ap_csigma')
    
    remove_lines_with_value(path_model + 'mp_p3_wrapper_canam.f90', 'use tracers_info_mod, only: iso4, iiwc, ilwc,')
    remove_lines_with_value(path_model + 'mp_p3_wrapper_canam.f90', 'icn,irc,irn')
    remove_lines_with_value(path_model + 'mp_p3_wrapper_canam.f90', 'iiwc1,iirc1,iin1,iirv1,iiz1,iilc1')
    
    remove_lines_with_value(path_model + 'microphy_p3.f90', 'iparam =')
    remove_lines_with_value(path_model + 'microphy_p3.f90','integer :: iparam')
    # insert_code_after_line( path_model + 'microphy_p3.f90', ' use input_kinematic_model\n', match_line=' subroutine get_cloud_dsd2(qc,nc,mu_c,rho,nu,dnu,lamc,cdist,cdist1,iSCF) ! FRED DEBUG remove CF dependency')

def insert_code_after_line_from_file(target_file, insert_file, match_line='555    continue'):
    """
    Insert content of insert_file after the line in target_file that contains match_line.
    If output_file is None, it will overwrite target_file.
    """
    with open(target_file, 'r') as tf:
        lines = tf.readlines()

    with open(insert_file, 'r') as inf:
        insert_lines = inf.readlines()

    new_lines = []
    inserted = False
    for line in lines:
        new_lines.append(line)
        if match_line in line and not inserted:
            new_lines.extend(insert_lines)
            inserted = True

    if not inserted:
        print(f"Warning: match line '{match_line}' not found.")

    with open(target_file, 'w') as outf:
        outf.writelines(new_lines)
        
def insert_code_after_line(target_file, insert_lines, match_line='555    continue'):
    """
    Insert content of insert_file after the line in target_file that contains match_line.
    If output_file is None, it will overwrite target_file.
    """
    with open(target_file, 'r') as tf:
        lines = tf.readlines()

    new_lines = []
    inserted = False
    for line in lines:
        new_lines.append(line)
        if match_line in line and not inserted:
            new_lines.extend(insert_lines)
            inserted = True

    if not inserted:
        print(f"Warning: match line '{match_line}' not found.")

    with open(target_file, 'w') as outf:
        outf.writelines(new_lines)       

def replace_line(target_file, insert_lines, match_line='555    continue'):
    """
    Insert content of insert_file after the line in target_file that contains match_line.
    If output_file is None, it will overwrite target_file.
    """
    with open(target_file, 'r') as tf:
        lines = tf.readlines()

    new_lines = []
    inserted = False
    for line in lines:
        if match_line in line and not inserted:
            new_lines.extend(insert_lines)
            inserted = True
        else:
            new_lines.append(line)

    if not inserted:
        print(f"Warning: match line '{match_line}' not found.")

    with open(target_file, 'w') as outf:
        outf.writelines(new_lines)       

def replace_string_in_file(path, old, new):
    with open(path, 'r', encoding='utf-8') as f:
        text = f.read()
    text = text.replace(old, new)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(text)

def remove_lines_with_value(input_file, value_to_remove):
    """
    Remove lines from input_file that exactly match value_to_remove.
    If output_file is None, overwrite input_file.
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    cleaned_lines = [line for line in lines if (value_to_remove not in line) ]

    with open(input_file, 'w') as f:
        f.writelines(cleaned_lines)




def reading_outputs(path_model):
    '''
    Reads the output from the column model
    Outputs four pandas dataframe containing the outputs    
    '''
    print('reading1')
    fort9000_1cat = pd.read_csv(path_model+'fort.9000',header=0,skipinitialspace=True,delimiter=r"\s+")
    fort8000_1cat = pd.read_csv(path_model+'fort.8000',header=0,skipinitialspace=True,delimiter=r"\s+")
    sfcprec_1cat = pd.read_csv(path_model+'sfcprec.dat',header=0,skipinitialspace=True,delimiter=r"\s+")
    outQcld_1cat = pd.read_csv(path_model+'outQcld.dat',header=0,skipinitialspace=True,delimiter=r"\s+")
    outNcld_1cat = pd.read_csv(path_model+'outNcld.dat',header=0,skipinitialspace=True,delimiter=r"\s+")
    # cond11_outputs =  pd.read_csv(path_model+'cond11_outputs.dat',header=0,skipinitialspace=True,delimiter=r"\s+")
    # xtrals3d_outputs =  pd.read_csv(path_model+'xtrals3d_outputs.dat',header=0,skipinitialspace=True,delimiter=r"\s+")
    # xtrals1d_outputs =  pd.read_csv(path_model+'xtrals1d_outputs.dat',header=0,skipinitialspace=True,delimiter=r"\s+")

    print('reading2')

    
    columns_to_divide0 = ['qrmlt','qccon','qcevp','qcnuc','qcacc','qcaut','qimlt','qrcol','qccol','qisub','qidep','qrcon','qrevp','qwgrth','qinuc','qcheti', 'qchetc', 'qrhetc', 'qrheti', \
'qrmul','qcshd','qimlt2','qrcol2','qccol2','nrcol2','qisub2','qidep2','qwgrth2','qinuc2', 'qcheti2', 'qchetc2', 'qrhetc2', 'qrheti2','qrmul2','qcshd2','qicol1_2',\
'qrmlt3','qimlt3','qrcol3','qccol3','nrcol3','qisub3','qidep3','qwgrth3','qinuc3', 'qcheti3', 'qchetc3', 'qrhetc3', 'qrheti3','qrmul3','qcshd3','qicol1_3','qicol2_3',\
'qrmlt4','qimlt4','qrcol4','qccol4','nrcol4','qisub4','qidep4','qwgrth4','qinuc4', 'qcheti4', 'qchetc4', 'qrhetc4', 'qrheti4','qrmul4','qcshd4','qicol1_4','qicol2_4','qicol3_4',\
                         ]
    columns_to_divide1 = ['Qc','Qr','Qg','Qi','Qv','Ql','Qg2','Qi2','Ql2','Qg3','Qi3','Ql3','Qg4','Qi4','Ql4','Qsedi_i1','Qsedi_i2','Qsedi_r','Qfluxi1','Qfluxi2','Qfluxr']
    columns_to_divide2 = ['Bg','Bg2','Bg3','Bg4']
    print('OKOKOK')
    for col in columns_to_divide0:
        if col in fort8000_1cat.columns:
            print(col)
            fort8000_1cat[col] = fort8000_1cat[col]/10**8
    
    for col in columns_to_divide1:
        if col in outQcld_1cat.columns:
            print(col)
            outQcld_1cat[col] = outQcld_1cat[col]/10**8
    for col in columns_to_divide2:
        if col in outNcld_1cat.columns:
            outNcld_1cat[col] = outNcld_1cat[col]/10**8
    

    # for col in cond11_outputs.columns:
    #     cond11_outputs[col] = pd.to_numeric(oCond11)
     
    k_to_z = outQcld_1cat[outQcld_1cat.time==outQcld_1cat.time.values[0]][['k','z']]
    
        
    if len(fort8000_1cat)>1:#input_dict['switch_cond11']=='.false.':
        fort8000_1cat['z'] = fort8000_1cat.apply(lambda row: k_to_z[k_to_z.k==row.k].z.values[0],axis=1 )
        fort9000_1cat['z'] = fort9000_1cat.apply(lambda row: k_to_z[k_to_z.k==row.k].z.values[0],axis=1 )
    
        #Delete unused columns
        for col in fort8000_1cat.columns:
            if np.isnan(fort8000_1cat[col]).all():
                del fort8000_1cat[col]
        for col in fort9000_1cat.columns:
            if np.isnan(fort9000_1cat[col]).all():
                del fort9000_1cat[col]
    print('reading3')
    return fort8000_1cat,fort9000_1cat, sfcprec_1cat, outQcld_1cat, outNcld_1cat

def to_netcdf(df):
    # Extract the coordinates and data variables
    coords_orig = {'time': df['time'].value_counts().keys().sort_values(), 'level': df['k'].value_counts().keys().sort_values()}
    coords  = {'time': df['time'].value_counts().keys().sort_values(),     'level': df.iloc[:len(coords_orig['level'])]['k'].values }
    
    data_vars = {key:(('time','level'), df[key].values.reshape(len(coords['time']),len(coords['level']))) for key in df.columns if key not in ['time', 'k']}

    # Create xarray Dataset
    ds = xr.Dataset(data_vars, coords=coords)
    if 'Qi' in df.columns:
        ds = assign_attr_to_outQ(ds)
    elif 'Ni' in df.columns:
        ds = assign_attr_to_outN(ds)
    elif 'qrmlt' in df.columns:
        ds = assign_attr_to_fort8000(ds)
    elif 'nislf' in df.columns:
        ds = assign_attr_to_fort9000(ds)
    elif 'rhc_cs_cond11' in df.columns:
        ds = assign_attr_to_cond11(ds)
    return ds


def assign_attr_to_outQ(ncQ):
    #Add a function when multiple categories is used
    ncQ['time']       = ncQ['time'].assign_attrs(units="s", description="Time since simulation start")
    ncQ['z']           = ncQ['z'].assign_attrs(units="m", description="Height MSL")
    ncQ['level']          = ncQ['level'].assign_attrs(units="", description="Model level number (in kinematic model)")

    ncQ['w']          = ncQ['w'].assign_attrs(units="m s$^{-1}$", description="Vertical air velocity")
    ncQ['p']          = ncQ['p'].assign_attrs(units="Pa", description="Air pressure")
    ncQ['rho']        = ncQ['rho'].assign_attrs(units=r"kg m$^{-3}$", description="Air pressure")
    ncQ['T']          = ncQ['T'].assign_attrs(units="°C", description="Temperature")
    ncQ['tf1']        = ncQ['tf1'].assign_attrs(units="°C", description="Frost point")
    ncQ['Td']         = ncQ['Td'].assign_attrs(units="°C", description="Dew point")
    ncQ['v_ice']      = ncQ['v_ice'].assign_attrs(units=r"kg m$^{-3}$", description="Mass-weighted mean fallspeed, ice category x")

    Qs          = ['Qc','Qr','Qg','Qi','Qv','Ql']
    description = ['Cloud','Rain','Rime mass','Ice','Vapor','Liquid on ice']
    for var,desc in zip(Qs,description):
        ncQ[var] = ncQ[var].assign_attrs(units=r"kg kg$^{-1}$", description=desc+" mixing ratio")
    ncQ = ncQ.transpose('time','level')
    ncQ['Qcond'] = ncQ['Qc']+ncQ['Qi'] 
    
    return ncQ

def assign_attr_to_outN(ncN):
    #Add a function when multiple categories is used
    ncN['time']   = ncN['time'].assign_attrs(units="s", description="Time since simulation start")
    ncN['z'] = ncN['z'].assign_attrs(units="m", description="Height MSL")
    ncN['level']      = ncN['level'].assign_attrs(units="", description="Model level number")

    ncN['Nc']    = ncN['Nc'].assign_attrs(units=r"#$kg^{-1}$", description="Cloud droplet number mixing ratio")
    ncN['Nr']    = ncN['Nr'].assign_attrs(units=r"#$kg^{-1}$",  description="Rain drops number mixing ratio")
    ncN['Ni']    = ncN['Ni'].assign_attrs(units=r"#$kg^{-1}$",  description="Ice particles in category 1 number mixing ration")
    ncN['Bg']    = ncN['Bg'].assign_attrs(units="m$^3$ kg$^{-1}$", description="Rime volume mixing ratio")
    ncN['diag_di']    = ncN['diag_di'].assign_attrs(units="m", description="Mean diameter of ice ")
    ncN['diag_rhoi']  = ncN['diag_rhoi'].assign_attrs(units=r"m s$^{-1}$",description='mass-weighted mean density, ice category x')
    ncN['diag_rhoi2']  = ncN['diag_rhoi2'].assign_attrs(units=r"m s$^{-1}$",description='mass-weighted mean density, ice category 2')
    
    ncN = ncN.transpose('time','level')
    return ncN

def assign_attr_to_fort8000(nc8):
    #fort8000 dataset attributes
    qvar = ['qimlt', 'qrmlt', 'qrcol', 'qccol', 'nrshdr', 'nrcol',
           'qisub', 'qidep', 'qrcon', 'qrevp', 'qwgrth', 'qinuc', 'qcheti',
           'qchetc', 'qrhetc', 'qrheti', 'qrmul', 'qcshd']
    for var in qvar:
        nc8[var] = nc8[var].assign_attrs(units=r"kg s^${-1}$kg$^{-1}$", description=var+" process rate")

    nc8['nimul']  = nc8['nimul'].assign_attrs(units=r"s^${-1}$kg$^{-1}$", description='Number of secondary ice')
    nc8['time']   = nc8['time'].assign_attrs(units="s", description="Time since simulation start")
    nc8['z'] = nc8['z'].assign_attrs(units="m", description="Height MSL")
    nc8['level']      = nc8['level'].assign_attrs(units="", description="Model level number (in P3)")
    nc8 = nc8.transpose('time','level')
    return nc8

def assign_attr_to_fort9000(nc9):
    #fort9000 dataset attributes
    qvar = ['ninuc', 'nimlt', 'nisub', 'nislf', 'nrhetc', 'nrheti',
           'nchetc', 'ncheti', 'nimul', 'nlevp']
    for var in qvar:
        nc9[var] = nc9[var].assign_attrs(units=r"# s^${-1}$kg$^{-1}$", description=var+" number process rate")

    nc9['time']   = nc9['time'].assign_attrs(units="s", description="Time since simulation start")
    nc9['z'] = nc9['z'].assign_attrs(units="m", description="Height MSL")
    nc9['level']      = nc9['level'].assign_attrs(units="", description="Model level number (in P3)")
    nc9 = nc9.transpose('time','level')  
    return nc9

def assign_attr_to_cond11(nccond11):
    nccond11['time']   = nccond11['time'].assign_attrs(units="s", description="Time since simulation start")
    nccond11['z']      = nccond11['z'].assign_attrs(units="m", description="Height MSL")
    nccond11['level']  = nccond11['level'].assign_attrs(units="", description="Model level number (in P3)")
    # nccond11           = nccond11.transpose('level','time')  
    return nccond11

def write_input_file(path_model,input_dict):
    input_dict_default = {'dtadv':60,'switch_cond11':'.false.','isubg':1,'almx_in':600.,'AMPA':2,'AMPB':5,\
                         'wprofile':1,'ttotmin':90,'qitop':0,'qltop':0,'qgtop':0,'bgtop':0,'nitop':0,'n_iceCat':1,'switch_micro':'P3',\
                         'nk_input':49,'dt_p3':30.,'outfreq':60,'skip_advection':'.true.','iparam':1,'switch_int_so4':'.false.','log_predictNc_1d':'.false.',\
                             'liqFrac':'.false.','bimm':2.,'triple_moment':'.false.'}
    for inpt in input_dict:
        input_dict_default[inpt] = input_dict[inpt]    
    
    f = open(path_model+"input_kinematic_model.f90", "w")
    f.write("module input_kinematic_model\n")
    f.write("implicit none\n") 

    f.write("real, parameter :: dtadv = %f \n" % input_dict_default['dtadv'])
    f.write("logical, parameter   :: switch_cond11 = %s \n" % input_dict_default['switch_cond11'])
    f.write("logical, parameter   :: skip_advection = %s \n" % input_dict_default['skip_advection'])
    f.write("character(len=20), parameter   :: switch_micro = \"%s\" \n" % input_dict_default['switch_micro'])
    
    #Turn off statcld5
    # f.write("integer, parameter   :: icvsg = %i \n" % input_dict_default['icvsg'])
    f.write("integer, parameter   :: isubg = %i \n" % input_dict_default['isubg'])
    
    f.write("real, parameter   :: almx_in = %f \n" % input_dict_default['almx_in'])
    f.write("real, parameter   :: AMPA = %f \n" % input_dict_default['AMPA'])
    f.write("real, parameter   :: AMPB = %f \n" % input_dict_default['AMPB'])

    f.write("integer, parameter   :: wprofile = %i \n" % input_dict_default['wprofile'])
    f.write("integer, parameter   :: ttotmin = %i \n" % input_dict_default['ttotmin'])
    f.write("integer, parameter   :: n_iceCat = %i \n" % input_dict_default['n_iceCat'])
    
    f.write("real, parameter :: qitop = %f \n" % input_dict_default['qitop'])
    f.write("real, parameter :: qltop = %f \n" % input_dict_default['qltop'])
    f.write("real, parameter :: qgtop = %f \n" % input_dict_default['qgtop'])
    f.write("real, parameter :: bgtop = %f \n" % input_dict_default['bgtop'])
    f.write("real, parameter :: nitop = %f \n" % input_dict_default['nitop'])
    f.write("integer, parameter :: nk_input = %i \n" % input_dict_default['nk_input'])
    f.write("integer, parameter :: lev_cond11 = %i \n" % (1+input_dict_default['nk_input']))
    f.write("real, parameter :: dt_p3 = %f \n" % input_dict_default['dt_p3'])
    f.write("integer, parameter :: outfreq = %i \n" % input_dict_default['outfreq'])
    f.write("integer, parameter :: iparam = %i \n" % input_dict_default['iparam'])
    f.write("logical, parameter :: switch_int_so4 = %s \n" % input_dict_default['switch_int_so4'])
    f.write("logical, parameter :: log_predictNc_1d = %s \n" % input_dict_default['log_predictNc_1d'])
    f.write("logical, parameter :: liqFrac = %s \n" % input_dict_default['liqFrac'])
    f.write("real, parameter :: bimm = %s \n" % input_dict_default['bimm'])
    f.write("logical, parameter :: trplMomIce = %s \n" % input_dict_default['triple_moment'])

    
    
    f.write("end module input_kinematic_model\n")   
    f.close()

def write_tprofile_file(df_temp_prof,path='/home/rml000/Scripts/kin1d_matlach/soundings/temp_to_use.dat'):
    '''
    Writes the sounding 
    '''
    
    f = open(path, "w")
    f.write("%i\n" % len(df_temp_prof))
    f.write("XX\n")
    f.write('-----------------------------------------------------------------------------\n')
    f.write('PRES   HGHT   TEMP   QVIN   RELH  DWPT   USPD  VSPD     MIXR   DRCT   SKNT   THTA   THTE   THTV\n')
    f.write('    hPa     m      C   kg/kg   %    C     m/s    m/s       g/kg    deg   knot     K      K      K\n')
    f.write('-----------------------------------------------------------------------------\n')

    if ('USPD' not in df_temp_prof.columns):
        df_temp_prof['USPD'] = 0.
        df_temp_prof['VSPD'] = 0.

    for i in range(len(df_temp_prof)):
        temp_prof_i =  df_temp_prof.iloc[i] 
        f.write('  %4.1f   %5.0f   %4.1f  %4.6f   %4.1f   %4.1f    %4.1f    %4.1f\n' % (temp_prof_i.PRES, temp_prof_i.HGHT,\
                                     temp_prof_i.TEMP, temp_prof_i.QV, temp_prof_i.RELH, temp_prof_i.DWPT, temp_prof_i.USPD, temp_prof_i.VSPD))

    f.close()

    
    
def model_run(path_model, input_dict, compilation=True, print_output=True):
    
    write_input_file(path_model,input_dict)
    
    ncat = 1 #Remove this line when it becomes an option
    if compilation==True:
        
        make_process = subprocess.run(['make', 'clean', '-C', path_model])
        make_process = subprocess.run(['make', 'execld', '-C', path_model])
    
    
    make_process = subprocess.Popen('./execld',cwd=path_model,    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True)
    output, _ = make_process.communicate()
    if print_output:
        print('Compiling..'+path_model)
        print('Running..'+path_model)
        print(output)
        print('Run completed')
    fort8000_1cat,fort9000_1cat, sfcprec_1cat, outQcld_1cat, outNcld_1cat  = reading_outputs(path_model)
    
    return fort8000_1cat,fort9000_1cat, sfcprec_1cat, outQcld_1cat, outNcld_1cat






def rh_to_td(t,rh):
    
    #from https://en.wikipedia.org/wiki/Dew_point
    a=6.112 
    b=17.67
    c=243.5
    d=234.5
    gamma_m = np.log(rh/100.*np.exp( (b-t/d)*(t/(c+t)) ))

    
    td = c*gamma_m/(b-gamma_m)
    return td
    
def read_wintremix_sounding(IOP_number, sounding_number):
    '''
    Reads wintremix sounding that were measured in Sorel. The inputs are the IOP # and sounding #.
    Check directory soundings/wintremix/ to see all the available soundings
    '''
    
    path_wintremi_soundings= '/home/rml000/Scripts/kin1d_matlach/soundings/wintremix/'

    def listdir_nohidden(path):
        f_return = []
        for f in os.listdir(path):
            if not f.startswith('.'):
                f_return.append(f)
        return f_return

    list_s = listdir_nohidden(path_wintremi_soundings+'IOP'+str(IOP_number).zfill(2)+'/master/')
    
    temp_prof = pd.read_csv(path_wintremi_soundings+'IOP'+str(IOP_number).zfill(2)+'/master/'+list_s[sounding_number],skiprows=[1,2],header=0, encoding='ISO-8859-1',sep='\s+')
    print(list_s[sounding_number])
    
    temp_prof['PRES'] = temp_prof['Press']
    temp_prof['TEMP'] = temp_prof['Temp']
    temp_prof['DWPT'] = rh_to_td(temp_prof['Temp'], temp_prof['RelHum'])
    temp_prof['HGHT'] = temp_prof['GPM_MSL']

    temp_prof = temp_prof[::20]
    temp_prof = temp_prof[temp_prof.PRES>100.]
    
    return temp_prof



def skew_t(temp_prof):
    # assign units.

    p = temp_prof['PRES'].values * units.hPa
    T = temp_prof['TEMP'].values * units.degC
    Td = temp_prof['DWPT'].values * units.degC
    # wind_speed = df['speed'].values * units.knots
    # wind_dir = df['direction'].values * units.degrees
    # u, v = mpcalc.wind_components(wind_speed, wind_dir)

    # Calculate the LCL
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])

    print(lcl_pressure, lcl_temperature)

    # Calculate the parcel profile.
    parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')

    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(4, 5))
    skew = SkewT(fig)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    skew.plot(p, T, 'r', linewidth=2)
    skew.plot(p, Td, 'g', linewidth=2)
    
    if ('USPD' in  temp_prof):
        skew.plot_barbs(p[::10], temp_prof['USPD'][::10]*units.ms, temp_prof['VSPD'][::10]*units.ms)
    
    skew.ax.set_xlim(-50, 40)

    # Plot LCL temperature as black dot
    # skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

    # Plot the parcel profile as a black line
    # skew.plot(p, parcel_prof, 'k', linewidth=2)

    # Shade areas of CAPE and CIN
    skew.shade_cin(p, T, parcel_prof, Td)
    skew.shade_cape(p, T, parcel_prof)
    return fig

def figure_comparison(ncQ,ncQ_p3,sfcprec,sfcprec_p3):
    fig,ax=plt.subplots(3,1,figsize=(7,6))
    
    '''Forcing''' 
    isotherm = np.arange(-60,60,20)
    ncQ['T'].T.plot(ax=ax[0],cbar_kwargs={'pad':0.18},cmap='Reds')
    # ct=ncQ['T'].T.plot.contour(levels=isotherm,colors='C0',ax=ax[0],label='Cond11')
    # plt.clabel(ct,fmt='%2.0f C')
    # ct=ncQ_p3['T'].T.plot.contour(levels=isotherm,colors='k',ax=ax[0],label='P3')
    # plt.clabel(ct,fmt='%2.0f C')
    ax[0].plot(1,1,color='k',label='P3')
    ax[0].plot(1,1,color='C0',label='Cond11')
    ax[0].legend(loc='upper right',bbox_to_anchor= (1.26,1),fontsize=8)
    
    ax[0].set_xlabel('')
    ax[0].set_title('Temperature forcing')
    
    to_plot1 =  ( (ncQ.Qc+ncQ.Qi)*ncQ.rho*1000. )
    to_plot2 = ( (ncQ_p3.Qc+ncQ_p3.Qr+ncQ_p3.Qi)*ncQ_p3.rho*1000.  )
    maxvar = np.max([to_plot1.max(),to_plot2.max()])


    ( (ncQ.Qc+ncQ.Qi)*ncQ.rho*1000. ).T.plot(cmap='Blues',ax=ax[1],cbar_kwargs={'pad':0.18,'label':'Total water content\n'+r'[g m$^{-3}$]'},vmin=0,vmax=maxvar)
    tax = ax[1].twinx()
    tax.plot(sfcprec.index,sfcprec.prt_liq+sfcprec.prt_sol,'k',label='Total')
    tax.plot(sfcprec.index,sfcprec.prt_liq,'green',label='Rain')
    tax.plot(sfcprec.index,sfcprec.prt_sol,'blue',label='Snow')
    tax.legend(loc='upper right',fontsize=8,ncol=3)
    tax.set_ylabel('Precipitation rate\n[mm h$^{-1}$]')
    tax.set_ylim(0,1.1*np.max( [np.max(sfcprec_p3.prt_liq+sfcprec_p3.prt_sol),  np.max(sfcprec.prt_liq+sfcprec.prt_sol)] ))
    

    ax[1].set_xlabel('')
    ax[1].set_title('cond11')
    
    ax[1].text(0.01, 0.9, 'PR tot: %2.1e mm' % sfcprec.PR.iloc[-1],
        horizontalalignment='left',
        verticalalignment='center',
        transform=ax[1].transAxes)


    ( (ncQ_p3.Qc+ncQ_p3.Qr+ncQ_p3.Qi)*ncQ_p3.rho*1000.  ).T.plot(cmap='Blues',ax=ax[2],cbar_kwargs={'pad':0.18,'label':'Total water content\n'+r'[g m$^{-3}$]'},vmin=0,vmax=maxvar)

    tax = ax[2].twinx()
    tax.plot(sfcprec_p3.index,sfcprec_p3.prt_liq+sfcprec_p3.prt_sol,'k',label='Total')
    tax.plot(sfcprec_p3.index,sfcprec_p3.prt_liq,'green',label='Rain')
    tax.plot(sfcprec_p3.index,sfcprec_p3.prt_sol,'blue',label='Snow')
    tax.legend(loc='upper right',fontsize=8,ncol=3)
    tax.set_ylabel('Precipitation rate\n[mm h$^{-1}$]')
    tax.set_ylim(0,1.1*np.max( [np.max(sfcprec_p3.prt_liq+sfcprec_p3.prt_sol),  np.max(sfcprec.prt_liq+sfcprec.prt_sol)] ))

    
    ax[2].text(0.01, 0.9, 'PR tot: %2.1e mm' % sfcprec_p3.PR.iloc[-1],
        horizontalalignment='left',
        verticalalignment='center',
        transform=ax[2].transAxes)

    ax[0].set_ylim([49,0])
    ax[1].set_ylim([49,0])
    ax[2].set_ylim([49,0])

    ax[2].set_title('P3')
    plt.subplots_adjust(hspace=0.38)
    
    # plt.savefig('/home/rml000/Scripts/figures/'+savename,  format='png',dpi=200, bbox_inches = 'tight')

def iwc_lwc(ncQ,ncQ_p3,titles=['Def.','Mod.']):

    isotherm = np.arange(-60,60,20)

    fig,ax=plt.subplots(4,1)
    
    ncQ    = ncQ.assign_coords(level=ncQ.z.mean(axis=1).values)
    ncQ_p3 = ncQ_p3.assign_coords(level=ncQ_p3.z.mean(axis=1).values)

    liquid_content = ((ncQ.Qc+ncQ.Qr)*ncQ.rho*1000.)
    liquid_content = liquid_content.where(liquid_content>10**-2)
    solid_content  = ((ncQ.Qi)*ncQ.rho*1000.)
    solid_content  = solid_content.where(solid_content>10**-2)
    
    
    liquid_content_p3 = ((ncQ_p3.Qc+ncQ_p3.Qr)*ncQ_p3.rho*1000.)
    liquid_content_p3 = liquid_content_p3.where(liquid_content_p3>10**-2)
    solid_content_p3  = ((ncQ_p3.Qi)*ncQ_p3.rho*1000.)
    solid_content_p3 = solid_content_p3.where(solid_content_p3>10**-2)
    
    
    vmin = 0
    vmax = 1
    # vmax = np.max([np.max(liquid_content),np.max(solid_content),np.max(liquid_content_p3),np.max(solid_content_p3)])
    
    
    liquid_content.plot(cmap='jet',ax=ax[0],cbar_kwargs={'pad':0.02,'label':titles[0]+' LWC\n'+r'[g m$^{-3}$]'},vmin=vmin,vmax=vmax,norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax),extend='both')
    solid_content.plot(cmap='jet',ax=ax[1],cbar_kwargs={'pad':0.02,'label':titles[0]+' IWC\n'+r'[g m$^{-3}$]'},vmin=vmin,vmax=vmax,norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax),extend='both')
    liquid_content_p3.plot(cmap='jet',ax=ax[2],cbar_kwargs={'pad':0.02,'label':titles[1]+' LWC\n'+r'[g m$^{-3}$]'},vmin=vmin,vmax=vmax,norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax),extend='both')
    solid_content_p3.plot(cmap='jet',ax=ax[3],cbar_kwargs={'pad':0.02,'label':titles[1]+' IWC\n'+r'[g m$^{-3}$]'},vmin=vmin,vmax=vmax,norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax),extend='both')
    
    ct=ncQ['T'].plot.contour(levels=isotherm,colors='k',ax=ax[0],linewidths=0.5)
    plt.clabel(ct,fmt='%2.0f C')
    ct=ncQ['T'].plot.contour(levels=isotherm,colors='k',ax=ax[1],linewidths=0.5)
    plt.clabel(ct,fmt='%2.0f C')
    ct=ncQ_p3['T'].plot.contour(levels=isotherm,colors='k',ax=ax[2],linewidths=0.5)
    plt.clabel(ct,fmt='%2.0f C')
    ct=ncQ_p3['T'].plot.contour(levels=isotherm,colors='k',ax=ax[3],linewidths=0.5)
    plt.clabel(ct,fmt='%2.0f C')


    ax[0].set_ylim(0,13000)
    ax[1].set_ylim(0,13000)
    ax[2].set_ylim(0,13000)
    ax[3].set_ylim(0,13000)

    ax[0].set_yticks([0,5000,10000])
    ax[1].set_yticks([0,5000,10000])
    ax[2].set_yticks([0,5000,10000])
    ax[3].set_yticks([0,5000,10000])

    ax[0].set_xlabel('')
    ax[1].set_xlabel('')
    ax[2].set_xlabel('')

    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    ax[2].set_xticklabels([])

    plt.subplots_adjust(hspace=0.2)

    # plt.savefig('/home/rml000/Scripts/figures/'+savename,  format='png',dpi=200, bbox_inches = 'tight')

def plot_diff(nctot1,nctot2,var,reldiff=True,titles=['nc1','nc2']):

    fig,ax = plt.subplots(3,1,figsize=(5,7))
    
    maxvar = np.max([nctot1[var].max(),nctot2[var].max()])
    
    nctot1[var].T.plot(ax=ax[0],vmin=0,vmax=maxvar,cmap='jet')
    ax[0].set_title(titles[0])
    nctot2[var].T.plot(ax=ax[1],vmin=0,vmax=maxvar,cmap='jet')
    ax[1].set_title(titles[1])
    if reldiff:
        ((nctot2[var]-nctot1[var])/nctot1[var]*100.).T.plot(ax=ax[2],cmap='bwr',vmin=-20,vmax=20)
        ax[2].set_title('100*('+titles[1]+' - ' + titles[0] + ')/P3 ')
    else:
        (nctot2[var]-nctot1[var]).T.plot(ax=ax[2],cmap='bwr',vmin=-np.abs(nctot2[var]-nctot1[var] ).max(),vmax=np.abs(nctot2[var]-nctot1[var] ).max())
        ax[2].set_title(titles[1]+' - '+titles[0])        

    ax[0].set_ylim([49,0])
    ax[1].set_ylim([49,0])
    ax[2].set_ylim([49,0])
    
    plt.subplots_adjust(hspace=0.45)


def ew(T):
    # '''
    # Calcul de la pression de vapeur saturante par rapport à l'eau liquide
    # from https://www.npl.co.uk/resources/q-a/how-do-i-convert-between-units-of-dew-point-and-re
    # '''
    # A0= - 6096.9385
    # A1=   21.2409642
    # A2= - 2.711193*10**(-2)
    # A3=   1.673952*10**(-5)
    # A4=   2.433502

    # return A0*T**(-1) + A1 + A2*T + A3*T**2 +  A4*np.log(T) 

    # As in CanAM
    rw1 = 53.67957
    rw2 = - 6743.769
    rw3 = - 4.8451    
    esw = np.exp(rw1 + rw2 / T) * T**rw3
    return esw

def ei(T):
    # '''
    # Calcul de la pression de vapeur saturante par rapport à la glace
    # from https://www.npl.co.uk/resources/q-a/how-do-i-convert-between-units-of-dew-point-and-re
    # '''    
    # B0 = - 6024.5282
    # B1 =   29.32707
    # B2 =   1.0613868*10**(-2)
    # B3 = - 1.3198825*10**-5
    # B4 = - 0.49382577
    
    # return B0*T**(-1) + B1 + B2*T + B3*T**2 +  B4*np.log(T) 
    ri1 = 23.33086
    ri2 = - 6111.72784
    ri3 = 0.15215
    esi = np.exp(ri1 + ri2 / T) * T**ri3
    return esi

def rh_to_ice(rh,T,nanifneg=True):
    if ((T>273.15) and nanifneg):
        return np.nan
    else:
        return rh * ew(T)/ei(T)



def calcul_de_Td_sature_par_rapport_a_la_glace(T, p):
    '''
    Calcule Td lorsque T<0 et que RH_ice=100%
    Par definition :
    ws(Td)=w(T)
    ew(Td)/(p-ew) = e(T)/(p-e)
    A saturation par rapport a la glace, e=ei.
    => ew(Td)/(p-ew) = ei(T)/(p-ei)
    => 0=ew(Td)*(A+1)-A*p ou A=ei/(p-ei) 
    Td est donc la température qui résout l'inéquation ci-haut.
    '''

    A = ei(T)/(p-ei(T))
    def fonction_zero(T):
        return ew(T)-A*p/(1+A)

    Td = bisect(fonction_zero, 173, 373)

    return Td


def modify_temp_saturated_ice_liq(T,p):
    
    # if T>=273.15:
        # Saturated liqs
    Td = T
    # else:
    #     Td = calcul_de_Td_sature_par_rapport_a_la_glace(T, p)
    
    return Td
     

def qv_from_rh(rh,T,pd):
    esw = microphy_p3_qv_sat.microphy_qv_sat.polysvp1(T,0)

    return  0.622*np.max([0.99,rh])*esw/pd #rh/(1./esw*pd*Rv/Rd-rh)
    
def dewpoint_from_rh(T, RH):
    # FROM CHATGPT
    """
    Dew point from temperature and relative humidity.

    Parameters
    ----------
    T : float or array
        Air temperature (K)
    RH : float or array
        Relative humidity (0-1)

    Returns
    -------
    Td : float or array
        Dew point temperature (K)
    """
    Tc = T - 273.15
    gamma = np.log(RH) + (17.67 * Tc) / (243.5 + Tc)
    Td = 243.5 * gamma / (17.67 - gamma)
    return Td + 273.15

def pr_to_N_Q_2(prec_rate_mmh,all_drop_same_size=True):
    prec_rate = prec_rate_mmh/3600./1000.
    v    = 1
    rhow = 997
    p    = 600 * 100
    T    = -5
    T0K  = 273.15
    rho_air = p/287.058/(T+T0K) 

    # Gunn and Mashal (1958), p. 59 de Pruppacher and Klett
    # lambda_D = 25.5 * 1000 * prec_rate_mmh**(-0.48) #mm^-1 -> m^-1
    # if prec_rate_mmh>0.:
    #     prec_rate_mmh_forced=2
    #     lambda_D = 25.5 * 1000 * prec_rate_mmh_forced**(-0.48) #mm^-1 -> m^-1
    # else:
    #     prec_rate_mmh_forced = 1
    #     lambda_D = 0.
    
    rho_i = 200 # kg/m^3     
    if all_drop_same_size:
        Qi = prec_rate*rhow/rho_air/v
        Di   = 0.001 # Targeted mean size for initialized snow
        Ni = Qi * 6 / (Di**3. * rho_i * np.pi)
    else:
        if prec_rate_mmh>0.:
            # Gunn and Mashal (1958), p. 59 de Pruppacher and Klett
            rho_i = 1000.
            N0 = 3.8 * 10**3 * prec_rate_mmh**(-0.87) * 1000. #m^-3 mm^-1 -> m^-4
            lambda_D = 25.5  * prec_rate_mmh**(-0.48) * 100. #cm^-1 -> m^-1 ### IMPORTANT: there is an error in Pruppacher and Klett (1978) p. 59, the units are cm^-1
            # Integral of N0*rho_i*4*pi/3*(D/2)^3 * exp(-lambda*D) dD from 0 to inf = N0*rho_i*pi/6 * 6/lambda^4 = N0*rho_i/lambda^4
            Qi = N0 * rho_i * np.pi / lambda_D**4  / rho_air
            Ni = N0 / lambda_D  /  rho_air
        else:
            Qi = 0.
            Ni = 0.
    
    return Qi, Ni

def run_with_prof(prate, qgtop_ratio=0,input_dict_mod={},path_model='p31d_sip-ffd-mod-notempdep',all_drop_same_size=True):
    # Write the initial vertical profile into soundings/temp_to_use.dat file
    '''
    Available inputs:
    dt :           Time step in seconds
    switch_cond11: If True, cond11 (CanESM microphysics scheme is run instead of P3)
    icvsg, isubg : Cond11 inputs to control convection
    AMPA         : Amplitude of the updraft
    '''
    qitop,nitop = pr_to_N_Q_2(prate,all_drop_same_size=all_drop_same_size)
    input_dict = {'dt_p3':5,'switch_cond11':'.false.','AMPA': 0.,'wprofile':2,'ttotmin':180,'qitop':qitop,'nitop':nitop,'qgtop':qgtop_ratio*qitop,'n_iceCat':1,'nk_input':102,'outfreq':1} # 1 for sinus in time, 2 for constant in time
    for input in input_dict_mod:
        input_dict[input] = input_dict_mod[input]
    path_model = path_model
    fort8000_p3,fort9000_p3, sfcprec_p3, outQ_p3, outN_p3 = model_run(path_model,input_dict) # Compile 

    ncQ_p3 = to_netcdf(outQ_p3)     # Output containing mass mixing ratio of the different species, and atmosphere conditions profile
    ncN_p3 = to_netcdf(outN_p3)     # Output containing number mixing ratio of the different species, as well as riming and ice density and ice size 
    nc8_p3 = to_netcdf(fort8000_p3) # Output containing mass rates of processes
    nc9_p3 = to_netcdf(fort9000_p3) # Output containing number rates of processes
    return sfcprec_p3,ncQ_p3,ncN_p3,nc8_p3,nc9_p3

