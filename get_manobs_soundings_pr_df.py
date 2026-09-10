import pandas as pd
import numpy as np 
import glob
import matplotlib.pyplot as plt
from metpy.plots import Hodograph, SkewT
import metpy.calc as mpcalc
from metpy.units import units, concatenate
import math
from datetime import time
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import matplotlib.ticker as ticker
import xarray as xr
from scipy.ndimage import gaussian_filter1d
from scipy import integrate
import os
import matplotlib.patches as mpatches
# import seaborn as sns
from scipy.interpolate import interp1d
import re
from collections import defaultdict
from matplotlib.lines import Line2D
import copy
import soundings_functions
import MetData_functions
import parsivel_functions
import sys
import kin1d_functions


def add_interp_sounding_stats_to_manobs(man_obs_with_two_files,df_all_events,use_Tw=True):
    man_obs_with_two_files = man_obs_with_two_files.sort_values(by=['site','Date & time (UTC)'])
    last_site='blou'
    # Find the closest sounding before and after manobs
    for cnt in range(len(man_obs_with_two_files)):#np.arange(1450,1452):
        two_files =  'None'
        row = man_obs_with_two_files.iloc[cnt]
        if row.site!=last_site:
            last_site = row.site
    
        df_all_events_for_site =  df_all_events[df_all_events.site==row['site']]
        if df_all_events_for_site.empty:
            print('No soundings for : ' + row['site'])
            continue
        try:
            closest_sounding_index = (row['Date & time (UTC)'] - df_all_events_for_site['dates']).abs().argmin()
        except:
            print(f"Error occurred while processing row: {row['Date & time (UTC)']}")
            print(row)
            raise Exception('ERROR')
        row_closest_sounding   =  df_all_events_for_site.iloc[closest_sounding_index]
        
        if ((row_closest_sounding['dates']<=row['Date & time (UTC)']) & ((closest_sounding_index+1)<len(df_all_events_for_site))):
            next_row = df_all_events_for_site.iloc[closest_sounding_index+1]
            if ( (next_row.site==row_closest_sounding.site) & (np.abs(next_row['dates']-row_closest_sounding['dates'])<timedelta(hours=6) ) ):
                two_files = [row_closest_sounding.file, next_row.file]
            else:
                two_files = 'None'
        elif (row_closest_sounding['dates']>row['Date & time (UTC)']):
            last_row = df_all_events_for_site.iloc[closest_sounding_index-1]
            if ( (last_row.site==row_closest_sounding.site) & (np.abs(last_row['dates']-row_closest_sounding['dates'])<timedelta(hours=6) ) ):
                two_files = [last_row.file,               row_closest_sounding.file]
            else:
                two_files = 'None'
        
        row_for_interp = pd.DataFrame({'dates':[row['Date & time (UTC)']],'file':[two_files]},index=[0])
        row_for_interp = row_for_interp.iloc[0]
        if two_files != 'None':
            file1_long, file2_long = two_files[0],two_files[1]
            
            file1 = 'soundings_wintremix/'+ file1_long.split('/')[-1]
            file2 = 'soundings_wintremix/'+ file2_long.split('/')[-1]
    
            if "Sorel" in file1:
                df1 = soundings_functions.read_soundings_data_Sorel(file1)
                df2 = soundings_functions.read_soundings_data_Sorel(file2)
                # df  = interpolate_between_soundings(row_for_interp,df1,df2,df_all_events)
            elif "Gault" in file1:
                df1 = soundings_functions.read_soundings_data_Gault(file1)
                df2 = soundings_functions.read_soundings_data_Gault(file2)
                # df  = interpolate_between_soundings(row_for_interp,df1,df2,df_all_events)
            elif "Trois" in file1 : 
                df1 = soundings_functions.read_soundings_data_TR(file1)
                df2 = soundings_functions.read_soundings_data_TR(file2)
                # df  = interpolate_between_soundings(row_for_interp,df1,df2,df_all_events)
            elif "Albany" in file1:
                df1 = soundings_functions.read_soundings_data_Albany(file1)
                df2 = soundings_functions.read_soundings_data_Albany(file2)
                # df  = interpolate_between_soundings(row_for_interp,df1,df2,df_all_events)
            elif "CU" in file1:
                df1 = soundings_functions.read_soundings_data_DOW(file1)
                df2 = soundings_functions.read_soundings_data_DOW(file2)
                # df  = interpolate_between_soundings(row_for_interp,df1,df2,df_all_events)


            df  = interpolate_between_soundings(row_for_interp,df1,df2,df_all_events)
            df['RelHumIce'] = df.apply(lambda row: 100.*kin1d_functions.rh_to_ice(row['RelHum']/100., row['Temp']+273.15,nanifneg=False), axis=1)
            
            man_obs_with_two_files.loc[ man_obs_with_two_files.index[cnt],'first_temp'] = df.Temp[1]
            
            output_sounding_analysis = stats_on_soundings(df,use_Tw=use_Tw)
            for output in output_sounding_analysis:
                man_obs_with_two_files.loc[ man_obs_with_two_files.index[cnt],output] = output_sounding_analysis[output]
                
            man_obs_with_two_files.loc[ man_obs_with_two_files.index[cnt],'two_files'] = file1+','+file2
            
    return man_obs_with_two_files
    
def interpolate_between_soundings(row,df1,df2,df_all_events):
    # #Interpolate
    f_int_t = interp1d(df2.height, df2.Temp,bounds_error=False, fill_value='extrapolate')
    f_int_tw = interp1d(df2.height, df2.Tw,bounds_error=False, fill_value='extrapolate')
    f_int_rh = interp1d(df2.height, df2.RelHum,bounds_error=False, fill_value='extrapolate')
    df2_inter = pd.DataFrame({'height':df1.height,'Temp':f_int_t(df1.height),'Tw':f_int_tw(df1.height),\
        'RelHum':f_int_rh(df1.height)})
    
    
    time1 = df_all_events[df_all_events.file==row.file[0]].dates.values[0]
    time_to_inter = row.dates
    time2 = df_all_events[df_all_events.file==row.file[1]].dates.values[0]


    a = (time_to_inter-time1)/(time2-time1)
    b = (time2-time_to_inter)/(time2-time1)

    df = b*df1[['Temp','height','Tw','RelHum']] + a*df2_inter[['Temp','height','Tw','RelHum']]

    return df
    

def stats_on_soundings(df_profile,use_Tw=True):
    df_profile = df_profile.drop( df_profile[np.isnan(df_profile['height'])].index,axis=0)
    df_profile = df_profile.sort_values(by='height')
    
    if use_Tw:
        T_var = 'Tw'
    else:
        T_var = 'Temp'
    
    heights = np.array(df_profile['height'])
    temps   = np.array(df_profile[T_var])
    rhs     = np.array(df_profile['RelHum'])
    rhis    = np.array(df_profile['RelHumIce'])
    if heights[0] < heights[-1]:
        heights = heights[::-1]
        temps = temps[::-1]
        rhs = rhs[::-1]
        rhis = rhis[::-1]
        
    dz = 5
    znew = np.arange(np.ceil(np.min(heights)), 5001, dz)
    #Interpolate
    f_pres = interp1d(heights, temps,bounds_error=False)
    temps = f_pres(znew)
    f_pres_rh = interp1d(heights, rhs,bounds_error=False)
    rhs = f_pres_rh(znew) 
    f_pres_rhis = interp1d(heights, rhis,bounds_error=False)
    rhis = f_pres_rhis(znew)
    
    heights = znew[::-1]
    temps = temps[::-1]
    rhs = rhs[::-1]
    rhis = rhis[::-1]
    
    
    tops_WL_idx = []
    bots_WL_idx = []
    tops_CL_idx = []
    bots_CL_idx = []
    last_temp = temps[0]
    if (temps[0]>0.):
        print('First temp positif, stop analysis') #Don't treat anything if the first temperature is positive
    else:
        for idx, (hgt, temp) in enumerate(zip(heights, temps)):
            if ((last_temp<=0.) and (temp>0.)):
                tops_WL_idx.append(idx)
                if len(tops_CL_idx)>len(bots_CL_idx): # We only add a bottom if a top exists
                    bots_CL_idx.append(idx-1)
            elif ((last_temp>0.) and (temp<=0.)): 
                if len(tops_WL_idx)>len(bots_WL_idx):
                    bots_WL_idx.append(idx-1) 
                tops_CL_idx.append(idx)
            last_temp = temp
            if (idx==(len(temps)-1)):# Ground is reached
                if (temp<=0.):
                    bots_CL_idx.append(idx)
                elif (temp>0.): 
                    bots_WL_idx.append(idx)
    
    if ((len(bots_WL_idx)==len(tops_WL_idx)>=1) & (len(bots_CL_idx)==len(tops_CL_idx)>=1 ) ):     
        depths_wl = []
        for i in range(len(tops_WL_idx)):
            depths_wl.append(heights[tops_WL_idx[i]]-heights[bots_WL_idx[i]])
        depths_cl = []
        for i in range(len(tops_CL_idx)):
            depths_cl.append(heights[tops_CL_idx[i]]-heights[bots_CL_idx[i]])
        
        idx_deepest_wl = np.argmax(depths_wl)
        idx_deepest_cl = np.argmax(depths_cl )         
        # if tops_CL_idx[idx_deepest_cl]
        
        
        wl_indices = np.arange(tops_WL_idx[idx_deepest_wl], bots_WL_idx[idx_deepest_wl]+1 )
        cl_indices = np.arange(tops_CL_idx[idx_deepest_cl], bots_CL_idx[idx_deepest_cl]+1 )
        
        wl_depth_hght = heights[tops_WL_idx[idx_deepest_wl]] - heights[bots_WL_idx[idx_deepest_wl]]
        cl_depth_hght = heights[tops_CL_idx[idx_deepest_cl]] - heights[bots_CL_idx[idx_deepest_cl]]
        max_WL_temp = np.max(temps[wl_indices])
        min_CL_temp = np.min(temps[cl_indices])
        WL_max_hght = np.max(heights[wl_indices])
        CL_min_hght = np.min(heights[cl_indices])
        cl_base_hght = np.min(heights[bots_CL_idx[idx_deepest_cl]])
        WL_mean_temp = np.mean(temps[wl_indices])
        CL_mean_temp = np.mean(temps[cl_indices])
        bot_cl_hght = heights[bots_CL_idx[idx_deepest_cl]]
        top_cl_hght = heights[tops_CL_idx[idx_deepest_cl]]
        bot_WL_hght = heights[bots_WL_idx[idx_deepest_wl]]
        top_WL_hght = heights[tops_WL_idx[idx_deepest_wl]]
    
        WL_integral = np.sum(temps[wl_indices]*dz)
        # print(np.max(temps[wl_indices]-273.15))
        # print(len(wl_indices))
        # plt.plot(temps[wl_indices],heights[wl_indices],'k--')
        CL_integral = np.sum(temps[cl_indices]*dz)
        # print(len(cl_indices))
        # plt.plot(temps[cl_indices],heights[cl_indices],'b--')
        # plt.show()
        
        CL_max_rhi = np.nanmax(rhis[cl_indices])
        CL_max_rh = np.nanmax(rhs[cl_indices])
        
        
    else:
        wl_depth_hght = 0.0
        cl_depth_hght = np.nan
        max_WL_temp = np.nan
        min_CL_temp = np.nan
        WL_max_hght = np.nan
        CL_min_hght = np.nan
        cl_base_hght = np.nan
        WL_mean_temp = np.nan
        CL_mean_temp = np.nan
        wl_indices = np.nan
        cl_indices = np.nan
        bot_cl_hght = np.nan
        top_cl_hght = np.nan
        bot_WL_hght = np.nan
        top_WL_hght = np.nan
    
        WL_integral = np.nan
        CL_integral = np.nan
        CL_max_rhi = np.nan
        CL_max_rh = np.nan
        
    df_to_return =  {'wl_depth_hght':wl_depth_hght, 'cl_depth_hght':cl_depth_hght, 'cl_base_hght':cl_base_hght,\
                'max_WL_temp':max_WL_temp, 'min_CL_temp':min_CL_temp, 'CL_min_hght':CL_min_hght,\
                'WL_max_hght':WL_max_hght,'WL_depth_1':np.nan,'WL_depth_2':np.nan,\
                'WL_mean_temp':WL_mean_temp,'CL_mean_temp': CL_mean_temp,\
                'bot_cl_hght':bot_cl_hght, 'top_cl_hght':top_cl_hght, 'bot_WL_hght':bot_WL_hght,'top_WL_hght':top_WL_hght,\
                'WL_integral':WL_integral,'CL_integral':CL_integral,'CL_max_rhi':CL_max_rhi,'CL_max_rh':CL_max_rh}
    return df_to_return



def get_all_event_df():
    
    path = 'man_obs_wintremix/'
    
    # Liste pour stocker tous les chemins de fichiers
    files = []
    
    # Utiliser os.walk pour parcourir tous les sous-répertoires
    for root, dirs, files_in_dir in os.walk(path):
        for nom_fichier in files_in_dir:
            chemin_complet = os.path.join(root, nom_fichier)
            # Ajouter le fichier à la liste si ce n'est pas un fichier système (comme .DS_Store)
            if ((".DS_Store" not in chemin_complet) & ('checkpoint' not in chemin_complet)):
                files.append(chemin_complet)
    
    #Open and concat all files
    dfs = []
    for file in files : 
        
        # Trouver le nom du site 
        site = file.split('.')[-2]
        if site not in ['ESSX', 'COW-CAN']:
            IOP = file.split('.')[2]
            df = pd.read_csv(file)
            
            #Ajouter une colonne avec le nom du site
            df['site'] = site
            df['IOP']  = IOP
            #Sélectionner les colonnes pertinentes 
            columns = ['Date & time (UTC)', 'Primary p-type', 'Secondary p-type', 'site','IOP']
            df = df[columns]
        
            #Mettre la colonne de Date en Datetime
            df['Date & time (UTC)'] = pd.to_datetime(df['Date & time (UTC)'], errors='coerce')
        
            dfs.append(df)  # Add the DataFrame to the list
        
    # # Concatenate all DataFrames in the list into one
    man_obs = pd.concat(dfs, ignore_index=True)
    man_obs = man_obs.replace('TR', 'Trois')
    # Optional: If you want to reset the index
    man_obs.reset_index(drop=True, inplace=True)
    man_obs['site_for_title'] = man_obs['site'] 
    
    
    replacements = {
        "UQAM-Sorel": "Sorel",
        "McGill-Gault": "Gault",
        "DOW-US-Plattsburgh": "Albany",
        "DOW-US-N": "Albany",
        "Albany-ESSX": "Albany",
        "CU-JEAN": "JEAN",
        "DOW-CAN-SE": "CU_DOW",
        "DOW-CAN-S":  "CU_DOW",
        "DOW-CAN-N":  "CU_DOW",
        "UQAM-Trois-Rivieres": "Trois"
    }

    man_obs["site"] = man_obs["site_for_title"].replace(replacements)
    
    
    man_obs = man_obs.replace('IP', 'PL')
    man_obs = man_obs.replace('SG', 'SN')
    man_obs = man_obs.applymap(lambda x: None if pd.isna(x) else x)
    man_obs = man_obs.replace('FZDZ', 'FZRA')
    man_obs = man_obs.replace('FRDZ', 'FZRA')
    man_obs = man_obs.replace('DZ', 'RA')
    man_obs = man_obs.replace('GS', 'SN')
    man_obs = man_obs.replace('Graupel', 'SN')
    man_obs = man_obs.replace('IC', None)
    man_obs = man_obs.replace('RA or FZRA', 'FZRA')
    man_obs = man_obs.replace('No ob',None)
    man_obs = man_obs.replace('IC/FZRA','FZRA')
    man_obs = man_obs.replace('FZRA/IC','FZRA')
    man_obs = man_obs.replace('Unknown or Uncertain',None) 
    man_obs = man_obs.replace('Unknown',None) 
    man_obs = man_obs.replace('unknown or uncertain',None) 
    man_obs = man_obs.replace('SN, FZDZ','SN') 
    man_obs = man_obs.replace('PL,FZDZ','PL') 
    
    # Correct dates
    man_obs['Date & time (UTC)'] =\
    man_obs['Date & time (UTC)'].apply(
        lambda d: d.replace(year=2022)
    )    
            
    
    files = os.listdir('soundings_wintremix/')#.remove('.DS_Store')
    files = [x for x in files if x != '.DS_Store']
    # files = files.remove('.DS_Store')
    site_list = [x.split('.')[3] for x  in files] 
    dates     = [ datetime.strptime(x.split('.')[2], "%Y%m%d%H%M") for x  in files] 

    df_all_events = pd.DataFrame({'site_complete_name':site_list,'file':files,'dates':dates})

    replacements = {
        "UQAM-Sorel": "Sorel",
        "McGill-Gault": "Gault",
        "Albany_DOW-US_Plattsburgh": "Albany",
        "Albany_DOW-US_N": "Albany",
        "Albany-ESSX": "Albany",
        "CU_JEAN": "JEAN",
        "CU_DOW-CAN_SE": "CU_DOW",
        "CU_DOW-CAN_S":  "CU_DOW",
        "CU_DOW-CAN_N":  "CU_DOW",
        "UQAM-Trois-Rivieres": "Trois"
    }

    df_all_events["site"] = df_all_events["site_complete_name"].replace(replacements, regex=False)

    # path_new = 'output_filename_Tw_summary3'
    # file = 'output_filename_Tw_summary.xlsx'
    # df_all_events = pd.read_excel(file)
    df_all_events = df_all_events.sort_values(by=['site','dates'])
    df_all_events = df_all_events[~((df_all_events['dates'] == '2022-03-06 14:00:00') & (df_all_events['site'] == 'Sorel'))]
    df_all_events = df_all_events[~((df_all_events['dates'] == '2022-03-06 14:00:00') & (df_all_events['site'] == 'Gault'))]
    df_all_events = df_all_events[~((df_all_events['dates'] == '2022-03-06 14:15:00') & (df_all_events['site'] == 'Sorel'))]
    df_all_events = df_all_events[~((df_all_events['dates'] == '2022-02-23 02:00:00') & (df_all_events['site'] == 'Trois'))]
    df_all_events = df_all_events[~((df_all_events['dates'] == '2022-02-18 05:00:00') & (df_all_events['site'] == 'Albany'))]

    # The soundings takes around 15 minutes to reach 4000 m. Hence we add 10 minutes to the launch to match with man obs.
    df_all_events['rounded_time'] = (df_all_events['dates'] + pd.Timedelta('10min')).dt.round('10min')
    df_all_events['site_date']    = df_all_events['site'].astype(str) + df_all_events['rounded_time'].astype(str)

    df_all_events_temp = df_all_events.copy(deep=True)
        
    '''Read soundings to get surface temperature'''
    for index, row in df_all_events.iterrows():
        # print(df_all_events.iloc[i].file.split('/')[-1])
        print(row.file)

        file = 'soundings_wintremix/'+ row.file.split('/')[-1]        
        
        if "Sorel" in file:
            df = soundings_functions.read_soundings_data_Sorel(file)
            df['site'] = 'Sorel'
            # df['date'] = date_datetime
            df['file'] = file
                
        elif "Gault" in file:
            df = soundings_functions.read_soundings_data_Gault(file)
            df['site'] = 'Gault'
            # df['date'] = date_datetime
            df['file'] = file
        elif "Trois" in file : 
            df = soundings_functions.read_soundings_data_TR(file)
            df['site'] = 'Trois'
            # df['date'] = date_datetime
            df['file'] = file
        elif "Albany" in file:
            df = soundings_functions.read_soundings_data_Albany(file)
            df['site'] = 'Albany'
            # df['date'] = date_datetime
            df['file'] = file
        elif "CU" in file:
            df = soundings_functions.read_soundings_data_DOW(file)
            df['site'] = 'CU'
            df['file'] = file
        
        df_all_events.at[index,'first_temp'] = df.Temp[1]
        df_all_events_temp.at[index,'first_temp'] = df.Temp[1]

        df['RelHumIce'] = df.apply(lambda row: 100.*kin1d_functions.rh_to_ice(row['RelHum']/100., row['Temp']+273.15,nanifneg=False), axis=1)
        
       
        
        output_sounding_analysis = stats_on_soundings(df,use_Tw=True)
        for output in output_sounding_analysis:
            df_all_events.loc[index, output] = output_sounding_analysis[output]
            
        output_sounding_analysis_temp = stats_on_soundings(df,use_Tw=False)
        for output in output_sounding_analysis_temp:
            df_all_events_temp.loc[index, output] = output_sounding_analysis_temp[output]

    man_obs_with_two_files_Tw   = man_obs.copy(deep=True) 
    man_obs_with_two_files_Tw   = add_interp_sounding_stats_to_manobs(man_obs_with_two_files_Tw,df_all_events,use_Tw=True)
    man_obs_with_two_files_Temp = man_obs.copy(deep=True) 
    man_obs_with_two_files_Temp = add_interp_sounding_stats_to_manobs(man_obs_with_two_files_Temp,df_all_events,use_Tw=False)
        
    def add_manobs_with_mix(row):
        if ((row['Secondary p-type']=='None') | (row['Secondary p-type']==None) | (row['Secondary p-type']==row['Primary p-type'])):
            row['manobs_with_mix'] = row['Primary p-type']
        elif  ((row['Primary p-type']=='PL') and (row['Secondary p-type']=='FZRA')\
          or (row['Primary p-type']=='FZRA') and (row['Secondary p-type']=='PL' )):
            row['manobs_with_mix'] = 'FZRA and PL'
        elif  ((row['Primary p-type']=='SN') and (row['Secondary p-type']=='PL')\
          or (row['Primary p-type']=='PL') and (row['Secondary p-type']=='SN' )):
            row['manobs_with_mix'] = 'SN and PL'
        elif  ((row['Primary p-type']=='RA') and (row['Secondary p-type']=='PL')\
          or (row['Primary p-type']=='PL') and (row['Secondary p-type']=='RA' )):
            row['manobs_with_mix'] = 'RA and PL'
        elif  ((row['Primary p-type']=='RA') and (row['Secondary p-type']=='FZRA')\
          or (row['Primary p-type']=='FZRA') and (row['Secondary p-type']=='RA' )
          or (row['Primary p-type']=='FZRA') and (row['Secondary p-type']=='FZRA' )):
            row['manobs_with_mix'] = 'FZRA'     
            row['Primary p-type'] = 'FZRA'
            row['Secondary p-type'] = None
        elif  ((row['Primary p-type']=='RA') and (row['Secondary p-type']=='SN')\
          or (row['Primary p-type']=='SN') and (row['Secondary p-type']=='RA' )):
            row['manobs_with_mix'] = 'SN and RA'
        elif  ((row['Primary p-type']=='FZRA') and (row['Secondary p-type']=='SN')\
          or (row['Primary p-type']=='SN') and (row['Secondary p-type']=='FZRA' )):
            row['manobs_with_mix'] = 'SN and FZRA'
        elif  ((row['Primary p-type']=='RA') and (row['Secondary p-type']=='RA')):
            row['manobs_with_mix'] = 'RA'
        elif  ((row['Primary p-type']=='PL') and (row['Secondary p-type']=='PL')):
            row['manobs_with_mix'] = 'PL'
        return row
    
    man_obs_with_two_files_Tw = man_obs_with_two_files_Tw.apply(add_manobs_with_mix,axis=1)
    man_obs_with_two_files_Temp = man_obs_with_two_files_Temp.apply(add_manobs_with_mix,axis=1)

    man_obs_with_two_files_Tw['rounded_time'] = (man_obs_with_two_files_Tw['Date & time (UTC)']-timedelta(minutes=2)).dt.round('10min')
    man_obs_with_two_files_Tw = man_obs_with_two_files_Tw.drop_duplicates(subset=['rounded_time','site_for_title'], keep='first')

    man_obs_with_two_files_Temp['rounded_time'] = (man_obs_with_two_files_Temp['Date & time (UTC)']-timedelta(minutes=2)).dt.round('10min')
    man_obs_with_two_files_Temp = man_obs_with_two_files_Temp.drop_duplicates(subset=['rounded_time','site_for_title'], keep='first')

    man_obs_with_two_files_Tw['site_date'] = man_obs_with_two_files_Tw['site'].astype(str) + man_obs_with_two_files_Tw['rounded_time'].astype(str)
    man_obs_with_two_files_Temp['site_date'] = man_obs_with_two_files_Temp['site'].astype(str) + man_obs_with_two_files_Temp['rounded_time'].astype(str)

    
    # Make sure that observations concurrent with soundings are kept. This avoids losing obs that are not between two soundings.
    vars_to_copy = ['wl_depth_hght', 'cl_depth_hght', 'cl_base_hght', 'max_WL_temp', 'min_CL_temp', 'CL_min_hght', 'WL_max_hght', 'WL_depth_1', 'WL_depth_2',\
                            'WL_mean_temp', 'CL_mean_temp', 'bot_cl_hght', 'top_cl_hght', 'bot_WL_hght', 'top_WL_hght', 'WL_integral', 'CL_integral','file',\
                                'CL_max_rhi','CL_max_rh']
    man_obs_with_two_files_Tw['Sounding'] = False    
    for _, row in df_all_events.iterrows():
        mask = man_obs_with_two_files_Tw.site_date == row.site_date
        if mask.any():
            for varname in vars_to_copy:
                man_obs_with_two_files_Tw.loc[man_obs_with_two_files_Tw.site_date==row.site_date, varname] = row[varname]    
            man_obs_with_two_files_Tw.loc[man_obs_with_two_files_Tw.site_date==row.site_date, 'Sounding'] = True     
        else:
            # create new row
            new_row = {col: pd.NA for col in man_obs_with_two_files_Tw.columns}
            new_row['site_date'] = row.site_date
            new_row['Sounding'] = True
            new_row['rounded_time'] = row['rounded_time']
            for varname in vars_to_copy:
                new_row[varname] = row[varname]

            man_obs_with_two_files_Tw = pd.concat([man_obs_with_two_files_Tw, pd.DataFrame([new_row])],ignore_index=True)
    
    man_obs_with_two_files_Temp['Sounding'] = False    
    for _, row in df_all_events_temp.iterrows():
        mask = man_obs_with_two_files_Temp.site_date == row.site_date
        if mask.any():
            for varname in vars_to_copy:
                man_obs_with_two_files_Temp.loc[man_obs_with_two_files_Temp.site_date==row.site_date, varname] = row[varname]    
            man_obs_with_two_files_Temp.loc[man_obs_with_two_files_Temp.site_date==row.site_date, 'Sounding'] = True     
        else:
            # create new row
            new_row = {col: pd.NA for col in man_obs_with_two_files_Temp.columns}
            new_row['site_date'] = row.site_date
            new_row['Sounding'] = True
            new_row['rounded_time'] = row['rounded_time']
            for varname in vars_to_copy:
                new_row[varname] = row[varname]

            man_obs_with_two_files_Temp = pd.concat([man_obs_with_two_files_Temp, pd.DataFrame([new_row])],ignore_index=True)

    man_obs_with_two_files_Tw   = man_obs_with_two_files_Tw.sort_values(['site_date','rounded_time'])
    man_obs_with_two_files_Temp = man_obs_with_two_files_Temp.sort_values(['site_date','rounded_time'])

    return man_obs_with_two_files_Tw,man_obs_with_two_files_Temp,df_all_events,df_all_events_temp


def read_soundings_data(file):
    
    if "Sorel" in file:
        df = soundings_functions.read_soundings_data_Sorel(file)
        df['site'] = 'Sorel'
        # df['date'] = date_datetime
        df['file'] = file
            
    elif "Gault" in file:
        df = soundings_functions.read_soundings_data_Gault(file)
        df['site'] = 'Gault'
        # df['date'] = date_datetime
        df['file'] = file
    elif "Trois" in file : 
        df = soundings_functions.read_soundings_data_TR(file)
        df['site'] = 'Trois'
        # df['date'] = date_datetime
        df['file'] = file
    elif "Albany" in file:
        df = soundings_functions.read_soundings_data_Albany(file)
        df['site'] = 'Albany'
        # df['date'] = date_datetime
        df['file'] = file
    elif "CU" in file:
        df = soundings_functions.read_soundings_data_DOW(file)
        df['site'] = 'CU'
        df['file'] = file
    return df


def plot_soundings_row(row,df_all_events):
    fig,ax = plt.subplots(1,figsize=(6,6))
    if (len(row.file)==2):
        file1 = 'soundings_wintremix/'+ row.file[0].split('/')[-1]
        file2 = 'soundings_wintremix/'+ row.file[1].split('/')[-1]
        two_files_bool=True
        file = file1
    else:
        file = 'soundings_wintremix/'+ row.file.split('/')[-1]
        two_files_bool=False
    
    if "Sorel" in file:
        if two_files_bool:
            df1 = soundings_functions.read_soundings_data_Sorel(file1)
            df2 = soundings_functions.read_soundings_data_Sorel(file2)
            df  = interpolate_between_soundings(row,df1,df2,df_all_events)
        else:
            df = soundings_functions.read_soundings_data_Sorel(file)
        df['site'] = 'Sorel'
        # df['date'] = date_datetime
        df['file'] = file
            
    elif "Gault" in file:
        if two_files_bool:
            df1 = soundings_functions.read_soundings_data_Gault(file1)
            df2 = soundings_functions.read_soundings_data_Gault(file2)
            df  = interpolate_between_soundings(row,df1,df2,df_all_events)
        else:
            df = soundings_functions.read_soundings_data_Gault(file)
        df['site'] = 'Gault'
        # df['date'] = date_datetime
        df['file'] = file
    elif "Trois" in file : 
        if two_files_bool:
            df1 = soundings_functions.read_soundings_data_TR(file1)
            df2 = soundings_functions.read_soundings_data_TR(file2)
            df  = interpolate_between_soundings(row,df1,df2,df_all_events)
        else:
            df = soundings_functions.read_soundings_data_TR(file)
        df['site'] = 'Trois'
        # df['date'] = date_datetime
        df['file'] = file
    elif "Albany" in file:
        if two_files_bool:
            df1 = soundings_functions.read_soundings_data_Albany(file1)
            df2 = soundings_functions.read_soundings_data_Albany(file2)
            df  = interpolate_between_soundings(row,df1,df2,df_all_events)
        else:
            df = soundings_functions.read_soundings_data_Albany(file)
        df['site'] = 'Albany'
        # df['date'] = date_datetime
        df['file'] = file
    elif "CU" in file:
        df = soundings_functions.read_soundings_data_DOW(file)
        df['site'] = 'CU'
        df['file'] = file


    # df_all_events.at[index,'first_temp'] = df.Temp[1]
    
    
    types_to_plot = {'RA':0, 'SN':1, 'FZRA':2, 'PL':3,'SN and PL':4,'SN and RA':5,'FZRA and PL':6,'RA and PL':7}


    # #Interpolate
    heights_composite = np.arange(50,5000,10)

    df = df[1:]
    f_int_t = interp1d(df.height, df.Temp,bounds_error=False, fill_value=np.nan)
    f_int_tw = interp1d(df.height, df.Tw, bounds_error=False, fill_value=np.nan)
    df_inter = pd.DataFrame({'height':heights_composite,'Temp':f_int_t(heights_composite),'Tw':f_int_tw(heights_composite)})

    ax.plot(df_inter.Tw,df_inter.height/1000.,color='C1',alpha=0.3)
    ax.plot(df_inter.Temp,df_inter.height/1000.,color='k',alpha=0.3)


    ax.set_xlim([-15,10])
    ax.set_ylim([0,4])

    ax.set_title(row['site_date']+ ' - '+ str(row['manobs_with_mix']))

    ax.set_xlabel('Tw [C]')
    ax.grid('on',linewidth=0.3)

    ax.set_ylabel('MSL [km]')
    ax.set_ylabel('MSL [km]')

    plt.show()
    
    return df



def add_pr_t_rh(df_manobs):
    
    
    nc_gault = xr.open_dataset('surface_measurements_wintremix/Gault/CFI_Sentinels_Gault_Parsivel_Disdrometer_data_WINTRE-MIX_03.nc')
    
    
    df_all_events_pr = copy.copy(df_manobs)
    df_all_events_pr['rounded_time'] = pd.to_datetime(df_all_events_pr['rounded_time'])
    df_all_events_pr['Date & time (UTC)'] = pd.to_datetime(df_manobs['Date & time (UTC)'])

    # IOP4 Plattsburgh
    # conds = (df_all_events_pr.site=='Albany') & (df_all_events_pr.IOP=='IOP4')
    # df_for_plot = copy.copy(df_all_events_pr[conds])
    # start,end = pd.to_datetime(df_for_plot['rounded_time'].values[0]),pd.to_datetime(df_for_plot['rounded_time'].values[-1])
    # parsi_to_merge = nc_chazy_iop4['Prcp_Intensity'].sel(time=slice(start,end)).resample(time='10min').mean().to_dataframe().rename(columns={'Prcp_Intensity':'pr'})
    # subset = df_all_events_pr.loc[conds].copy()
    # subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    # merged = pd.merge(subset['rounded_time'], parsi_to_merge['pr'], left_on='rounded_time', right_on='time', how='left')
    # df_all_events_pr.loc[conds, merged.columns] = merged.values
    # print('IOP4 Plattsburgh')
    # print(df_all_events_pr.columns)

    # IOP4 DOW-CAN-S
    # conds=(df_all_events_pr.site=='CU_DOW') & (df_all_events_pr.IOP=='IOP4')
    # df_for_plot = copy.copy(df_all_events_pr[conds])
    # start,end = pd.to_datetime(df_for_plot['rounded_time'].values[0]),pd.to_datetime(df_for_plot['rounded_time'].values[-1])
    # parsi_to_merge = nc_chazy_iop4['Prcp_Intensity'].sel(time=slice(start,end)).resample(time='10min').mean().to_dataframe().rename(columns={'Prcp_Intensity':'pr'})
    # subset = df_all_events_pr.loc[conds].copy()
    # subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    # merged = pd.merge(subset['rounded_time'], parsi_to_merge['pr'], left_on='rounded_time', right_on='time', how='left')
    # df_all_events_pr.loc[conds, 'pr'] = merged.pr.values
    # print('IOP4 DOW-CAN-S')
    # print(df_all_events_pr.columns)

    # IOP4 Gault
    conds = (df_manobs.site=='Gault') & (df_manobs.IOP=='IOP4')
    df_for_plot = copy.copy(df_manobs[conds])
    start,end = pd.to_datetime(df_for_plot['rounded_time'].values[0]),pd.to_datetime(df_for_plot['rounded_time'].values[-1])
    parsi_to_merge = nc_gault['Prcp_Intensity'].sel(time=slice(start,end)).to_dataframe().shift(freq='5min').resample('10min').mean().rename(columns={'Prcp_Intensity':'pr'})
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], parsi_to_merge['pr'], left_on='rounded_time', right_index=True, how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values
    print('IOP4 Gault')

    # IOP5 Gault
    conds = (df_manobs.site=='Gault') & (df_manobs.IOP=='IOP5')
    df_for_plot = copy.copy(df_manobs[conds])
    start,end = pd.to_datetime(df_for_plot['rounded_time'].values[0]),pd.to_datetime(df_for_plot['rounded_time'].values[-1])
    parsi_to_merge = nc_gault['Prcp_Intensity'].sel(time=slice(start,end)).to_dataframe().shift(freq='5min').resample('10min').mean().rename(columns={'Prcp_Intensity':'pr'})
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], parsi_to_merge['pr'], left_on='rounded_time', right_on='time', how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values
    print('IOP5 Gault')

    # IOP5 Sorel
    conds = (df_manobs.site=='Sorel') & (df_manobs.IOP=='IOP5')
    df_for_plot = copy.copy(df_manobs[conds])
    start,end = pd.to_datetime(df_for_plot['Date & time (UTC)'].values[0]),pd.to_datetime(df_for_plot['Date & time (UTC)'].values[-1])
    output_sorel = MetData_functions.import_parsi(start,end)[0].shift(freq='5min').resample('10min').mean().rename(columns={'intensite':'pr'})
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], output_sorel['pr'], left_on='rounded_time', right_index=True, how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values
    print('IOP5 Sorel')

    # IOP5 Trois
    conds=(df_manobs.site=='Trois') & (df_manobs.IOP=='IOP5')
    df_for_plot = copy.copy(df_manobs[conds])
    start,end = pd.to_datetime(df_for_plot['Date & time (UTC)'].values[0]),pd.to_datetime(df_for_plot['Date & time (UTC)'].values[-1])
    output_trois = parsivel_functions.parsi_import_plist_file(start,end)[0].shift(freq='5min').resample('10min').mean().rename(columns={'intensite':'pr'})
    output_trois.columns = ['pr']
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], output_trois['pr'], left_on='rounded_time', right_index=True, how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values
    print('IOP5 Trois')


    # IOP8 Sorel
    conds = (df_manobs.site=='Sorel') & (df_manobs.IOP=='IOP8')
    df_for_plot = copy.copy(df_manobs[conds])
    start,end = pd.to_datetime(df_for_plot['Date & time (UTC)'].values[0]),pd.to_datetime(df_for_plot['Date & time (UTC)'].values[-1])
    output_sorel = MetData_functions.import_parsi(start,end)[0].shift(freq='5min').resample('10min').mean().rename(columns={'intensite':'pr'})
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], output_sorel['pr'], left_on='rounded_time', right_index=True, how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values
    print('IOP8 Sorel')

    # # IOP8 Trois
    conds=(df_manobs.site=='Trois') & (df_manobs.IOP=='IOP8')
    df_for_plot = copy.copy(df_manobs[conds])
    start,end = pd.to_datetime(df_for_plot['Date & time (UTC)'].values[0]),pd.to_datetime(df_for_plot['Date & time (UTC)'].values[-1])
    output_trois = parsivel_functions.parsi_import_plist_file(start,end)[0].shift(freq='5min').resample('10min').mean()
    output_trois.columns = ['pr']
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], output_trois['pr'], left_on='rounded_time', right_index=True, how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values
    print('IOP8 Trois')

    # # IOP9 Gault
    conds = (df_manobs.site=='Gault') & (df_manobs.IOP=='IOP9')
    df_for_plot = copy.copy(df_manobs[conds])
    start,end = df_for_plot['Date & time (UTC)'].values[0],df_for_plot['Date & time (UTC)'].values[-1]
    nc_gault_to_plot = nc_gault['Prcp_Intensity'][178828:178828+7*60+50].to_dataframe().shift(freq='5min').resample('10min').mean().rename(columns={'Prcp_Intensity':'pr'})
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_gault_to_plot['pr'], left_on='rounded_time', right_on='time', how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values
    print('IOP9 Gault')


    '''Add T and RH'''
    nc_sentinels = xr.open_dataset('surface_measurements_wintremix/sentinels_metdata.nc')

    # IOP4 Gault
    conds = (df_all_events_pr.site=='Gault') & (df_all_events_pr.IOP=='IOP4')
    df_for_plot = copy.copy(df_all_events_pr[conds]) 
    start,end = df_for_plot['Date & time (UTC)'].values.min(),df_for_plot['Date & time (UTC)'].values.max()
    nc_sentinels_to_plot = nc_sentinels.sel(station='GAUL',time=slice(start,end))[['temp_2m','relative_humidity']].to_dataframe()[['temp_2m','relative_humidity']].shift(freq='5min').resample('10min').mean()
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_sentinels_to_plot, left_on='rounded_time', right_on='time', how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values

    # # IOP5 Gault
    conds = (df_all_events_pr.site=='Gault') & (df_all_events_pr.IOP=='IOP5')
    df_for_plot = copy.copy(df_all_events_pr[conds])
    start,end = df_for_plot['Date & time (UTC)'].values.min(),df_for_plot['Date & time (UTC)'].values.max()
    nc_sentinels_to_plot = nc_sentinels.sel(station='GAUL',time=slice(start,end))[['temp_2m','relative_humidity']].to_dataframe()[['temp_2m','relative_humidity']].shift(freq='5min').resample('10min').mean()
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_sentinels_to_plot, left_on='rounded_time', right_on='time', how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values

    # # IOP5 Sorel
    conds= (df_all_events_pr.site=='Sorel') & (df_all_events_pr.IOP=='IOP5')
    df_for_plot = copy.copy(df_all_events_pr[conds])
    start,end = df_for_plot['Date & time (UTC)'].values.min(),df_for_plot['Date & time (UTC)'].values.max()
    output_sorel = MetData_functions.import_wxt(pd.to_datetime(start),pd.to_datetime(end))
    nc_sentinels_to_plot = output_sorel[['temperature','RH']].shift(freq='5min').resample('10min').mean().rename(columns={'temperature':'temp_2m','RH':'relative_humidity'})
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_sentinels_to_plot, left_on='rounded_time', right_index=True, how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values

    #IOP5 Trois
    conds = (df_all_events_pr.site=='Trois') & (df_all_events_pr.IOP=='IOP5')
    df_for_plot = copy.copy(df_all_events_pr[conds])
    start,end = df_for_plot['Date & time (UTC)'].values.min(),df_for_plot['Date & time (UTC)'].values.max()
    nc_sentinels_to_plot = nc_sentinels.sel(station='TROI',time=slice(start,end))[['temp_2m','relative_humidity']].to_dataframe()[['temp_2m','relative_humidity']].shift(freq='5min').resample('10min').mean()
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_sentinels_to_plot, left_on='rounded_time', right_on='time', how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values

    # IOP8 Sorel
    conds = (df_all_events_pr.site=='Sorel') & (df_all_events_pr.IOP=='IOP8')
    df_for_plot = copy.copy(df_all_events_pr[conds])
    start,end = df_for_plot['Date & time (UTC)'].values.min(),df_for_plot['Date & time (UTC)'].values.max()
    output_sorel = MetData_functions.import_wxt(pd.to_datetime(start),pd.to_datetime(end))
    nc_sentinels_to_plot = output_sorel[['temperature','RH']].shift(freq='5min').resample('10min').mean().rename(columns={'temperature':'temp_2m','RH':'relative_humidity'})
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_sentinels_to_plot, left_on='rounded_time', right_index=True, how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values

    # IOP8 Trois
    conds = (df_all_events_pr.site=='Trois') & (df_all_events_pr.IOP=='IOP8')
    df_for_plot = copy.copy(df_all_events_pr[conds])
    start,end = df_for_plot['Date & time (UTC)'].values.min(),df_for_plot['Date & time (UTC)'].values.max()
    nc_sentinels_to_plot = nc_sentinels.sel(station='TROI',time=slice(start,end))[['temp_2m','relative_humidity']].to_dataframe()[['temp_2m','relative_humidity']].shift(freq='5min').resample('10min').mean()
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_sentinels_to_plot, left_on='rounded_time', right_on='time', how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values

    # IOP9 Gault
    conds =  (df_all_events_pr.site=='Gault') & (df_all_events_pr.IOP=='IOP9')
    df_for_plot = copy.copy(df_all_events_pr[conds])
    start,end = df_for_plot['Date & time (UTC)'].values.min(),df_for_plot['Date & time (UTC)'].values.max()
    nc_sentinels_to_plot = nc_sentinels.sel(station='GAUL',time=slice(start,end))[['temp_2m','relative_humidity']].to_dataframe()[['temp_2m','relative_humidity']].shift(freq='5min').resample('10min').mean()
    subset = df_all_events_pr.loc[conds].copy()
    subset['rounded_time'] = pd.to_datetime(subset['rounded_time'])
    merged = pd.merge(subset['rounded_time'], nc_sentinels_to_plot, left_on='rounded_time', right_on='time', how='left')
    df_all_events_pr.loc[conds, merged.columns] = merged.values



    return df_all_events_pr

