import pandas as pd
import numpy as np 
import glob
import matplotlib.pyplot as plt
# from metpy.plots import Hodograph, SkewT
# import metpy.calc as mpcalc
# from metpy.units import units, concatenate
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

def wet_bulb_tw(T, RH):
    a = T * np.arctan(0.151977 * (RH + 8.313659) ** 0.5)
    b =  np.arctan(T + RH) - np.arctan(RH - 1.676331)
    c = 0.00391838 * RH ** (1.5) *  np.arctan(0.023101 * RH)
    d = -4.686035
    Tw = a + b + c + d
    return Tw

def pot_t(p,T):
    T=T+273.15 #Conversion de T en K
    pot = T*(1000/p)**(2/7)
#     pot = pot - 273.15
    return pot

def read_soundings_data_Sorel(path_soundings):
    df = pd.read_table(path_soundings,sep=r'\s+',
                       skiprows=(1, 2), 
                       encoding='unicode_escape')

    if df.iloc[0].isna().any():
        df = df.iloc[1:].reset_index(drop=True)

    elif df['UTC_Time'].isin(['AM', 'PM']).any():
        columns_to_keep = list(df.columns[:7]) + [df.columns[-1]]
        df = df[columns_to_keep]

        # Ajouter un index
        df = df.set_index('FltTime')
        df = df.reset_index(drop=False)
        df = df.rename_axis('Flt Time')

        # Changer le nom des colonnes
        nouveaux_noms_colonnes = ['Press', 'Temp', 'RelHum', 'WSpeed', 'WDirn', 'UTC_Date', 'UTC_Time', 'GPM_MSL']
        df.columns = nouveaux_noms_colonnes

        # Calcul de Tw
    df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
    df.rename(columns={'GPM_MSL': 'height'}, inplace=True)
    return df

def read_soundings_data_Gault(path_soundings):
    df = pd.read_table(path_soundings, sep=r'\s+', skiprows=(1, 2), encoding='unicode_escape')
    df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
    df.rename(columns={'GPM_AGL': 'height'}, inplace=True)
    return df
    
def read_soundings_data_DOW(path_soundings):       
    df = pd.read_table(path_soundings, sep=',', skiprows=(1, 2), encoding='unicode_escape')
    if df.columns[0] == 'ï»¿/Row/@Altitude':
        # Sélectionner les bonnes colonnes
        columns_to_keep = ['ï»¿/Row/@Altitude', '/Row/@Pressure', '/Row/@Temperature', '/Row/@Humidity']
        df = df[columns_to_keep]

        # Renommer les colonnes
        nouveaux_noms_colonnes = {
            'ï»¿/Row/@Altitude': 'GPM_MSL',
            '/Row/@Pressure': 'Press',
            '/Row/@Temperature': 'Temp',
            '/Row/@Humidity': 'RelHum'
        }

        # Utiliser la méthode rename pour renommer les colonnes
        df = df.rename(columns=nouveaux_noms_colonnes)

        # Convertir Temp de K en C
        df.Temp = df.Temp - 273.15
        df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
        df.rename(columns={'GPM_MSL': 'height'}, inplace=True)

    elif df.columns[0] == '/Row/@Altitude':

        # Sélectionner les bonnes colonnes
        columns_to_keep = ['/Row/@Altitude', '/Row/@Pressure', '/Row/@Temperature', '/Row/@Humidity']
        df = df[columns_to_keep]

        # Renommer les colonnes
        nouveaux_noms_colonnes = {
            '/Row/@Altitude': 'GPM_MSL',
            '/Row/@Pressure': 'Press',
            '/Row/@Temperature': 'Temp',
            '/Row/@Humidity': 'RelHum'
        }

        # Utiliser la méthode rename pour renommer les colonnes
        df = df.rename(columns=nouveaux_noms_colonnes)

        # Convertir Temp de K en C
        df.Temp = df.Temp - 273.15
        df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
        df.rename(columns={'GPM_MSL': 'height'}, inplace=True)

    elif df.columns[0] == 'ï»¿/Row/@Altitude':
        # Sélectionner les bonnes colonnes
        columns_to_keep = ['/Row/@Altitude', '/Row/@Pressure', '/Row/@Temperature', '/Row/@Humidity']
        df = df[columns_to_keep]

        # Renommer les colonnes
        nouveaux_noms_colonnes = {
            '/Row/@Altitude': 'GPM_MSL',
            '/Row/@Pressure': 'Press',
            '/Row/@Temperature': 'Temp',
            '/Row/@Humidity': 'RelHum'
        }

        # Utiliser la méthode rename pour renommer les colonnes
        df = df.rename(columns=nouveaux_noms_colonnes)

        # Convertir Temp de K en C
        df.Temp = df.Temp - 273.15
        df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
        df.rename(columns={'GPM_MSL': 'height'}, inplace=True)

    elif df.columns[0] == 'ï»¿/Row/@Height':
        # Sélectionner les bonnes colonnes
        columns_to_keep = ['ï»¿/Row/@Height', '/Row/@Pressure', '/Row/@Temperature', '/Row/@Humidity']
        df = df[columns_to_keep]

        # Renommer les colonnes
        nouveaux_noms_colonnes = {
            'ï»¿/Row/@Height': 'GPM_MSL',
            '/Row/@Pressure': 'Press',
            '/Row/@Temperature': 'Temp',
            '/Row/@Humidity': 'RelHum'
        }

        # Utiliser la méthode rename pour renommer les colonnes
        df = df.rename(columns=nouveaux_noms_colonnes)

        # Convertir Temp de K en C
        df.Temp = df.Temp - 273.15
        df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
        df.rename(columns={'GPM_MSL': 'height'}, inplace=True)
    else:
        columns_to_keep = ['Temp','GPM_MSL','RelHum']
        df = df[columns_to_keep]
        #df.Temp = df.Temp - 273.15
        df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
        df.rename(columns={'GPM_MSL': 'height'}, inplace=True)
    
    return df

def read_soundings_data_TR(path_soundings):
    df = pd.read_table(path_soundings, 
                       sep=r'\s+', 
                       skiprows=(1, 2), 
                       encoding='unicode_escape')
    df.reset_index(inplace=True)

    if df['UTC_Time'].isin(['AM', 'PM']).any():
        columns_to_keep = ['Press', 'Temp', 'RelHum', 'GPM_MSL']
        df = df[columns_to_keep]
        df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
        df.rename(columns={'GPM_MSL': 'height'}, inplace=True)

    else:
        df = pd.DataFrame(columns=['Press', 'Temp', 'RelHum', 'GPM_MSL', 'Tw'])  # Return an empty DataFrame with columns
        df.rename(columns={'GPM_MSL': 'height'}, inplace=True)

    return df

def read_soundings_data_Albany(path_soundings):
    # df=pd.read_table(path_soundings, sep=r'\s+', 
    #              skiprows=(1,2),
    #              encoding= 'unicode_escape')

    df = pd.read_csv(path_soundings,
                    sep='\s+',
                    skiprows=5,  # Ignorer les 6 premières lignes
                    encoding='unicode_escape',
                    names=[
                        "UTC_Date", "UTC_Time", "FltTime", "Ascent", "GPM_AGL", "GPM_MSL",
                        "Alt_AGL", "Alt_MSL", "Press", "Temp", "RelHum", "Mix_Rat", "DP",
                        "WSpeed", "WDirn", "Long/E", "Lat/N"
                        ]  # Définir manuellement les noms de colonnes
                )

    df.reset_index(inplace=True)
    
    if df['UTC_Time'].isin(['AM', 'PM']).any():        
        columns_to_keep = [ 'Press' ,'Temp','RelHum','GPM_MSL']
        df = df[columns_to_keep]
        
    # Calcul de Tw
    df['Tw'] = wet_bulb_tw(df['Temp'], df['RelHum'])
    df.rename(columns={'GPM_MSL': 'height'}, inplace=True)
    return df


def read_soundings_data_JEAN(path_soundings):    
    df = pd.read_csv(path_soundings, sep=',', skiprows=(1, 2), encoding='unicode_escape')

    if df.columns[0] == 'ï»¿/Row/@Altitude':
        # Sélectionner les bonnes colonnes
        columns_to_keep = ['ï»¿/Row/@Altitude', '/Row/@Pressure', '/Row/@Temperature', '/Row/@Humidity']

        df = df[columns_to_keep]

        # Renommer les colonnes
        nouveaux_noms_colonnes = {
            'ï»¿/Row/@Altitude': 'height',
            '/Row/@Pressure': 'Press',
            '/Row/@Temperature': 'Temp',
            '/Row/@Humidity': 'RelHum'
        }

        # Utiliser la méthode rename pour renommer les colonnes
        df = df.rename(columns=nouveaux_noms_colonnes)

        # Convertir Temp de K en C
        df['Temp'] = df['Temp'] - 273.15
        
    elif df.columns[0] == '/Row/@Altitude':
        # Sélectionner les bonnes colonnes
        columns_to_keep = ['/Row/@Altitude', '/Row/@Pressure', '/Row/@Temperature', '/Row/@Humidity']
        df = df[columns_to_keep]
        
        # Renommer les colonnes
        nouveaux_noms_colonnes = {
            '/Row/@Altitude': 'height',
            '/Row/@Pressure': 'Press',
            '/Row/@Temperature': 'Temp',
            '/Row/@Humidity': 'RelHum'
        }

        # Utiliser la méthode rename pour renommer les colonnes
        df = df.rename(columns=nouveaux_noms_colonnes)

        # Convertir Temp de K en C
        df['Temp'] = df['Temp'] - 273.15
        
    else : 
        columns_to_keep = ['GPM_MSL', 'Press', 'Temp', 'RelHum']
        df = df[columns_to_keep]  
            # Renommer les colonnes
        nouveaux_noms_colonnes = {
            'GPM_MSL': 'height'
        }

        # Utiliser la méthode rename pour renommer les colonnes
        df = df.rename(columns=nouveaux_noms_colonnes)
    
    return df 

def process_soundings_files(files):
    """
    Cette fonction lit les fichiers de ballons sondes, les traite et les stocke dans des listes de DataFrames,
    avec des informations supplémentaires sur chaque fichier (site, date, etc.).
    
    Arguments :
    files -- Liste des fichiers à traiter
    soundings -- Module contenant les fonctions de lecture des données spécifiques pour chaque site
    
    Retourne :
    dfs_ballons -- Liste des DataFrames traités
    dfs_ballons_Sorel, dfs_ballons_Albany, dfs_ballons_TROI, dfs_ballons_Gault, dfs_ballons_DOW, dfs_ballons_JEAN -- Listes spécifiques aux sites
    df_names -- Liste des noms de fichiers
    dates -- Liste des dates extraites des noms de fichiers
    dates_datetime -- Liste des dates au format datetime
    sites -- Liste des sites extraits des noms de fichiers
    treated_files -- Liste des fichiers traités
    """
    
    # Initialiser les listes pour chaque site et les informations globales
    dfs_ballons_Sorel = []
    dfs_ballons_Albany = []
    dfs_ballons_TROI = []
    dfs_ballons_Gault = []
    dfs_ballons_DOW = []
    dfs_ballons_JEAN = []
    dfs_ballons = [] 
    df_names = []
    dates = []
    dates_datetime = []
    sites = []
    treated_files = []

    for file in files:        
        # Extraire la date du fichier
        date = file.split('.')[2]
        dates.append(date)
        
        # Extraire la date au format datetime
        year = int(date[0:4])
        month = int(date[4:6])
        day = int(date[6:8])
        hour = int(date[8:10])
        minute = int(date[10:12])
        date_datetime = datetime(year, month, day, hour, minute, 0)
        dates_datetime.append(date_datetime)
        
        # Extraire le site à partir du nom de fichier
        site = file.split('.')[3]
        sites.append(site)
        
        # Ajouter à la liste des fichiers traités
        treated_files.append(file)
        
        # Processus en fonction du site
        if "Sorel" in file:
            df = read_soundings_data_Sorel(file)
            df['site'] = 'Sorel'
            df['date'] = date_datetime
            df['file'] = file
            dfs_ballons_Sorel.append(df)
            dfs_ballons.append(df)

        elif "Gault" in file:
            df = read_soundings_data_Gault(file)
            df['site'] = 'Gault'
            df['date'] = date_datetime
            df['file'] = file
            dfs_ballons_Gault.append(df)
            dfs_ballons.append(df)

        elif "Trois" in file:
            df = read_soundings_data_TR(file)
            df['site'] = 'Trois'
            df['date'] = date_datetime
            df['file'] = file
            df.rename(columns={"height": "Height"}, inplace=True)
            dfs_ballons_TROI.append(df)
            dfs_ballons.append(df)

        elif "Albany_DOW-US_Plattsburgh.txt" in file:
            df = read_soundings_data_Albany(file)
            df['site'] = 'Albany'
            df['date'] = date_datetime
            df['file'] = file
            df.rename(columns={"height": "Height"}, inplace=True)
            dfs_ballons_Albany.append(df)
            dfs_ballons.append(df)
        
        # Ajouter le nom du fichier à la liste
        df_names.append(file)

    # Retourner les résultats
    return (dfs_ballons, dfs_ballons_Sorel, dfs_ballons_Albany, dfs_ballons_TROI, 
            dfs_ballons_Gault, dfs_ballons_DOW, dfs_ballons_JEAN, 
            df_names, dates, dates_datetime, sites, treated_files)
