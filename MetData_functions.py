import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import re


'''
Parsivel
'''
def import_parsi(start,end):
    timelist = []
    mat= []
    intensite = []

    current=start
    end_to_import = datetime(end.year,end.month,end.day,0,0)
    while ((timelist==[]) or (timelist[-1]<end_to_import )):
        path = '/home/rml000/sitestore8/UQAM/projet_margaux/measurements_papier_margaux/Sorel/rubber_daqy/' + str(current.year)+'_'+str(current.month).zfill(2) + '/raw/'
        f = open(path+'Parsivel'+ str(current.year)[2:]+str(current.month).zfill(2)+str(current.day).zfill(2)+'.txt','r')
        for line in f:
            if ((line =='\n') or (line=='HEADER: item1, item2, item3\n')):
                yo=1
            else:    
                linesplit = line.split(',')

                bins1024=np.array(linesplit[-1025:-1]).astype(float)

                mat.append(bins1024)
                time = datetime.strptime(linesplit[0][0:19],'%Y/%m/%d %H:%M:%S')
                timelist.append(time)
                intensite.append(np.array(linesplit[0][21:]).astype(float) )
                if (timelist[-1]>=end):
                    break
        f.close()
        current=current+timedelta(days=1)


    df_parsivel = pd.DataFrame(index = timelist, data=np.array(mat))
    prec_tab = pd.DataFrame(index=timelist,data=np.array(intensite),columns=['intensite'])
    
    return prec_tab,df_parsivel

def time_series_size(start,end, appended_data, ax):

    filtred_data = appended_data[(appended_data.index >= start) & (appended_data.index <= end)]

    Dbins=[0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17,19,21.5,24.5]
    Dspread = [0.125 ,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.250, 0.250, 0.500, 0.500, 0.500, 0.500, 0.500, 1.000, 1.000, 1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 2.000, 2.000, 3.000, 3.000]
    Dspread2D = np.array([Dspread]*len(filtred_data.index))

    resample_data = np.sum(np.array(filtred_data).reshape(len(filtred_data),32,32),axis=-2)
    resample_data = np.ma.masked_where(0,resample_data)

    pc = ax.pcolormesh(filtred_data.index, Dbins, np.log10(resample_data/Dspread2D).T, cmap="rainbow",vmin=0,vmax=3)

    ax.set_ylim([0,8])
    ax.set_ylabel('D [mm]')

    yo2 = plt.colorbar(pc,ticks=[0,1,2,3],ax=ax, pad=0.02)
    yo2.ax.set_yticklabels([r'$10^0$',r'$10^1$',r'$10^2$',r'$10^3$'])
    yo2.set_label(r'Count [min$^{-1}$ mm$^{-1}$]')


def time_series_speed(start,end, appended_data, ax):

    filtred_data = appended_data[(appended_data.index >= start) & (appended_data.index <= end)]

    vTbins=[0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.10,1.30,1.5,1.7,1.9,2.2,2.6,3,3.4,3.8,4.4,5.2,6,6.8,7.6,8.8,10.4,12,13.6,15.2,17.6,20.80]
    Vspread = [0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.200, 0.200, 0.200, 0.200, 0.200,0.400, 0.400, 0.400, 0.400, 0.400, 0.800, 0.800, 0.800, 0.800, 0.800, 1.600, 1.600, 1.600, 1.600, 1.600, 3.200, 3.200]
    Vspread2D = np.array([Vspread]*len(filtred_data.index))

    resample_data = np.sum(filtred_data.values.reshape(len(filtred_data),32,32),axis=-1)
    resample_data = np.ma.masked_where(0,resample_data)

    pc = ax.pcolormesh(filtred_data.index, vTbins, np.log10(resample_data/Vspread2D).T, cmap="rainbow",vmin=0,vmax=3)

    ax.set_ylim([0,13])
    ax.set_ylabel('v [m/s]')

    yo2 = plt.colorbar(pc,ticks=[0,1,2,3],ax=ax, pad=0.02)
    yo2.ax.set_yticklabels([r'$10^0$',r'$10^1$',r'$10^2$',r'$10^3$'])
    yo2.set_label(r'Count [min$^{-1}$ (m/s)$^{-1}$]')

def parsi_spectrum_plot(start,end,appended_data):

    # filtre des donnees sur une periode specifique
    filtred_data = appended_data[(appended_data.index >= start) & (appended_data.index <= end)]

    Dbins=[0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17,19,21.5,24.5]
    vTbins=[0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.10,1.30,1.5,1.7,1.9,2.2,2.6,3,3.4,3.8,4.4,5.2,6,6.8,7.6,8.8,10.4,12,13.6,15.2,17.6,20.80]
    Dspread = [0.125 ,0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.250, 0.250,          0.500, 0.500, 0.500, 0.500, 0.500, 1.000, 1.000, 1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 2.000, 2.000, 3.000, 3.000]
    Vspread = [0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.200, 0.200, 0.200, 0.200, 0.200,0.400, 0.400, 0.400, 0.400, 0.400, 0.800, 0.800, 0.800, 0.800, 0.800, 1.600, 1.600, 1.600, 1.600, 1.600, 3.200, 3.200]

    Dspread2D = np.array([Dspread]*32)
    Vspread2D = np.array([Vspread]*32).T

    dgrid=np.zeros(len(Dbins)+1)
    vgrid=np.zeros(len(vTbins)+1)

    for i in range(len(dgrid)):
        if i==len(dgrid)-1:
            dgrid[i] = Dbins[i-1]+Dspread[i-1]/2
            vgrid[i] = vTbins[i-1]+Dspread[i-1]/2
        else:
            dgrid[i] = Dbins[i]-Dspread[i]/2
            vgrid[i] = vTbins[i]-Dspread[i]/2


    resample_data = filtred_data.sum(axis = 0)  #somme des particules au 10 minutes 
    resample_data = resample_data.replace(0, np.nan)
    b = np.array(resample_data).reshape(32,32)
    b = np.ma.masked_where(np.isnan(b),b)


    fig, ax = plt.subplots(1,1)
    X,Y = np.meshgrid(dgrid,vgrid)

    if np.max(b/Dspread2D/Vspread2D/((end-start).seconds/60))>100:
        pc = ax.pcolormesh(X,Y,np.log10(b/Dspread2D/Vspread2D/((end-start).seconds/60)) , cmap="rainbow", vmin=0, vmax=3)
        yo2 = fig.colorbar(pc,ticks=[0,1,2,3], fraction=0.046, pad=0.04)
        yo2.ax.set_yticklabels([r'$1$',r'$10$',r'$100$',r'$1000$'])
    else:
        pc = ax.pcolormesh(X,Y,np.log10(b/Dspread2D/Vspread2D/((end-start).seconds/60)) , cmap="rainbow", vmin=0, vmax=2)
        yo2 = fig.colorbar(pc,ticks=[0,1,2], fraction=0.046, pad=0.04)
        yo2.ax.set_yticklabels([r'$1$',r'$10$',r'$100$'])

    yo2.set_label('Raw counts, Normalized by\n Bin size [(m s$^{-1}$ mm)$^{-1}$]')

    for vT in vgrid:
        ax.axhline(y=vT, xmin=0., xmax=18, linewidth=0.5, color = 'gray')
    for D in dgrid:
        ax.axvline(x=D, ymin=0.0, ymax = 10, linewidth=0.5, color='gray')

    Dlist = np.linspace(0,0.006,100)
    ax.plot(Dlist*1000,4854*Dlist**1*np.exp(-195*Dlist),label='THOM Rain')
    ax.plot(Dlist*1000,442*Dlist**0.89,label='THOM Graupel')
    ax.plot(Dlist*1000,40*Dlist**0.55*np.exp(-100*Dlist),label='THOM snow')

    ax.set_title(str(start)[0:16]+'-'+str(end)[11:16] + ' UTC')
    ax.axis((0,Dbins[31],0,vTbins[30]))


    ax.set_xlabel('Diameters (mm)')
    ax.set_ylabel('Speed (m/s)')
    ax.set_xlim([0,5])
    ax.set_ylim([0,10])

    ax.legend(loc='upper left')

'''
Pluvio
'''

def import_pluvio(start,end):
    
    timelist = []
    mat= []
    accu1 = []
    accu2 = []
    accu3 = []
    
    current=start
    end_to_import = datetime(end.year,end.month,end.day,0,0)
    while ((timelist==[]) or (timelist[-1]<end_to_import )):
        path = '/home/rml000/sitestore8/UQAM/projet_margaux/measurements_papier_margaux/Sorel/rubber_daqy/' + str(current.year)[2:]+str(current.month).zfill(2) + '/'
        filename = 'Pluvio2' + str(current.year)[2:]+str(current.month).zfill(2)+str(current.day).zfill(2)+ '.txt'
        f = open(path+filename,'r')
        for line in f:
            if ((line =='\n') or (line=='HEADER: item1, item2, item3\n') or ('Error' in line)):
                yo=1
            else:    
                try:
                    linesplit = line.split(';')
		
                    accu1.append(np.array(linesplit[3]).astype(float))
                    accu2.append(np.array(linesplit[4]).astype(float))
                    accu3.append(np.array(linesplit[5]).astype(float))

                    time = datetime.strptime(linesplit[0][0:19],'%Y/%m/%d %H:%M:%S')
                    timelist.append(time)
                    if (timelist[-1]>=end):
                        break
                except:
                    yo=2

        f.close()
        current=current+timedelta(days=1)

    avg_accu = np.mean([np.array(accu2)-np.array(accu2)[0],np.array(accu3)-np.array(accu3)[0]],axis=0 )
    df_pluvio = pd.DataFrame(data=avg_accu,index=timelist,columns=['meanaccu'])
    df_pluvio_resa = df_pluvio.resample('10min').mean()

    intensite = np.gradient(df_pluvio_resa['meanaccu'])*6.
    intensite[intensite<0] = 0
    df_pluvio_resa['intensite'] = intensite

    return df_pluvio_resa

def RH_to_TD(RH, T, P):
    #T in celcius, P in hPa
#     a = 6.1121 mbar, b = 18.678, c = 257.14 C, d = 234.5C
    Rd = 286.9
    Rw = 461.5  # http://www.newworldencyclopedia.org/entry/Gas_constant
    A = 0.61121
    B = 18.678
    C = 234.5
    D = 257.14
    
    Gamma_m = np.log(RH/100. * np.exp( (B-T/C)*(T/(D+T)) ))
    
    TDP = D*Gamma_m/(B-Gamma_m)
    return TDP   

'''
WXT
'''
def import_wxt(start,end):
    timelist = []
    winddir  = []
    windspd  = []
    pressure = []
    temperature = []
    RH  = []
    intensite = []

    def str_to_float(a):
        return np.array(re.sub('[^-0123456789.]','', a))

    current=start
    end_to_import = datetime(end.year,end.month,end.day,0,0)
    while ((timelist==[]) or (timelist[-1]<end_to_import )):


        cnt = 0
        path = '/home/rml000/sitestore8/UQAM/projet_margaux/measurements_papier_margaux/Sorel/rubber_daqy/' + str(current.year)+'_'+str(current.month).zfill(2) + '/raw/'
        filename = 'WXT520' + str(current.year)[2:]+str(current.month).zfill(2)+str(current.day).zfill(2)+ '.txt'
        f = open(path+filename,'r')
        for line in f:
            try:
                if ((line =='\n') or (line=='HEADER: item1, item2, item3\n') or ('Error' in line) or (line=='')):
                    print('skip')
                elif (np.mod(cnt,2)==1):
                    time = datetime.strptime(line[0:19],'%Y/%m/%d %H:%M:%S')
                    timelist.append(time)
                else:    
                    linesplit = line.split(',')
                    winddir.append(str_to_float(linesplit[2]).astype(float))
                    windspd.append(str_to_float(linesplit[5]).astype(float))
                    pressure.append(str_to_float(linesplit[9]).astype(float))
                    temperature.append(str_to_float(linesplit[7]).astype(float))
                    RH.append(str_to_float(linesplit[8]).astype(float))
                    intensite.append(str_to_float(linesplit[12]).astype(float))
                
                    
                    if (timelist[-1]>=end):
                        break
                cnt =cnt+1
            except:		
                cnt = cnt+ 1

        current = current + timedelta(days=1)
	
        f.close()
        
    winddir = np.array(winddir)
    windspd = np.array(windspd)
    pressure = np.array(pressure)
    temperature = np.array(temperature)
    RH = np.array(RH)
    intensite = np.array(intensite)
    TD = RH_to_TD(RH, temperature, pressure)

    df_wxt = pd.DataFrame(data={'winddir':winddir,'windspd':windspd,'pressure':pressure,'temperature':temperature,'RH':RH,'intensite':intensite,'TD':TD},index=timelist) 
    df_wxt = df_wxt[slice(start,end)]
    return df_wxt 


