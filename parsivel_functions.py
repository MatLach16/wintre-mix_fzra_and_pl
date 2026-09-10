import pandas as pd
import datetime
from datetime import datetime,  timedelta
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches 
import matplotlib.dates as mdates
import os


def parsi_import(data_file, start_date, end_date, interval = 'day', ext_str='.dat'): 
    
    Dbins=[0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.25,3.75,4.25,4.75,5.5,6.5,7.5,8.5,9.5,11,13,15,17,19,21.5,24]     
    vTbins=[0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.10,1.30,1.5,1.7,1.9,2.2,2.6,3,3.4,3.8,4.4,5.2,6,6.8,7.6,8.8,10.4,12,13.6,15.2,17.6,20.80]

    dgrid=np.zeros(len(Dbins))
    vgrid=np.zeros(len(vTbins))

    dgrid[0] = 0+ Dbins[0]/2.
    vgrid[0] = 0+ vTbins[0]/2.
    for i in range(0,len(Dbins)):
        if i<(len(Dbins)-1):
            dgrid[i+1]=Dbins[i]+(Dbins[i+1]-Dbins[i])/2.
            vgrid[i+1]=vTbins[i]+(vTbins[i+1]-vTbins[i])/2.

    dgrid,vgrid = np.meshgrid(dgrid,vgrid)

    current_day = start_date
    end_day     = end_date    
    appended_data = []
    mat = []
    intensite=[]
    temps=[]
    current_in_file = start_date
    
    fichier_to_open = data_file+str(current_day.year)+str(current_day.month).zfill(2)+str(current_day.day).zfill(2)+ext_str
     

    while ((current_in_file<=end_date) & (os.path.exists(fichier_to_open))) :

        fichierIn = open(fichier_to_open,'rb')

        i=0
        for line in  fichierIn:
            i+=1 
            ligne = str(line)

            if i==1 :
                titleTab = ligne.split(";")               
            elif (i-1)%3==0 :
                    continue
            elif i%3==0 :
                    continue

            else :
                ch = ligne.split("<SPECTRUM>")
                chDesc = ch[0]
                tabDesc = chDesc.split(";") 
                tabDesc[0] = tabDesc[0].replace(".","-")                               # extraction de la date et du temps
                date_time_str = tabDesc[0] + ' ' + tabDesc[1] 

                date_time_obj = datetime.strptime(date_time_str[2:], '%d-%m-%Y %H:%M:%S')
                current_in_file = date_time_obj
                
                
                
                if (current_in_file>=start_date) & (current_in_file<=end_date):
                    temps.append(date_time_obj)
                    nbpart = tabDesc[10]
                    intensity = tabDesc[2]
                    intensite.append(float(intensity))

                    chSpec = ch[1]
                    chSpec3 = chSpec.replace("</SPECTRUM>","")    
                    if "ZERO" in chSpec3:

                        tabSpec_N = np.zeros(1024)
                    else:
                        tabSpec = chSpec3.split(";")
                        tabSpec_N = np.zeros(1024)
                        for bin in range(1024):
                            if tabSpec[bin] == '':
                                tabSpec_N[bin] = 0
                            else:
                                try:
                                    tabSpec_N[bin] = int(tabSpec[bin])
                                except:
                                    continue

                    tmp=pd.DataFrame(tabSpec_N).T
                    tmp['Date'] = date_time_obj        
                    tmp.set_index('Date', inplace=True)        
                    appended_data.append(tmp)

                    ar = np.array([[int(nbpart), float(intensity)]])
                    df = pd.DataFrame(ar, columns=['nb_particules', 'Intensity'], index = [date_time_obj])
                    mat.append(df)
                    
        current_day = date_time_obj + timedelta(1)
        fichier_to_open = data_file+str(current_day.year)+str(current_day.month).zfill(2)+str(current_day.day).zfill(2)+ext_str
        
    prec_tab = pd.DataFrame(intensite,index = temps)
    if len(mat)>1:
        mat =  pd.concat(mat)         
    if len(appended_data)>1:
        appended_data = pd.concat(appended_data)   
    return prec_tab, mat, appended_data

def parsi_import_plist_file(start_date, end_date):
    
    times = []
    Int   = []
    Npart = []
    bins  = []
    mat   = []
    print(start_date)
    current = datetime(start_date.year,start_date.month,start_date.day,0)
    while (current < end_date):
        fichier_to_open = f'surface_measurements_wintremix/Trois/parsivel/{current.year}_{current.month:02}/raw/{current.year}_{current.month:02}_{current.day:02}_Parsivel_LOP_Data.txt'
        print(fichier_to_open)
        fichierIn = open(fichier_to_open,'rb')
        line = fichierIn.readline()
        line = line[3:]
        lines = line.split(b"\\r\\n'][b'")
        for linei in lines:
            # print(linei)
            linei_split = linei.split(b';')
            times.append(datetime.strptime(str(linei_split[0]+linei_split[1]), "b'%d.%m.%Y%H:%M:%S'"))
            Int.append(linei_split[5])
            Npart.append(linei_split[6])
            bins.append(linei_split[7:(1024+7)])
        
        current = current+timedelta(days=1)
    print(times[-1])
    prec_tab = pd.DataFrame(np.array(Int).astype(float),index = times)
    appended_data = pd.DataFrame(np.array(bins).astype(float),index=times)
    #prec_tab = pd.to_numeric(prec_tab)
    #appended_data = pd.to_numeric(appended_data)
    return prec_tab,mat, appended_data


def parsi_import_plist_file_hadleigh(data_file, start_date, end_date, interval = 'day', ext_str='.dat'):
    header = ['Date','Time','Station','Sensor temperature','Interval','Intensity of precipitation (mm/h)','Number of detected particle', \
          'V0D0', 'V0D1', 'V0D2', 'V0D3', 'V0D4', 'V0D5', 'V0D6', 'V0D7', 'V0D8', 'V0D9', 'V0D10', 'V0D11', 'V0D12', 'V0D13', 'V0D14', 'V0D15', 'V0D16', 'V0D17', 'V0D18', 'V0D19', 'V0D20', 'V0D21', 
          'V0D22', 'V0D23', 'V0D24', 'V0D25', 'V0D26', 'V0D27', 'V0D28', 'V0D29', 'V0D30', 'V0D31', \
          'V1D0', 'V1D1', 'V1D2', 'V1D3', 'V1D4', 'V1D5', 'V1D6', 'V1D7', 'V1D8', 'V1D9', 'V1D10', 'V1D11', 'V1D12', 'V1D13', 'V1D14', 'V1D15', 'V1D16', 'V1D17', 'V1D18', 'V1D19', 'V1D20', 'V1D21', 
          'V1D22', 'V1D23', 'V1D24', 'V1D25', 'V1D26', 'V1D27', 'V1D28', 'V1D29', 'V1D30', 'V1D31', \
          'V2D0', 'V2D1', 'V2D2', 'V2D3', 'V2D4', 'V2D5', 'V2D6', 'V2D7', 'V2D8', 'V2D9', 'V2D10', 'V2D11', 'V2D12', 'V2D13', 'V2D14', 'V2D15', 'V2D16', 'V2D17', 'V2D18', 'V2D19', 'V2D20', 'V2D21', 
          'V2D22', 'V2D23', 'V2D24', 'V2D25', 'V2D26', 'V2D27', 'V2D28', 'V2D29', 'V2D30', 'V2D31', \
          'V3D0', 'V3D1', 'V3D2', 'V3D3', 'V3D4', 'V3D5', 'V3D6', 'V3D7', 'V3D8', 'V3D9', 'V3D10', 'V3D11', 'V3D12', 'V3D13', 'V3D14', 'V3D15', 'V3D16', 'V3D17', 'V3D18', 'V3D19', 'V3D20', 'V3D21', 
          'V3D22', 'V3D23', 'V3D24', 'V3D25', 'V3D26', 'V3D27', 'V3D28', 'V3D29', 'V3D30', 'V3D31', \
          'V4D0', 'V4D1', 'V4D2', 'V4D3', 'V4D4', 'V4D5', 'V4D6', 'V4D7', 'V4D8', 'V4D9', 'V4D10', 'V4D11', 'V4D12', 'V4D13', 'V4D14', 'V4D15', 'V4D16', 'V4D17', 'V4D18', 'V4D19', 'V4D20', 'V4D21', 
          'V4D22', 'V4D23', 'V4D24', 'V4D25', 'V4D26', 'V4D27', 'V4D28', 'V4D29', 'V4D30', 'V4D31', \
          'V5D0', 'V5D1', 'V5D2', 'V5D3', 'V5D4', 'V5D5', 'V5D6', 'V5D7', 'V5D8', 'V5D9', 'V5D10', 'V5D11', 'V5D12', 'V5D13', 'V5D14', 'V5D15', 'V5D16', 'V5D17', 'V5D18', 'V5D19', 'V5D20', 'V5D21', 
          'V5D22', 'V5D23', 'V5D24', 'V5D25', 'V5D26', 'V5D27', 'V5D28', 'V5D29', 'V5D30', 'V5D31', \
          'V6D0', 'V6D1', 'V6D2', 'V6D3', 'V6D4', 'V6D5', 'V6D6', 'V6D7', 'V6D8', 'V6D9', 'V6D10', 'V6D11', 'V6D12', 'V6D13', 'V6D14', 'V6D15', 'V6D16', 'V6D17', 'V6D18', 'V6D19', 'V6D20', 'V6D21', 
          'V6D22', 'V6D23', 'V6D24', 'V6D25', 'V6D26', 'V6D27', 'V6D28', 'V6D29', 'V6D30', 'V6D31', \
          'V7D0', 'V7D1', 'V7D2', 'V7D3', 'V7D4', 'V7D5', 'V7D6', 'V7D7', 'V7D8', 'V7D9', 'V7D10', 'V7D11', 'V7D12', 'V7D13', 'V7D14', 'V7D15', 'V7D16', 'V7D17', 'V7D18', 'V7D19', 'V7D20', 'V7D21', 
          'V7D22', 'V7D23', 'V7D24', 'V7D25', 'V7D26', 'V7D27', 'V7D28', 'V7D29', 'V7D30', 'V7D31', \
          'V8D0', 'V8D1', 'V8D2', 'V8D3', 'V8D4', 'V8D5', 'V8D6', 'V8D7', 'V8D8', 'V8D9', 'V8D10', 'V8D11', 'V8D12', 'V8D13', 'V8D14', 'V8D15', 'V8D16', 'V8D17', 'V8D18', 'V8D19', 'V8D20', 'V8D21', 
          'V8D22', 'V8D23', 'V8D24', 'V8D25', 'V8D26', 'V8D27', 'V8D28', 'V8D29', 'V8D30', 'V8D31', \
          'V9D0', 'V9D1', 'V9D2', 'V9D3', 'V9D4', 'V9D5', 'V9D6', 'V9D7', 'V9D8', 'V9D9', 'V9D10', 'V9D11', 'V9D12', 'V9D13', 'V9D14', 'V9D15', 'V9D16', 'V9D17', 'V9D18', 'V9D19', 'V9D20', 'V9D21', 
          'V9D22', 'V9D23', 'V9D24', 'V9D25', 'V9D26', 'V9D27', 'V9D28', 'V9D29', 'V9D30', 'V9D31', \
          'V10D0', 'V10D1', 'V10D2', 'V10D3', 'V10D4', 'V10D5', 'V10D6', 'V10D7', 'V10D8', 'V10D9', 'V10D10', 'V10D11', 'V10D12', 'V10D13', 'V10D14', 'V10D15', 'V10D16', 'V10D17', 'V10D18', 'V10D19', 'V10D20', 'V10D21', 
          'V10D22', 'V10D23', 'V10D24', 'V10D25', 'V10D26', 'V10D27', 'V10D28', 'V10D29', 'V10D30', 'V10D31', \
          'V11D0', 'V11D1', 'V11D2', 'V11D3', 'V11D4', 'V11D5', 'V11D6', 'V11D7', 'V11D8', 'V11D9', 'V11D10', 'V11D11', 'V11D12', 'V11D13', 'V11D14', 'V11D15', 'V11D16', 'V11D17', 'V11D18', 'V11D19', 'V11D20', 'V11D21', 
          'V11D22', 'V11D23', 'V11D24', 'V11D25', 'V11D26', 'V11D27', 'V11D28', 'V11D29', 'V11D30', 'V11D31', \
          'V12D0', 'V12D1', 'V12D2', 'V12D3', 'V12D4', 'V12D5', 'V12D6', 'V12D7', 'V12D8', 'V12D9', 'V12D10', 'V12D11', 'V12D12', 'V12D13', 'V12D14', 'V12D15', 'V12D16', 'V12D17', 'V12D18', 'V12D19', 'V12D20', 'V12D21', 
          'V12D22', 'V12D23', 'V12D24', 'V12D25', 'V12D26', 'V12D27', 'V12D28', 'V12D29', 'V12D30', 'V12D31', \
          'V13D0', 'V13D1', 'V13D2', 'V13D3', 'V13D4', 'V13D5', 'V13D6', 'V13D7', 'V13D8', 'V13D9', 'V13D10', 'V13D11', 'V13D12', 'V13D13', 'V13D14', 'V13D15', 'V13D16', 'V13D17', 'V13D18', 'V13D19', 'V13D20', 'V13D21', 
          'V13D22', 'V13D23', 'V13D24', 'V13D25', 'V13D26', 'V13D27', 'V13D28', 'V13D29', 'V13D30', 'V13D31', \
          'V14D0', 'V14D1', 'V14D2', 'V14D3', 'V14D4', 'V14D5', 'V14D6', 'V14D7', 'V14D8', 'V14D9', 'V14D10', 'V14D11', 'V14D12', 'V14D13', 'V14D14', 'V14D15', 'V14D16', 'V14D17', 'V14D18', 'V14D19', 'V14D20', 'V14D21', 
          'V14D22', 'V14D23', 'V14D24', 'V14D25', 'V14D26', 'V14D27', 'V14D28', 'V14D29', 'V14D30', 'V14D31', \
          'V15D0', 'V15D1', 'V15D2', 'V15D3', 'V15D4', 'V15D5', 'V15D6', 'V15D7', 'V15D8', 'V15D9', 'V15D10', 'V15D11', 'V15D12', 'V15D13', 'V15D14', 'V15D15', 'V15D16', 'V15D17', 'V15D18', 'V15D19', 'V15D20', 'V15D21', 
          'V15D22', 'V15D23', 'V15D24', 'V15D25', 'V15D26', 'V15D27', 'V15D28', 'V15D29', 'V15D30', 'V15D31', \
          'V16D0', 'V16D1', 'V16D2', 'V16D3', 'V16D4', 'V16D5', 'V16D6', 'V16D7', 'V16D8', 'V16D9', 'V16D10', 'V16D11', 'V16D12', 'V16D13', 'V16D14', 'V16D15', 'V16D16', 'V16D17', 'V16D18', 'V16D19', 'V16D20', 'V16D21', 
          'V16D22', 'V16D23', 'V16D24', 'V16D25', 'V16D26', 'V16D27', 'V16D28', 'V16D29', 'V16D30', 'V16D31', \
          'V17D0', 'V17D1', 'V17D2', 'V17D3', 'V17D4', 'V17D5', 'V17D6', 'V17D7', 'V17D8', 'V17D9', 'V17D10', 'V17D11', 'V17D12', 'V17D13', 'V17D14', 'V17D15', 'V17D16', 'V17D17', 'V17D18', 'V17D19', 'V17D20', 'V17D21', 
          'V17D22', 'V17D23', 'V17D24', 'V17D25', 'V17D26', 'V17D27', 'V17D28', 'V17D29', 'V17D30', 'V17D31', \
          'V18D0', 'V18D1', 'V18D2', 'V18D3', 'V18D4', 'V18D5', 'V18D6', 'V18D7', 'V18D8', 'V18D9', 'V18D10', 'V18D11', 'V18D12', 'V18D13', 'V18D14', 'V18D15', 'V18D16', 'V18D17', 'V18D18', 'V18D19', 'V18D20', 'V18D21', 
          'V18D22', 'V18D23', 'V18D24', 'V18D25', 'V18D26', 'V18D27', 'V18D28', 'V18D29', 'V18D30', 'V18D31', \
          'V19D0', 'V19D1', 'V19D2', 'V19D3', 'V19D4', 'V19D5', 'V19D6', 'V19D7', 'V19D8', 'V19D9', 'V19D10', 'V19D11', 'V19D12', 'V19D13', 'V19D14', 'V19D15', 'V19D16', 'V19D17', 'V19D18', 'V19D19', 'V19D20', 'V19D21', 
          'V19D22', 'V19D23', 'V19D24', 'V19D25', 'V19D26', 'V19D27', 'V19D28', 'V19D29', 'V19D30', 'V19D31', \
          'V20D0', 'V20D1', 'V20D2', 'V20D3', 'V20D4', 'V20D5', 'V20D6', 'V20D7', 'V20D8', 'V20D9', 'V20D10', 'V20D11', 'V20D12', 'V20D13', 'V20D14', 'V20D15', 'V20D16', 'V20D17', 'V20D18', 'V20D19', 'V20D20', 'V20D21', 
          'V20D22', 'V20D23', 'V20D24', 'V20D25', 'V20D26', 'V20D27', 'V20D28', 'V20D29', 'V20D30', 'V20D31', \
          'V21D0', 'V21D1', 'V21D2', 'V21D3', 'V21D4', 'V21D5', 'V21D6', 'V21D7', 'V21D8', 'V21D9', 'V21D10', 'V21D11', 'V21D12', 'V21D13', 'V21D14', 'V21D15', 'V21D16', 'V21D17', 'V21D18', 'V21D19', 'V21D20', 'V21D21', 
          'V21D22', 'V21D23', 'V21D24', 'V21D25', 'V21D26', 'V21D27', 'V21D28', 'V21D29', 'V21D30', 'V21D31', \
          'V22D0', 'V22D1', 'V22D2', 'V22D3', 'V22D4', 'V22D5', 'V22D6', 'V22D7', 'V22D8', 'V22D9', 'V22D10', 'V22D11', 'V22D12', 'V22D13', 'V22D14', 'V22D15', 'V22D16', 'V22D17', 'V22D18', 'V22D19', 'V22D20', 'V22D21', 
          'V22D22', 'V22D23', 'V22D24', 'V22D25', 'V22D26', 'V22D27', 'V22D28', 'V22D29', 'V22D30', 'V22D31', \
          'V23D0', 'V23D1', 'V23D2', 'V23D3', 'V23D4', 'V23D5', 'V23D6', 'V23D7', 'V23D8', 'V23D9', 'V23D10', 'V23D11', 'V23D12', 'V23D13', 'V23D14', 'V23D15', 'V23D16', 'V23D17', 'V23D18', 'V23D19', 'V23D20', 'V23D21', 
          'V23D22', 'V23D23', 'V23D24', 'V23D25', 'V23D26', 'V23D27', 'V23D28', 'V23D29', 'V23D30', 'V23D31', \
          'V24D0', 'V24D1', 'V24D2', 'V24D3', 'V24D4', 'V24D5', 'V24D6', 'V24D7', 'V24D8', 'V24D9', 'V24D10', 'V24D11', 'V24D12', 'V24D13', 'V24D14', 'V24D15', 'V24D16', 'V24D17', 'V24D18', 'V24D19', 'V24D20', 'V24D21', 
          'V24D22', 'V24D23', 'V24D24', 'V24D25', 'V24D26', 'V24D27', 'V24D28', 'V24D29', 'V24D30', 'V24D31', \
          'V25D0', 'V25D1', 'V25D2', 'V25D3', 'V25D4', 'V25D5', 'V25D6', 'V25D7', 'V25D8', 'V25D9', 'V25D10', 'V25D11', 'V25D12', 'V25D13', 'V25D14', 'V25D15', 'V25D16', 'V25D17', 'V25D18', 'V25D19', 'V25D20', 'V25D21', 
          'V25D22', 'V25D23', 'V25D24', 'V25D25', 'V25D26', 'V25D27', 'V25D28', 'V25D29', 'V25D30', 'V25D31', \
          'V26D0', 'V26D1', 'V26D2', 'V26D3', 'V26D4', 'V26D5', 'V26D6', 'V26D7', 'V26D8', 'V26D9', 'V26D10', 'V26D11', 'V26D12', 'V26D13', 'V26D14', 'V26D15', 'V26D16', 'V26D17', 'V26D18', 'V26D19', 'V26D20', 'V26D21', 
          'V26D22', 'V26D23', 'V26D24', 'V26D25', 'V26D26', 'V26D27', 'V26D28', 'V26D29', 'V26D30', 'V26D31', \
          'V27D0', 'V27D1', 'V27D2', 'V27D3', 'V27D4', 'V27D5', 'V27D6', 'V27D7', 'V27D8', 'V27D9', 'V27D10', 'V27D11', 'V27D12', 'V27D13', 'V27D14', 'V27D15', 'V27D16', 'V27D17', 'V27D18', 'V27D19', 'V27D20', 'V27D21', 
          'V27D22', 'V27D23', 'V27D24', 'V27D25', 'V27D26', 'V27D27', 'V27D28', 'V27D29', 'V27D30', 'V27D31', \
          'V28D0', 'V28D1', 'V28D2', 'V28D3', 'V28D4', 'V28D5', 'V28D6', 'V28D7', 'V28D8', 'V28D9', 'V28D10', 'V28D11', 'V28D12', 'V28D13', 'V28D14', 'V28D15', 'V28D16', 'V28D17', 'V28D18', 'V28D19', 'V28D20', 'V28D21', 
          'V28D22', 'V28D23', 'V28D24', 'V28D25', 'V28D26', 'V28D27', 'V28D28', 'V28D29', 'V28D30', 'V28D31', \
          'V29D0', 'V29D1', 'V29D2', 'V29D3', 'V29D4', 'V29D5', 'V29D6', 'V29D7', 'V29D8', 'V29D9', 'V29D10', 'V29D11', 'V29D12', 'V29D13', 'V29D14', 'V29D15', 'V29D16', 'V29D17', 'V29D18', 'V29D19', 'V29D20', 'V29D21', 
          'V29D22', 'V29D23', 'V29D24', 'V29D25', 'V29D26', 'V29D27', 'V29D28', 'V29D29', 'V29D30', 'V29D31', \
          'V30D0', 'V30D1', 'V30D2', 'V30D3', 'V30D4', 'V30D5', 'V30D6', 'V30D7', 'V30D8', 'V30D9', 'V30D10', 'V30D11', 'V30D12', 'V30D13', 'V30D14', 'V30D15', 'V30D16', 'V30D17', 'V30D18', 'V30D19', 'V30D20', 'V30D21', 
          'V30D22', 'V30D23', 'V30D24', 'V30D25', 'V30D26', 'V30D27', 'V30D28', 'V30D29', 'V30D30', 'V30D31', \
          'V31D0', 'V31D1', 'V31D2', 'V31D3', 'V31D4', 'V31D5', 'V31D6', 'V31D7', 'V31D8', 'V31D9', 'V31D10', 'V31D11', 'V31D12', 'V31D13', 'V31D14', 'V31D15', 'V31D16', 'V31D17', 'V31D18', 'V31D19', 'V31D20', 'V31D21', 
          'V31D22', 'V31D23', 'V31D24', 'V31D25', 'V31D26', 'V31D27', 'V31D28', 'V31D29', 'V31D30', 'V31D31','toDelete','LOP']

    f_list = []
    current = datetime(start_date.year,start_date.month,start_date.day)
    while (current < end_date):

        f = pd.read_csv(data_file+ str(start_date.year)+'_'+str(start_date.month).zfill(2)+'_'+str(start_date.day).zfill(2) +ext_str ,names=header, sep=';', lineterminator=']')
    
        #Note, V0D0 is included unless line is uncommented
        f.drop(['toDelete'], axis = 1, inplace = True)
        #f['V0D0'] = np.nan
        f.reset_index(drop = False, inplace = True)
        f['DateTime'] = f['Date'] + ' ' + f['Time']
        f.drop(['Date','Time'], axis = 1, inplace = True)
        f.index = f['DateTime'].apply(lambda x: datetime.strptime(x, '[b\'%d.%m.%Y %H:%M:%S'))
        f.drop(['DateTime'], axis = 1, inplace = True)
        f.index.rename('Timestamp', inplace=True)
        f.sort_index(inplace=True)
        f.index = f.index.round('min')

        #Drop cols that are not required
        f.drop(['LOP'], axis = 1, inplace = True)
        f.drop(['index'], axis = 1, inplace = True)

        f_list.append(f)
    
    df = pd.concat(f_list)

    return f


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
    yo2.ax.set_yticklabels([r'$10^0$',r'$10^1$',r'$10^2$',r'$10^3$',r'$10^4$'])
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
    yo2.ax.set_yticklabels([r'$10^0$',r'$10^1$',r'$10^2$',r'$10^3$',r'$10^4$'])
    yo2.set_label(r'Count [min$^{-1}$ (m/s)$^{-1}$]')
def parsi_spectrum_plot(start,end,appended_data):
    
#     start = str(datetime.strptime(start_I,"%Y-%m-%d %H:%M:%S")  + timedelta(seconds=avg_time))
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

