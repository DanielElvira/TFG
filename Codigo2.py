# -*- coding: utf-8 -*-
"""
Segundo codigo 

@author: Daniel Elvira López
"""
import os
import glob
import Polarization as pol
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import linregress

# =============================================================================
# Change font sizes
# =============================================================================

SMALL_SIZE = 8
MEDIUM_SIZE = 11
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

###############################################################################


pthAn = r"C:\Users\USUARIO\Desktop\TFGEnd\DataPol\Analisis"+"\\"
#Funciones útiles

def sortAng(ang,others):
    """
    Ordena un array (ang) de menor a mayor y reordena los arrays de la lista
    others de forma que mantengas la correspondecia con el array ang inicial,
    ang[i] = sr_ang[j] --> other[i] = sr_other[j]
    """
    arr1inds = ang.argsort()
    sr_ang = ang[arr1inds]
    sr_others = []
    for i in range(len(others)):
        sr_others.append(others[i][arr1inds])
    return sr_ang,sr_others

def Number_Error(number,error):#│Los 0's del final no los pone. Ej: 4.698, 0.021 --> 4.7 +/- 0.02
    expNum = int(np.floor(np.log10(number))) #Tells you the exponent of number 10 requiered
    expErr = int(np.floor(np.log10(error)))
    numCie = np.round(number/10**expNum,abs(expErr-expNum))
    errCie = int(error//(10**expErr))
    if expErr==0:
        if expNum == 0:
            string = f"{int(numCie)} +/- {errCie}"
        else:
            string = f"{int(numCie)}E{expNum} +/- {errCie}"
    else:
        if abs(expNum) <= 3:
            Numstr = f"{round(numCie*10**expNum,abs(expErr))}"
        else:
            Numstr = f"{numCie}E{expNum}"
        if abs(expErr) <= 3:
            Errstr = f"{round(errCie*10**expErr,abs(expErr))}"
        else:
            Errstr = f"{errCie}E{expErr}"
        string = f"{Numstr} +/- {Errstr}"
    return [numCie*10**expNum,errCie*10**expErr],string


def add_subplot_axes(ax,rect,facecolor='w'): #Embedded subplots from : https://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],facecolor=facecolor)  # matplotlib 2.0+
    # subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


### Agrupación de los datos en el Summary

def MakeSummary(pathAnalisis):
    path = pathAnalisis+"*.txt" #https://www.geeksforgeeks.org/list-all-files-of-certain-type-in-a-directory-using-python/
    files = glob.glob(path) #name of txt files
    
    for i in range(len(files)):
        files[i] = files[i].replace(pathAnalisis+"Data","").replace(".txt","")
    
    ### Makes a list of the unique names (without taking into account what goes after the _) in the files 
    files_names = [] #List with all the unique names
    for file in files:
        unico = True
        name = file.split("_")[0]
        for i in range(len(files_names)):
            if name == files_names[i]:
                unico = False
        if unico:
            files_names.append(name)
        else: pass
    
    measures_names = []
    for file in files_names:
        name_measu = ''.join([i for i in file if not i.isdigit()])
        measures_names.append(name_measu)
        
    measures_names = list(dict.fromkeys(measures_names))
    print(measures_names)
        
    
    with open(pathAnalisis+"Resumenes\\Summary.txt", 'w') as f:
        for group in files_names:
            #Getting all the files of that measurement
            path2 = pathAnalisis+f"Data{group}_*.txt"
            files2 = glob.glob(path2) #All the paths of the files that share the same name. Ex.: Test_1, Test_2, Test_10
            
            #Title of each type of measurement
            f.writelines(f"##### {group} #####\n")
            
            #Obtaining and calculating each Stokes parameter and its average
            S_tot = np.zeros((len(files2),4))
            for i in range(len(files2)):
                file = files2[i]
                with open(file, 'r') as g:
                    g.readline()
                    segLinea = g.readline()
                    S_array = segLinea.split(",")
                    S_array = np.array([float(ele) for ele in S_array])
                    # print(S_array)
                    S_tot[i] = S_array
            S_avg = np.average(S_tot, axis=0)
            S_std = np.std(S_tot,axis=0)
            f.writelines("\n--- Stockes Parameters: Average and std ---\n")
            f.writelines(f"{S_avg[0]},{S_avg[1]},{S_avg[2]},{S_avg[3]}\n")
            f.writelines(f"{S_std[0]},{S_std[1]},{S_std[2]},{S_std[3]}\n")
            
            #Calculating the Ellipse
            A,B,ang,sentido = pol.Ellipse(S_avg)
            a = A/2 #Semieje mayor
            b = B/2 #Semieje menor
            e = np.sqrt(np.abs(1-b*b/a/a))
            
            f.writelines("\n--- Ellipse ---\n")
            f.writelines(f"Mayor Axis {A}, Minor Axis {B},\nAngle(º) {(ang*180/np.pi)%180}, Direction {sentido}\n")
            f.writelines(f"Excentricidad {e}\n")
            
            
            f.writelines(f"Files: {len(files2)}\n")
            f.writelines("\n")
    print("\n############# Summary created #############")



### Analisis

def ExcentricityAnalisis(name_analisis,pathAnalisis=pthAn):
    angle_array = np.array([])
    exc_array= np.array([])
    error_array = np.array([])
    angOut_array = np.array([])
    errang_array = np.array([])
    with open(pathAnalisis+"Resumenes\\Summary.txt", 'r') as fi:
        for ln in fi:
            if ln.startswith("#"):
                ### Read which meassure correspond to each line ###
                lne = ln.replace("#","")
                lne = lne.strip()
                # print(lne)
                name_measu = ''.join([i for i in lne if not i.isdigit()]) #Le quita los numeros a una string
                angle_in = lne.replace(name_measu,"")
                if angle_in == "": angle_in = 0 #SOLUCION TEMPORAL PARA LAS MEDIDAS ANTIGUAS
                angle_in = int(angle_in)
                # print(name_measu)
                # print(angle_in)
                if name_measu == name_analisis: #Check if it is the name you wanted
                    angle_array= np.append(angle_array,angle_in)
                    fi.readline()
                    fi.readline()
                    terLinea = fi.readline()
                    cuaLinea = fi.readline()
                    S_array = terLinea.split(",")
                    S_array = np.array([float(ele) for ele in S_array])
                    std_array = cuaLinea.split(",")
                    std_array = np.array([float(ele) for ele in std_array])
                    # print(S_array, std_array)
                    kk,kk2,angOut,kk3 = pol.Ellipse(S_array)
                    kk,kk2,e_angOut = pol.EllipseUncertainty(S_array,std_array)
                    exc,e_exc = pol.excentricity(S_array,std_array)
                    
                    exc_array = np.append(exc_array, exc)
                    error_array = np.append(error_array, e_exc)
                    angOut_array = np.append(angOut_array, angOut)
                    errang_array = np.append(errang_array, e_angOut)
        angle_array, sorted_list = sortAng(angle_array,[exc_array, error_array,angOut_array,errang_array])
        exc_array, error_array,angOut_array,errang_array = sorted_list[0],sorted_list[1],sorted_list[2],sorted_list[3]       
    return angle_array, exc_array, error_array,angOut_array,errang_array #ANGULOS OUT EN RADIANES

#Plots del analisis

def EllipseMeasure(name_measure,pathAnalisis=pthAn):
    with open(pathAnalisis+"Resumenes\\Summary.txt", 'r') as fi:
        for ln in fi:
            if ln.startswith("#"):
                ### Read which meassure correspond to each line ###
                lne = ln.replace("#","")
                lne = lne.strip()
                if lne == name_measure: #Check if it is the name you wanted
                    fi.readline()
                    fi.readline()
                    terLinea = fi.readline()
                    cuaLinea = fi.readline()
                    S_array = terLinea.split(",")
                    S_array = np.array([float(ele) for ele in S_array])
                    std_array = cuaLinea.split(",")
                    std_array = np.array([float(ele) for ele in std_array])
                    pol.ploteo(S_array,eS=std_array ,title=name_measure, pathSave=pathAnalisis+f"Resumenes\\Ellipse{name_measure}",printt=False)

def plotExcAnalisis(name_an,saveFig=0):
    ang,exc, e_exc, angOut, e_angOut = ExcentricityAnalisis(name_an)
    plt.figure(f"{name_an}_Excentricidad")
    if all(e_exc == np.zeros(len(e_exc))):
        plt.plot(ang,exc,marker="o",linestyle = "--",label="Nada")
    else:
        plt.errorbar(ang,exc,yerr=e_exc,capsize=10,marker="o",linestyle = "--",label="Nada")
    plt.xlabel("Ángulo $\\lambda$/2 (º)")
    plt.ylabel("Excentricidad")
    plt.grid()
    if saveFig==1:# I do not want titles in the downloaded figures
        plt.savefig(pthAn+"Resumenes\\Exc_"+name_an+".png",bbox_inches='tight')
    else:
        plt.title(f"{name_an} Excentricidad")
    plt.show()
    
    plt.figure(f"{name_an}_Excentricidad Modulo")
    if all(e_exc == np.zeros(len(e_exc))):
        plt.plot(ang%180,exc,"o",label="Nada")
    else:
        plt.errorbar(ang%180,exc,yerr=e_exc,capsize=10,marker="o",linestyle="",label="Nada")
    plt.xlabel("Ángulo $\\lambda$/2 (º)")
    plt.ylabel("Excentricidad")
    plt.grid()
    if saveFig==1:# I do not want titles in the downloaded figures
        plt.savefig(pthAn+"Resumenes\\ExcMod_"+name_an+".png",bbox_inches='tight')
    else:
        plt.title(f"{name_an}: Excentricidad Modulo")
    plt.show()

def plotAngAnalisis(name_an,saveFig=0):
    ang,exc, e_exc, angOut, e_angOut = ExcentricityAnalisis(name_an)
    
    pks = find_peaks(angOut)[0]
    angOut_ad = np.copy(angOut)
    for i in range(len(pks)):
        angOut_ad[pks[i]+1:] += np.pi
    
    angOut_ad = angOut_ad*180/np.pi
    result= linregress(ang,angOut_ad)
    m = result.slope
    b = result.intercept
    e_m = result.stderr
    e_b = result.intercept_stderr
    r2 = result.rvalue**2

    linea = b+m*ang
    angDif = angOut_ad-linea
    e_angDif = e_angOut # Desprecio el error que tiene la regresión ya que es muy precisa
    
    if saveFig!= 2:
        # Linear Regression
        plt.figure(f"Linear regression {name_an}")
        plt.errorbar(ang,angOut_ad,yerr=e_angOut,capsize=10,marker="o",linestyle="",label="",zorder=0)
        plt.plot(ang,linea,"-",label=f"({Number_Error(m, e_m)[1]})x+({Number_Error(b, e_b)[1]}), R^2 = {r2:.6f}",zorder=10)
        # English
        #plt.xlabel("$\\lambda$/2 Angle (º)")
        #plt.ylabel("Angle of the light (º)")
        # Español
        plt.xlabel("Ángulo $\\lambda$/2 (º)")
        plt.ylabel("Ángulo de polarización inicial (º)")  
        plt.grid()
        plt.legend()
        if saveFig==1:# I do not want titles in the downloaded figures
            plt.savefig(pthAn+"Resumenes\\LR"+name_an+".png",bbox_inches='tight')
        else:
            plt.title(f"Linear regression {name_an}")
        plt.show()
        
        # Angle Difference
        plt.figure(f"Deviation {name_an}")
        plt.errorbar(ang,angDif,yerr=e_angDif,capsize=10,marker="o",linestyle="--",label="")
        # English
        #plt.xlabel("$\\lambda$/2 Angle (º)")
        # plt.ylabel("Deviation from regression (º)")
        # Español
        plt.xlabel("Ángulo $\\lambda$/2 (º)")
        plt.ylabel("Desviación de la regresión (º)")
        plt.grid()
        if saveFig==1:# I do not want titles in the downloaded figures
            plt.savefig(pthAn+"Resumenes\\Deviation"+name_an+".png",bbox_inches='tight')
        else:
            plt.title(f"Deviation {name_an}")
        plt.show()
    else: #Figures Embeded
        fig = plt.figure(figsize=(6.4,5))
        ax = fig.add_subplot(111)
        # Main Plot
        ax.errorbar(ang,angDif,yerr=e_angDif,capsize=10,marker="o",linestyle="--",label="")
        ax.grid()
        # English
        # ax.set_xlabel("$\\lambda$/2 Angle (º)")
        # ax.set_ylabel("Deviation from regression (º)")
        # Español
        ax.set_xlabel("Ángulo $\\lambda$/2 (º)")
        ax.set_ylabel("Desviación de la regresión (º)")
        # Subplot
        rect = [0.3,0.7,0.3,0.3]
        subax = add_subplot_axes(ax,rect)
        subax.errorbar(ang,angOut_ad,yerr=e_angOut,markersize=4,capsize=3,marker="o",linestyle="",label="",zorder=0)
        subax.plot(ang,linea,"-", linewidth=1.5,label=f"({Number_Error(m, e_m)[1]})x+({Number_Error(b, e_b)[1]}), R^2 = {r2:.6f}",zorder=10)
        # English
        # subax.set_xlabel("$\\lambda$/2 Angle (º)")
        # subax.set_ylabel("Angle of the light (º)")
        # Español
        subax.set_xlabel("Ángulo $\\lambda$/2 (º)", fontsize=BIGGER_SIZE*rect[2]**0.5)
        subax.set_ylabel("Polarización inicial (º)", fontsize=BIGGER_SIZE*rect[3]**0.5)  
        subax.legend(loc="lower left",bbox_to_anchor=(-0.25, 1.),fontsize=BIGGER_SIZE*rect[2]**0.5)
        subax.grid()
        plt.savefig(pthAn+"Resumenes\\Embedded"+name_an+".png",bbox_inches='tight')
        plt.show()

#%%

pathAn = pthAn
MakeSummary(pathAn)


#%%
EllipseMeasure("COC")

#%%
plotAngAnalisis("SinNadaGiradoC", saveFig=(1))
plotAngAnalisis("Definitivo",saveFig=(1))
plotAngAnalisis("DefinitivoSinObjetivo",saveFig=(1))
plotAngAnalisis("DefinitivoSinObjetivo",saveFig=(2))

#%%
# plotExcAnalisis("Giro")
# plotExcAnalisis("SinNadaGiradoC(Sens)")
# plotExcAnalisis("ConObjetivo")
# plotExcAnalisis("ConObjetivoGirado")
# plotExcAnalisis("ConObjetivoGiradoB")
# plotExcAnalisis("ConObjetivoGiradoC") 
# plotExcAnalisis("FotodetectorA")
# plotExcAnalisis("Completo")

plotExcAnalisis("SinNadaGiradoC",saveFig=(1))
plotExcAnalisis("Definitivo",saveFig=(1))
plotExcAnalisis("DefinitivoSinObjetivo",saveFig=(1))
#%%
from scipy.fft import rfft,irfft


# =============================================================================
# Convolucion
# =============================================================================

def Objetivo(name):
    ang,exc, e_exc, angOut, e_angOut = ExcentricityAnalisis(name)

    arr1inds = ang.argsort()
    sr_ang = ang[arr1inds]
    sr_exc = exc[arr1inds]


    angL,excL, e_excL, angOutL, e_angOutL = ExcentricityAnalisis("SinNadaGiradoC")
    arr1inds = angL.argsort()
    sr_angL = angL[arr1inds]
    sr_excL = excL[arr1inds]

    P = rfft(sr_exc)
    F = rfft(sr_excL)
    G = P/F
    g_t = irfft(G)
    return sr_ang,sr_exc, g_t

angA,excA,gA = Objetivo("ConObjetivoGirado")
angB,excB,gB = Objetivo("ConObjetivoGiradoB")
angC,excC,gC = Objetivo("ConObjetivoGiradoC")
plt.figure("Todo")
plt.plot(angA-35,excA,"o--",label="20")
plt.plot(angB-25,excB,"o--",label="40")
plt.plot(angC+0,excC,"o--",label="90")
plt.legend()
plt.show()

plt.figure("Objectivo")
plt.plot(angA-35,gA,"o--",label="20")
plt.plot(angB-25,gB,"o--",label="40")
plt.plot(angC+0,gC,"o--",label="90")
plt.legend()
plt.show()


#%%

#Resumen Analizador 1

names = ["SOB","COB","SOC","COC"]
exc = [0.9976,0.997,0.999,0.997]
err = [9e-4,1e-3,0.003,8e-4]

plt.figure()
plt.errorbar(names,exc,yerr=err,capsize=10,marker="o",linestyle = "None")
plt.ylim((None,1))
plt.xlabel("Tipos de medidas")
plt.ylabel("Excentricidad")
plt.title("Resumen analizador 1")
plt.show()
