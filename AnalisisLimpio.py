# -*- coding: utf-8 -*-
"""
Code to analyze ang-V curves

@author: Daniel Elvira
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import Polarization as pol

# =============================================================================
# Change font sizes
# =============================================================================

SMALL_SIZE = 8
MEDIUM_SIZE = 11
BIGGER_SIZE = 16

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

###############################################################################


pathData = r"C:\Users\USUARIO\Desktop\TFGEnd\DataPol\Data"+"\\"
pathAnalisis = r"C:\Users\USUARIO\Desktop\TFGEnd\DataPol\Analisis"+"\\"


# Funciones útiles

def printIf(string, boolean):
    """
    print if boolean is true
    """
    if boolean:
        print(string)
    else:pass
        

def find_nearest(array, value):
    idx = np.argmin(np.abs(np.array(array)-value))
    return idx

def normalize(array):
    return array/np.sqrt(np.min(array)**2+np.max(array)**2)

def findMinMax(angle,voltage):
    """
    It finds the peaks and valleys of the data of polarization.
    Returns: 4 arrays the first 2 with the position and value of the valleys and the other two with the peaks
    """
    max_index = find_peaks(voltage)[0]
    min_index = find_peaks(-voltage)[0]
    return angle[min_index],voltage[min_index], angle[max_index],voltage[max_index]

def fix(angle, voltage, totalExtr, extreme):
    """
    Solves the problem of findMinMax that it gets false peaks or valley
    
    totalExt: Number of real peaks or valleys there should be
    extreme: Solves for: 0 -> Valleys ; 1 -> Peaks
    """
    if extreme == 0: 
        ext = np.argmax
    else: 
        ext = np.argmin
    angfixed = np.copy(angle)
    Vfixed = np.copy(voltage)
    while totalExtr<len(Vfixed):
        index_Vext = ext(Vfixed)
        angfixed = np.delete(angfixed,index_Vext)
        Vfixed = np.delete(Vfixed,index_Vext)
    return angfixed,Vfixed


### Lectura de datos ### 

def readData_Pol(path):
    """
    It reads data from LabVIEW KinesisMotor/Main
    
    Returns
    -------
    angle : array
    voltage : array
    """
    f = open (path,'r')
    rows = f.read().split("\n")
    f.close()
    N = len(rows)-1 #Last row is empty
    angle = np.zeros(N)
    voltage = np.zeros(N)
    for i in range(N):
        AllData = rows[i].split("\t")
        angle[i] = float(AllData[1].replace(",","."))
        voltage[i] = float(AllData[2].replace(",","."))
    
    return angle, voltage



### Funciones auxiliares para el analisis ###

## Cálculo

def Bias(voltage,pathNoise="C:\\Users\\Usuario\\Documents\\Daniel\\Data\\Background_1.dat",prompt=True):
    """
    Parameters
    ----------
    voltage : array 1xN
        Values of voltage from which you want to eliminate the background.
    pathNoise : string, optional
        path to the .dat that has the measurement of the background noise.
        The default is "C:\\Users\\Usuario\\Documents\\Daniel\\Data\\Background_1.dat".

    Returns
    -------
    array 1xN
        voltage array with the average voltage noise (from pathNoise) substacted.
    """
    kk, VBias = readData_Pol(pathNoise)
    VBiasAv = np.sum(VBias)/len(VBias)
    printIf(f"Background Noise: {VBiasAv} V",prompt)
    return voltage-VBiasAv

def fitVoltage(angle,S0,S1,S2,S3,analyzer,args):
    """
    Function that is used in the fit analysis.
    """
    S_input = np.array([S0,S1,S2,S3])
    Intensidad = np.zeros(len(angle))
    for i in range(len(angle)):
        Intensidad[i] = analyzer(angle[i], S_input,args) 
    return Intensidad

## Representación

def AjusteVueltas(angle,voltage):
    """
    Moves Experimental data to the angle interval (0º,360º)
    Returns 2 2D-lists: 
    angAjust: List of sublist where each sublist represents a lap and
    the sublists contain the angle data in the interval (0º,360º)
    
    Vvuelta: List of sublist where each sublist represents a lap and
    the sublists contain the voltage data
    """
    angVueltas = (angle//360).astype(int)
    angAjust = []
    Vvuelta = []
    for i in range(angVueltas[-1]+1):
        angAjust.append([])
        Vvuelta.append([])
    for i in range(len(angVueltas)):
        angAjust[angVueltas[i]].append(angle[i]%360)
        Vvuelta[angVueltas[i]].append(voltage[i])
    return angAjust,Vvuelta


def plotVueltas(angle,voltage,pathSave="",titleplot=""):
    """
    Plots Experimental Data but in the same angle interval (0,360) with its corresponding laps
    """
    angleAdjust, Vadjust = AjusteVueltas(angle,voltage)
    plt.figure(f"Vueltas {titleplot}")
    for i in range(len(angleAdjust)):
        plt.plot(angleAdjust[i],Vadjust[i],"-",label=f"Lap {i}, ang: ({i*360},{(i+1)*360})")
    plt.legend(loc="upper right")
    plt.xlabel("Ángulo Analizador (º)")
    plt.ylabel("Voltaje (V)")
    if pathSave!="":
        plt.savefig(pathSave,bbox_inches='tight')
    plt.show()


### Analisis ###

def PeakAnalisis(name, pathIn=pathData):
    print(f"##### {name} #####")
    path = pathIn+name+".dat"
    ang,Volt = readData_Pol(path)
    angMin,vMin,angMax, vMax = findMinMax(ang[:-1], Volt[:-1])
    checkMin = vMin<1. #peaks a veces detecta minimos que no lo son, de ahí esta comprobación
    angMin,vMin = angMin[checkMin],vMin[checkMin]
    angMinMod = angMin%90
    angMaxMod = angMax%90
    vMin_av = np.average(vMin)
    vMax_av = np.average(vMax)
    vMin_std = np.std(vMin)
    vMax_std = np.std(vMax)
    print("\nMinima")
    print("Angles")
    print(angMin)
    print(angMinMod)
    print("Voltage")
    print(vMin)
    print(f"avg: {vMin_av}, std: {vMin_std}")
    
    print("\nMaxima")
    print("Angles")
    print(angMax)
    print(angMaxMod)
    print("Voltage")
    print(vMax)
    print(f"avg: {vMax_av}, std: {vMax_std}")
    badness = np.sqrt(vMin_std**2+vMax_std**2)
    print(f"\n'Badness': {badness}")
    return vMin_av,vMax_av,vMin_std,vMax_std

def Analisis2(fileAnalisis,analyzer,args,fileBackground="Background_1",showFigs=1,saveData=0,prompt=True): #saveData: 0:Nothing, 1:All, 2:txt and Ellipse, 3:Only txt, 4:Only figures
    printIf("\n################"+fileAnalisis+"################\n",True)
    path = pathData+fileAnalisis+".dat"
    pathOut = pathAnalisis
    pathBackground = pathData+fileBackground+".dat"
    if saveData == 0 or saveData==3:
        pathVueltas = ""
        pathElipse = ""
    elif saveData == 1 or saveData == 4:
        pathVueltas = pathOut+"Vueltas"+fileAnalisis+".png"
        pathElipse = pathOut+"Elipse"+fileAnalisis+".png"
    elif saveData == 2:
        pathVueltas = ""
        pathElipse = pathOut+"Elipse"+fileAnalisis+".png"
    else: pass
    ### Read Data ###
    ang, V = readData_Pol(path) #WARNING: NOT IN RADIANS
    V = Bias(V,pathNoise=pathBackground,prompt=prompt)
    
    plotVueltas(ang,V,pathSave=pathVueltas,titleplot=fileAnalisis)
    ## Suponiendo Luz puramente polarizada (S0 = Ip)
    Sg0 = np.sqrt(np.max(V)**2+np.min(V)**2)
    
    ### Buscamos un buen Guess ###
    zero_ref = args[0] #EN RADIANES
    
    if analyzer[0] == "Polarizer":
        # Tomando solo la primera vuelta
        ang0 = find_nearest(ang, 0+zero_ref)
        I1 = V[ang0]
        Sg1 = 2*I1-Sg0
        ang45 = find_nearest(ang, 45+zero_ref)
        I2 = V[ang45]
        Sg2 = 2*I2-Sg0
        Sg3Check = Sg0**2-Sg1**2-Sg2**2
        if Sg3Check>=0:
            Sg3 = -np.sqrt(Sg3Check)
        else:
            Sg3 = 0
            printIf(f"Error en calculo de S3 en {fileAnalisis} (del guess): Sg3^2 = {Sg3Check}\nTomando S3g = 0",prompt)
    elif analyzer[0] == "Halfwave":
        ang1 = find_nearest(ang, ((np.pi/8-zero_ref/2)*180/np.pi)%90)
        I1 = V[ang1]
        Sg1 = 2*I1-Sg0
        ang2 = find_nearest(ang, (-zero_ref/2*180/np.pi)%90)
        I2 = V[ang2]
        Sg2 = 2*I2-Sg0
        Sg3Check = Sg0**2-Sg1**2-Sg2**2
        if Sg3Check>=0:
            Sg3 = -np.sqrt(Sg3Check)
        else:
            Sg3 = 0
            printIf(f"Error en calculo de S3 en {fileAnalisis} (del guess): Sg3^2 = {Sg3Check}\nTomando S3g = 0",prompt)
        
        
        
    else: pass
    Scutre = np.array([Sg0,Sg1,Sg2,Sg3])
    ### Ajustamos al mejor S ###
    
    parameters, covariance = curve_fit(lambda x, s1,s2: fitVoltage(x,Sg0,s1,s2,0,analyzer[1],args[1:]), ang/180*np.pi, V, p0=Scutre[1:-1])
    # parameters, covariance = curve_fit(lambda x, s0,s1,s2: fitVoltage(x,s0,s1,s2,0,analyzer,args), ang/180*np.pi, V, p0=Scutre[:-1])
    fS0 = Sg0 #parameters[0]
    fS1 = parameters[0]
    fS2 = parameters[1]
    fS3 = np.sqrt(fS0**2-fS1**2-fS2**2) #parameters[2]
    
    abs_sigma = np.sqrt(np.diag(covariance))
    errTot = np.linalg.norm(abs_sigma)
    fS = np.array([fS0,fS1,fS2,fS3])
    printIf(f"Parámetro de Stokes Guess: {np.round(Scutre,4)}, Diferencia: {Sg0**2-Sg1**2-Sg2**2-Sg3**2:.4f}",prompt)
    printIf(f"Parámetro de Stokes Fit: {np.round(fS,4)}, Diferencia: {fS0**2-fS1**2-fS2**2-fS3**2:.4f}",prompt)
    printIf(f"Error: {abs_sigma}, Error Total: {errTot}",prompt)
    fIntensidad = np.zeros(len(ang))
    for i in range(len(ang)):
        fIntensidad[i] = analyzer[1](ang[i]*np.pi/180, fS,(zero_ref,))
        
    plt.figure(f"Comparativa {fileAnalisis}")
    plt.plot(ang,V,"o",label="Experimental",markersize=1)
    plt.plot(ang,fIntensidad,"--",label="Ajuste")
    plt.legend(loc="upper right")
    plt.xlabel("Ángulo Analizador (º)")
    plt.ylabel("Voltaje (V)")
    if saveData == 1 or saveData == 4:
        plt.savefig(pathOut+"Comparativa"+fileAnalisis+".png",bbox_inches='tight')
    plt.show()
    
    A,B,angElipse,sentido = pol.ploteo(fS, title=fileAnalisis,pathSave=pathElipse,printt=prompt)
    if saveData == 1 or saveData == 2 or saveData==3:
        with open(pathOut+"Data"+fileAnalisis+".txt", 'w') as f:
            f.write(f"Parametros de Stokes:\n{fS[0]},{fS[1]},{fS[2]},{fS[3]}\n")
            f.write(f"Elipse de Polarizacion:\nSemieje Mayor: {A}, Semieje Menor: {B}, Excentricidad: {np.sqrt(np.abs(1-B*B/A/A))},\nAngulo(º): {(angElipse*180/np.pi)}, Sentido: {sentido}\n")
        
    if showFigs==0:
        plt.close("all")
    elif showFigs==2:
        plt.close(f"Vueltas {fileAnalisis}")
        plt.close(f"Comparativa {fileAnalisis}")
    elif showFigs==3:
        plt.close(f"Comparativa {fileAnalisis}")
    elif showFigs==4:
        plt.close(f"Vueltas {fileAnalisis}")
        plt.close(f"Eliptica {fileAnalisis}")
    else:
        pass
    return ang, V, fS, Scutre

def Badness_Analize(medida,numberExtremes ,analyzer,showFigs=1,saveFig=0):
    step = 10
    stdB_array = np.array([])
    stdP_array = np.array([])
    stdV_array = np.array([])
    for i in range(0,360,step):
        ang, Volt = readData_Pol(pathData+f"{medida}{i}_2.dat")
        angMin,vMin,angMax, vMax = findMinMax(ang[:-1], Volt[:-1])
        angMin,vMin = fix(angMin,vMin,numberExtremes,0)
        angMax,vMax = fix(angMax, vMax,numberExtremes,1)
        vMin_av = np.average(vMin)
        vMax_av = np.average(vMax)
        vMin_std = np.std(vMin)
        vMax_std = np.std(vMax)
        stdP_array = np.append(stdP_array,np.sqrt(vMax_std**2)) # Ranking by peaks 
        stdV_array = np.append(stdV_array,np.sqrt(vMin_std**2)) # Ranking by valleys
        Norm = np.sqrt(vMin_av**2+vMax_av**2)
        stdB_array = np.append(stdB_array,np.sqrt(vMax_std**2+vMin_std**2)/Norm) # Ranking by badness
    
    
    
    std_array2 = np.copy(stdB_array)
    data = np.arange(0,360,step,dtype=int)
    for i in range(1,len(stdB_array)+1):
        worst = np.argmax(std_array2)
        if i==1:
            Analisis2(f"{medida}{data[worst]}_2",analyzer,(0,0),showFigs=(4),saveData=(0),prompt=(False))
        ang, Volt = readData_Pol(pathData+f"{medida}{data[worst]}_2.dat")
        angMin,vMin,angMax, vMax = findMinMax(ang[:-1], Volt[:-1])
        print(f"\nPeor medida: Top {i}",data[worst])
        print(vMax)
        print(vMin)
        print(std_array2[worst])
        std_array2 = np.delete(std_array2, worst)
        data = np.delete(data,worst)
    
        
    ang_inc = np.arange(0,360,step,dtype=int)
    
    

    plt.figure(f"Peaks {medida}")
    plt.plot(ang_inc,stdP_array, "o--")
    # English
    plt.xlabel("$\\lambda$/2 Angle(º)")
    plt.ylabel("Peak error")
    # Español
    plt.xlabel("Ángulo $\\lambda$/2 (º)")
    plt.ylabel("Error en los picos")
    plt.grid()
    if saveFig == 1:
        plt.savefig(pathAnalisis+"Resumenes\\Peaks_"+medida+".png",bbox_inches='tight')
    else:
        plt.title(medida)
    plt.show()
    
    plt.figure(f"Valleys {medida}")
    plt.plot(ang_inc,stdV_array, "o--")
    # English
    plt.xlabel("$\\lambda$/2 Angle(º)")
    plt.ylabel("Valley error")
    # Español
    plt.xlabel("Ángulo $\\lambda$/2 (º)")
    plt.ylabel("Error en los valles")
    plt.grid()
    if saveFig == 1:
        plt.savefig(pathAnalisis+"Resumenes\\Valleys_"+medida+".png",bbox_inches='tight')
    else:
        plt.title(medida)
    plt.show()
    
    
    plt.figure(f"Badness {medida}")
    plt.plot(ang_inc,stdB_array, "o--")
    # English
    plt.xlabel("$\\lambda$/2 Angle(º)")
    plt.ylabel("Badness")
    # Español
    plt.xlabel("Ángulo $\\lambda$/2 (º)")
    plt.ylabel("Error Extremos")
    plt.grid()
    if saveFig == 1 or saveFig == 2:
        plt.savefig(pathAnalisis+"Resumenes\\Badness_"+medida+".png",bbox_inches='tight')
    else:
        plt.title(medida)
    plt.show()
    
    if showFigs == 2:
        plt.close(f"Valleys {medida}")
        plt.close(f"Peaks {medida}")
    elif showFigs==0:
        plt.close(f"Valleys {medida}")
        plt.close(f"Peaks {medida}")
        plt.close(f"Badness {medida}")
    else: pass
    return stdB_array


#%%

# Para ver una medida

Analisis2(f"Definitivo0_{1}",["Polarizer",pol.analyzer_Simple],(0,0),showFigs=(1),saveData=(0),prompt=(True))

# Para ver un conjunto de ellas
medida = "SinNadaGiradoC"
print(f"## {medida} ##")

for i in range(0,360,90):
    Analisis2(f"{medida}{i}_1",["Halfwave",pol.analyzer_Halfwave],(0,0),showFigs=(1),saveData=(0),prompt=(True))

#%%

# Analisis de los picos (Para una medida cualquiera)

Badness_Analize("DefinitivoSinObjetivo", 2, ["Polarizer",pol.analyzer_Simple],showFigs=2 ,saveFig=2)
Badness_Analize("Definitivo", 2, ["Polarizer",pol.analyzer_Simple], showFigs=2,saveFig=2)
Badness_Analize("SinNadaGiradoC", 4, ["Halfwave",pol.analyzer_Halfwave],showFigs=2, saveFig=2)


#%%

### Para guardar analisis

# for i in range(0,360,20):
#     Analisis2(f"ConObjetivo{i}_1",["Halfwave",pol.analyzer_Halfwave],(0,0),showFigs=(0),fileBackground=f"Background{i//20*20}_1",saveData=(1),prompt=(False))

# for i in range(0,360,10):
#     Analisis2(f"FotodetectorA{i}_1",["Halfwave",pol.analyzer_Halfwave],(0,np.pi/2),showFigs=(4),fileBackground=f"Background{i//20*20}_1",saveData=(1),prompt=(False))


# for j in range(1,6):
#     for i in range(0,360,10):
#         Analisis2(f"Completo{i}_{j}",["Halfwave",pol.analyzer_Halfwave],(0,0),showFigs=(0),fileBackground=f"Background{i//20*20}_1",saveData=(1),prompt=(False))

# for j in range(1,6):
#     for i in range(0,360,10):
#         Analisis2(f"SinNadaGiradoC{i}_{j}",["Halfwave",pol.analyzer_Halfwave],(0,0),showFigs=(0),saveData=(3),prompt=(False))

# for j in range(1,6):
#     for i in range(0,360,10):
#         Analisis2(f"Definitivo{i}_{j}",["Polarizer",pol.analyzer_Simple],(0,0),showFigs=(0),saveData=(3),prompt=(False))

# for j in range(2,7):
#     for i in range(0,360,10):
#         Analisis2(f"DefinitivoSinObjetivo{i}_{j}",["Polarizer",pol.analyzer_Simple],(0,0),showFigs=(0),saveData=(2),prompt=(False))






