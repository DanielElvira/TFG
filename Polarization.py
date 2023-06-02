# -*- coding: utf-8 -*-
"""
@author: Daniel Elvira Lopez
"""
import numpy as np
import matplotlib.pyplot as plt

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





def printIf(string, boolean):
    """
    print if the boolean is true
    """
    if boolean:
        print(string)
    else:pass

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


def LinearPolarizer(a): 
    """a: ángulo de rotación del polarizador"""
    s = np.sin(2*a)
    c = np.cos(2*a)
    return np.array([[1,c,s,0],[c,c*c,s*c,0],[s,s*c,s*s,0],[0,0,0,0]])/2


def LinearRetarder(a,d):
    """
    a: ángulo de rotación del eje rápido (Desde la vertical, creo)
    d: Diferencia de fase entre eje rápido y lento
    d = ¿+/-? pi/2 -> Lambda Cuartos ; d = pi -> Lambda Medios
    """
    s = np.sin(2*a)
    c = np.cos(2*a)
    sd = np.sin(d)
    cd = np.cos(d)
    return np.array([[1,0,0,0],[0,c*c+s*s*cd,c*s*(1-cd),s*sd],[0,c*s*(1-cd),c*c*cd+s*s,-c*sd],[0,-s*sd,c*sd,cd]])


def Ellipse(S): #Wikipedia: Stokes parameters (Define S3= V= +2ABh)
    L = S[1]+1j*S[2]
    A = np.sqrt((S[0]+np.abs(L))/2)
    B = np.sqrt((S[0]-np.abs(L))/2)
    theta = np.angle(L)/2 #L=0 devuelve 0 (Caso circular)
    h = np.sign(S[3])
    return A,B,theta,h


def EllipseUncertainty(S,eS):
    A,B,ang,sentido = Ellipse(S)
    e_A = np.sqrt(eS[0]**2+(S[1]**2*eS[1]**2+S[2]**2*eS[2]**2)/4/(S[1]**2+S[2]**2))/2/A
    e_B = np.sqrt(eS[0]**2+(S[1]**2*eS[1]**2+S[2]**2*eS[2]**2)/4/(S[1]**2+S[2]**2))/2/B
    e_ang = np.sqrt(S[2]**2*eS[1]**2+S[1]**2*eS[2]**2)/2/(S[1]**2+S[2]**2)
    return e_A,e_B,e_ang
    
def excentricity(S,eS):
    A,B,ang,sentido = Ellipse(S)
    e_A,e_B,e_ang = EllipseUncertainty(S,eS)
    exc = np.sqrt(np.abs(1-B*B/A/A))
    e_exc = B/A/A/exc*np.sqrt(B*B/A/A*e_A**2+e_B**2)
    return exc,e_exc


def ploteo(S, eS=np.array([None]),title="", pathSave = "",printt=True):
    A,B,ang,sentido = Ellipse(S)
    a = A/2 #Semieje mayor
    b = B/2 #Semieje menor
    e = np.sqrt(np.abs(1-b*b/a/a))
    printIf(f"Semieje Mayor {A:.4f}, Semieje Menor {B:.4f}, Excentricidad {e:.4f},\nAngulo(º) {(ang*180/np.pi)%180:.4f}, Sentido {sentido}",printt)
    if b != 0:
        theta = np.linspace(0,2*np.pi,1000)
        r = b*np.sqrt(1/(1-(e*np.cos(theta-ang))**2)) #Desde el centro de la elipse 
        #(https://math.stackexchange.com/questions/315386/ellipse-in-polar-coordinates)
        plt.figure(f"Eliptica {title}")
        plt.plot(r*np.cos(theta),r*np.sin(theta))
        plt.axis("equal")
        plt.xlabel("Ex")
        plt.ylabel("Ey")
        if all(eS==None):
            plt.title("A: {:.3f}, B: {:.3f}, e: {:.4f} ,Ang (º): {:.4f}, Sentido: {}".format(A,B,e,(ang*180/np.pi)%180,sentido))
        else:
            eA,eB,e_ang = EllipseUncertainty(S, eS)
            e, eS = excentricity(S, eS)
            plt.title("A: {} # B: {} # e: {}\nAng (º): {} # Sentido: {}".format(Number_Error(A, eA)[1],Number_Error(B, eB)[1],Number_Error(e, eS)[1],Number_Error(ang*180/np.pi%180, e_ang*180/np.pi)[1],sentido))
        if pathSave!="":
            plt.savefig(pathSave,bbox_inches='tight')
        plt.show()
    else:
        x = np.linspace(-a,a,100)
        y = np.tan(ang)*x
        plt.figure(f"Lineal {title}")
        plt.plot(x,y)
        lim = max(a,np.tan(ang)*a)
        plt.xlabel("Ex")
        plt.ylabel("Ey")
        plt.xlim(-lim,lim) #Los limites son para ver bien el angulo de la linea
        plt.ylim(-lim,lim) #Los limites son para ver bien el angulo de la linea
        plt.axis("equal")
        plt.title("A: {:.3f}, Ang (º): {:.2f}".format(A,(ang*180/np.pi)%180)) #Para tenerlo entre (0,180)
        if pathSave!="":
            plt.savefig(pathSave,bbox_inches='tight')
        plt.show()
    return A,B,ang,sentido



def analyzer_Simple(ang,S_input,args):
    """
    The most simple analyzer, just a linear polarizer and a photodetector.
    Note: With this analyzer S_input[3] is irrelevant
    Note 2: args is just to make all analyzers functions work
    (in this one it does not have a purpose)
    """
    Sf = np.matmul(LinearPolarizer(ang),S_input)
    Inte = Sf[0]
    return Inte

def analyzer_Halfwave(ang,S_input,args):
    """
    A more complex analyzer, a halfwave plate a linear polarizer (LP) and a photodetector.
    args[0] is the angle of the fixed LP
    Note: With this analyzer S_input[3] is irrelevant
    """
    polAngle = args[0]
    M = np.matmul(LinearPolarizer(polAngle),LinearRetarder(ang, np.pi))
    Sf = np.matmul(M,S_input)
    Inte = Sf[0]
    return Inte

def SimpleSystem(phi,S_input,args):
    theta,alpha,xi = args
    analyzer = np.matmul(LinearPolarizer(theta),LinearRetarder(phi, np.pi))
    other = np.matmul(LinearRetarder(alpha, np.pi),LinearPolarizer(xi+np.pi/2))
    everything = np.matmul(analyzer,other) 
    Sf = np.matmul(everything,S_input)
    return Sf[0]

def ErraticSystem(phi,S_input,args): #Halfwave plates have slight changes in retardance
    theta,alpha,xi = args
    delta = np.pi+np.abs((phi%np.pi))*1e-0
    analyzer = np.matmul(LinearPolarizer(theta),LinearRetarder(phi, delta))
    other = np.matmul(LinearRetarder(alpha, np.pi-0.2),LinearPolarizer(xi+np.pi/2))
    everything = np.matmul(analyzer,other)
    Sf = np.matmul(everything,S_input)
    return Sf[0]

def Simulation(S_input,system, args):
    """
    Simulation of what should we obtain in the angle-voltage curve for a specific incident
    light, S_input, going through a system
    """
    Angmax = 360
    Intensidad = np.zeros(Angmax)
    for i in range(Angmax):
        ang = i/180*np.pi
        Intensidad[i] = system(ang, S_input,args)
        
    plt.figure("Intensidad")
    plt.plot(np.arange(Angmax),Intensidad)
    plt.show()
    return Intensidad

if __name__ == "__main__":
    plt.close("all")
    # I1 = 1.
    # I2 = 1.
    # I3 = 1.
    # I0 = np.sqrt(I1**2+I2**2+I3**2)
    # # I0 = 1.
    # S_ini = np.array([I0,I1,I2,I3])
    # ploteo(S_ini)
    # Simulation(S_ini,SimpleSystem,(np.pi/2,np.pi/2,np.pi/2))
    
    ### Prueba de aumento teorico de elipticidad por lambda/medios --> No hay
    I1,I2,I3 = 1.,1.,0.2
    I0 = np.sqrt(I1**2+I2**2+I3**2)
    S_ini = np.array([I0,I1,I2,I3])
    ploteo(S_ini,title="Antes")
    Sf = np.matmul(LinearRetarder(np.pi,np.pi),S_ini)
    ploteo(Sf,title="Despues")

    ### Prueba de un cambio elipticidad por tener un lambda medios que su retarde depende del angulo
    I1,I2,I3 = 1.,1.,0.2
    I0 = np.sqrt(I1**2+I2**2+I3**2)
    S_ini = np.array([I0,I1,I2,I3])
    ploteo(S_ini,title="Entrada")
    Simulation(S_ini,ErraticSystem,(np.pi/2,np.pi/2,np.pi/2))
    
    ### Ejemplos de luz
    path = r"C:\Users\USUARIO\Desktop\TFGEnd\TFGImagenes"+"\\"
    S_a = np.array([1,1,0,0])
    ploteo(S_a,pathSave=path+"EjA.png",title="A")
    S_b = np.array([1,0,1,0])
    ploteo(S_b,pathSave=path+"EjB.png",title="B")
    S_c = np.array([1,0,0,1])
    ploteo(S_c,pathSave=path+"EjC.png",title="C")
    S_d = np.array([1,1/np.sqrt(3),1/np.sqrt(3),-1/np.sqrt(3)])
    ploteo(S_d,pathSave=path+"EjD.png",title="D")
    

