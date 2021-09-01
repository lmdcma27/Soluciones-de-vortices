
import numpy as np
import matplotlib.pyplot as plt


#vector de ecuacioes diferenciales
def v(p,phi1,a1_0,phi,a_0):    
    N = 1.0
    U = 2.0   
    w = np.array([((N-a_0)*(N-a_0)*phi)/(p*p) - U*(1-phi*phi)*phi/2 - phi1/p, a1_0/p - (N - a_0)*phi*phi, phi1, a1_0 ])
    return w

def iterar2(v_0,p,n_i):
#abrimos y archivo para guardar los datos        
    #definimos el paso
    h = 0.001    
    #iteramos con el metodo de RK4
    for i in range(n_i):          
        #construimos los vectores k1, k2, k3 y k4        
        k_1 = h*v(p, v_0[0],v_0[1],v_0[2],v_0[3])        
        k_2 = h*v(p - h/2,v_0[0] + k_1[0]/2,v_0[1] + k_1[1]/2 ,v_0[2] + k_1[2]/2,v_0[3] + k_1[3]/2)
        k_3 = h*v(p - h/2,v_0[0] + k_2[0]/2,v_0[1] + k_2[1]/2 ,v_0[2] + k_2[2]/2,v_0[3] + k_2[3]/2)
        k_4 = h*v(p - h,v_0[0] + k_3[0],v_0[1] + k_3[1], v_0[2] + k_3[2],v_0[3] + k_3[3])
        #calculamos la siguiente iteración
        v_1 = v_0 + (k_1 + 2*k_2 + 2*k_3 + k_4)/6
        v_0 = v_1     
        p = p - h                             
        #hacemos que p tome el valor del siguinte paso                                         
    return v_0
  
def iterar1(v_0,p, n_i):
#abrimos y archivo para guardar los datos 
    file = open("datos2.txt", "w")
    #guardamos el primer dato
    file.write( str(p) + " " + str(v_0[0]) + " " + str(v_0[1]) + " " +str(v_0[2]) +" "+ str(v_0[3]) + '\n')      
    h = 0.001    
    #iteramos
    for i in range(n_i):          
        #construimos los vectores k1, k2, k3 y k4        
        k_1 = h*v(p, v_0[0],v_0[1],v_0[2],v_0[3])        
        k_2 = h*v(p - h/2, v_0[0] + k_1[0]/2,v_0[1] + k_1[1]/2 ,v_0[2] + k_1[2]/2,v_0[3] + k_1[3]/2)
        k_3 = h*v(p - h/2, v_0[0] + k_2[0]/2,v_0[1] + k_2[1]/2 ,v_0[2] + k_2[2]/2,v_0[3] + k_2[3]/2)
        k_4 = h*v(p - h, v_0[0] + k_3[0],v_0[1] + k_3[1], v_0[2] + k_3[2],v_0[3] + k_3[3])
        #calculamos la siguiente iteración
        v_1 = v_0 + (k_1 + 2*k_2 + 2*k_3 + k_4)/6
        v_0 = v_1     
        p = p - h        
        #hacemos que r tome el valor del siguinte paso         
        file.write( str(p) + " " + str(v_0[0])+ " " +str(v_0[1])+ " "+ str(v_0[2])+ " "+ str(v_0[3])+ '\n')               
    print(v_0)
    #cerramos el achirvo     
    file.close()        
    graficar()


def graficar():
    file = open("datos2.txt", "r")
    x = []
    y = []
    z = []
    w = []
    u = []
    for line in file:        
        sp = line.split(" ")
        x.append( float(sp[0]) )
        y.append( float(sp[1]) )        
        z.append( float(sp[2]) )
        w.append( float(sp[3]) )        
        u.append( float(sp[4]) )        
    plt.plot(x,w)
    plt.plot(x,u)
    plt.savefig("Grafica para 9950 iteraciones con un paso h=0,001")
    plt.show()
    file.close()

#Implementación del shooting Method
def shooting(p, l, vr):

    #Llamamos al método numérico de cuarto orden
    #y le pasamos los valores necesarios        
    n_i = 9950
    m = []
    z = []
    #número de variables
    #la formula para la longitud puede variar, dependiendo
    # del número de segundas derivada     
    lon = len(l[0])//2    
    #iteramos para obtener los valores necesarios para
    #hacer una iterpolación
    for i in range(lon):                
        m.append(iterar2(l[i], p,n_i)[:lon] )                                        
    #Ahora interpolamos los valores más adecuados para las derivadas
    for j in range(lon):
        #Primero calculamos la pendiente            
        M = (m[j][1] - m[j][0])/(l[j][1] - l[j][0])        
        
       #guardamos los en la lista z
        z.append( m[j][0] + M*(vr[j] - l[j][0]) )      
    print(z)                  
#shooting(10, np.array([[-0.0000257,-0.000037,1.0,1.0], [-0.00013,-0.000054,1.0,1.0] ]), np.array([0,0]) )

#Las condiciones iniciales obtenidas son

v_0 = np.array([-0.000034189,-0.000038438,1.0,1.0] )
shooting(10, np.array([[-0.00000000257,-0.00000000037,1.0,1.0], [-0.000000013,-0.0000000054,1.0,1.0] ]), np.array([0,0]) )

iterar1(v_0, 10.0, 9950)


#np.array( [-0.000051,-0.0000062,1.0,1.0] )
#v_0 = np.array( [-0.00013,-0.000054,1.0,1.0] )
#np.array([-0.000034189,-0.000038438,1.0,1.0] )
#[-1.54312499e-02 -4.94110446e+00  9.53758241e-05  8.73509708e-05]
