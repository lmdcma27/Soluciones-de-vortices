import numpy as np
import matplotlib.pyplot as plt
 
 
#vector de ecuacioes diferenciales
def v(p,f,A,A_0,g):   
   N = 1.0
   U = 2.0  
   e = 1.0
   k = 1.0
   w = np.array([g, e*e*f*f*A_0/k - A/p, e*e*f*f*(A - 1/(e*p))/k, e*e*N*N*(A - 1/(e*p))*(A - 1/(e*p))*f - e*e*N*N*A_0*A_0*f - U*f*(1-f*f)*(3*f*f - 1)/2 - g/p])
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
      z.append( float(sp[4]) + float(sp[1])/float(sp[0]) )       
      w.append( float(sp[2]) )  
      u.append(0)             
   plt.plot(x,y, 'b')
   plt.plot(x,z, 'k')
   plt.plot(x,w, 'y')
   plt.plot(x,u, 'r')
   plt.savefig("Grafica para 9950 iteraciones con un paso h=0,001")
   #plt.savefig("tercero")
   plt.show()
   file.close()
 
 
#Las condiciones iniciales obtenidas son
v_0 = np.array([1.0,0.009,-0.009,-0.01037585])

#escalar
#v_0 = np.array([1.0,0.0386,0.0386,-0.00235] )
#v_0 = np.array([1.0,0.0377,0.0378,-0.00250020] )
#v_0 = np.array([1.0,0.009,-0.009,-0.01037585] )



iterar1(v_0, 10.0, 9990)
 
 
#np.array( [-0.000051,-0.0000062,1.0,1.0] )
#v_0 = np.array( [-0.00013,-0.000054,1.0,1.0] )
#np.array([-0.000034189,-0.000038438,1.0,1.0] )
#[-1.54312499e-02 -4.94110446e+00  9.53758241e-05  8.73509708e-05]
