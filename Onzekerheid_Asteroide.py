import numpy as np
import numpy.ctypeslib as ctl
import ctypes
import numpy
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import numpy.linalg as la


from Help_functions import *


libname="N_body.so"
libdir = './'
mylib=ctl.load_library(libname, libdir)




functie_botsing=mylib.perturbation_data_final
functie_botsing.restype = ctypes.c_longlong
functie_botsing.argtypes = [ctl.ndpointer(np.float64,flags='aligned, c_contiguous'),ctl.ndpointer(np.float64,flags='aligned, c_contiguous'),ctypes.c_int,ctypes.c_double, ctypes.c_bool, ctypes.c_bool]
functie_botsing_lyapunov=mylib.perturbation_data_lyapunov
functie_botsing_lyapunov.restype = ctypes.c_longlong
functie_botsing_lyapunov.argtypes = [ctl.ndpointer(np.float64,flags='aligned, c_contiguous'),ctl.ndpointer(np.float64,flags='aligned, c_contiguous'),ctypes.c_int,ctypes.c_double, ctypes.c_bool, ctypes.c_bool,ctl.ndpointer(np.float64,flags='aligned, c_contiguous'),ctl.ndpointer(np.float64,flags='aligned, c_contiguous')]


Time_iter=876


def get_Botsing_data(N_t, load_state, save_state):
    #Runs simulation for some fixed time T=30 days
    global End_day
    End_day=365
    global dt
    dt=End_day/N_t*24*60*60
    print("dt:",dt/(3600), "uur")
    # bevat x y en z coordinaat
    body_1=np.zeros(3*abs(N_t)+3,dtype=np.float64)
    body_2=np.zeros(3*abs(N_t)+3,dtype=np.float64)
    functie_botsing(body_1, body_2, abs(N_t), dt,load_state, save_state)
    return body_1, body_2


def get_Botsing_data_lyapunov(N_t, load_state, save_state):
    #Runs simulation for some fixed time T=30 days
    global End_day
    End_day=180
    global dt
    dt=End_day/N_t*24*60*60
    # N_t+=1
    # bevat x y en z coordinaat
    body_1=np.zeros(3*abs(N_t),dtype=np.float64)
    body_2=np.zeros(3*abs(N_t),dtype=np.float64)
    sep=-1000
    # sep=0.1
    init_diff_vector=np.array([sep,sep,sep,sep,sep,sep,sep,sep,sep,sep,sep,sep],dtype=np.float64)
    diff_vector=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],dtype=np.float64)
    functie_botsing_lyapunov(body_1, body_2, abs(N_t), dt,load_state, save_state,init_diff_vector,diff_vector)
    print(diff_vector)
    print("Python dt:", dt,"s", dt/(24*3600),"days")


    #Calculate lyapunov exponent
    T=End_day*24*60*60


    component=0
    lambda_1=(np.log(la.norm(diff_vector[component]))-np.log(la.norm(init_diff_vector[component])))/T
    print("lambda_1", lambda_1)
    return body_1, body_2





terugreken=True
# terugreken=False


if terugreken:
    bool_1=False
    bool_2=True
    Time_iter=-Time_iter
else:
    bool_1=True
    bool_2=False
# bool_1=False
# bool_2=False






import time


st = time.time()
x_pos_sun, x_pos_asteroid=get_Botsing_data(N_t=Time_iter, load_state=bool_1, save_state=bool_2)
# print(x_pos_sun[0::3])
# get the execution time
et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_pos_sun[0::3],x_pos_sun[1::3],x_pos_sun[2::3],label="earth")
ax.plot(x_pos_asteroid[0::3],x_pos_asteroid[1::3],x_pos_asteroid[2::3],label="asteroide")


a=1.5
sample=0
ax.set(xlim3d=(-1.496e11*a, a*1.496e11), xlabel='X')
ax.set(ylim3d=(-1.496e11*a, a* 1.496e11), ylabel='Y')
ax.set(zlim3d=(-1.496e11*a, a*1.496e11), zlabel='Z')
ax.view_init(elev=80, azim=-90)
plt.legend()
file_name="Asteroid_inconsistent_trajectory"+"_dt_"+ "%s" % str(dt/(24*3600)) +".eps"
# plt.savefig(file_name)
# plt.savefig("Asteroid_inconsistent_trajectory.png")


plt.show()












plt.figure()
# plt.plot(x_pos_sun[0::3],x_pos_sun[1::3],label="earth")
plt.plot(x_pos_asteroid[0::3],x_pos_asteroid[1::3],label="asteroide")
plt.legend()
a=a*10
plt.xlim([-1.496e11*a, a*1.496e11])
plt.ylim([-1.496e11*a, a* 1.496e11])
plt.show()




# max_arr=[]
# for h in range(0,3):
#     x_pos_sun,x_pos_asteroid=get_Botsing_data(N_t=Time_iter*2**h,load_state=bool_1,save_state=bool_2)
#     max_arr.append(x_pos_asteroid[0::3])
#
#
#
#
# print("test", np.arange(10)[0::2])
#
#
#
#
# calculate_richardson(max_arr)














t=0
def update_lines(num, walks, lines):
    global t
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[num-1:num, :2].T)
        line.set_3d_properties(walk[num-1:num, 2])
        plt.title("Days %s" % str(t/(2*24*3600)))
        t+=dt




    return lines




# Data: 40 random walks as (num_steps, 3) arrays


sample=1
x=x_pos_sun[0::3][0::sample]
y=x_pos_sun[1::3][0::sample]
z=x_pos_sun[2::3][0::sample]
num_steps = len(x)
walk = np.array(list(zip(x,y,z)))




x=x_pos_asteroid[0::3][0::sample]
y=x_pos_asteroid[1::3][0::sample]
z=x_pos_asteroid[2::3][0::sample]
walk_1 = np.array(list(zip(x,y,z)))




walks = np.array([walk, walk_1])


# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")


# Create lines initially without data
String_names_list=["earth","astroide"]
lines = [ax.plot([], [], [], "o",label=str(String_names_list[count]))[0] for count,_ in enumerate(walks)]


a=1
sample=-1
view_earth=0
ax.set(xlim3d=(-1.496e11*a+x_pos_sun[0::3][0]*view_earth,x_pos_sun[0::3][0]*view_earth+ a*1.496e11), xlabel='X')
ax.set(ylim3d=(-1.496e11*a+x_pos_sun[1::3][0]*view_earth,x_pos_sun[1::3][0]*view_earth+ a* 1.496e11), ylabel='Y')
ax.set(zlim3d=(-1.496e11*a+x_pos_sun[2::3][0]*view_earth,x_pos_sun[2::3][0]*view_earth+ a*1.496e11), zlabel='Z')
ax.view_init(elev=20, azim=-90)
ax.legend()


# Creating the Animation object
ani = animation.FuncAnimation(fig, update_lines, num_steps, fargs=(walks, lines), interval=1, repeat=False)


plt.show()




# print("debug", np.max(x_pos_asteroid[1::3]))
