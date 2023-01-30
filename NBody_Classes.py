import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
print(mpl.get_backend())
mpl.use("TkAgg")
import matplotlib.animation as animation
np.random.seed(2384724)

##################################################################################
                            ### EXCESS FUNCTIONS ###

def update_lines_w(num, walks, lines):
    '''Function used for plotting body with path history'''
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[:num, :2].T)
        line.set_3d_properties(walk[:num, 2])
    return lines

def update_lines_wo(num, walks, lines):
    '''Function used for plotting body without path history'''
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[num-1:num, :2].T)
        line.set_3d_properties(walk[num-1:num, 2])
    return lines

def calculate_richardson(max_arr):
    print("lengtes", len(max_arr[0]))
    print("lengtes", len(max_arr[1][0::2]))
    print("lengtes", len(max_arr[2][0::4]))

    max_arr[1] = max_arr[1][0::2]
    max_arr[2] = max_arr[2][0::4]
    rich_arr = (max_arr[1] - max_arr[0]) / (max_arr[2] - max_arr[1])

    # we houden alleen de positieve waarden
    rich_arr = rich_arr[rich_arr > 0]
    print("richardson error:", rich_arr, "fout orde:", np.log(np.nanmean(rich_arr)) / np.log(2))

    # if (len(max_arr[1][1::2])>len(max_arr[0])):
    #     max_arr[1]=np.delete(max_arr[1][1::2],-1)
    # else:
    #     max_arr[1]=max_arr[1][1::2]
    # if (len(max_arr[2][3::4])>len(max_arr[0])):
    #     max_arr[2]=np.delete(max_arr[2][3::4],-1)
    # else:
    #     max_arr[2]=max_arr[2][3::4]
    #
    #
    # rich_arr=(max_arr[1]-max_arr[0])/ (max_arr[2]-max_arr[1])
    # print("richardson:",np.nanmean(rich_arr), "fout orde:", np.log(np.nanmean(rich_arr))/np.log(2))

# def generate_random_pos():

##################################################################################
                            ### BODY CLASS ###

class Body(object):

    def __init__(self, name, mass, radius, pos_vec, vel_vec):
        self.name = name # a string
        self.mass = mass
        self.rad = radius
        self.pos = pos_vec
        self.vel = vel_vec
        self.acc = np.zeros(3)
        self.x_hist = np.array([self.pos[0]])
        self.y_hist = np.array([self.pos[1]])
        self.z_hist = np.array([self.pos[2]])
        self.in_vel = vel_vec
        #self.vx_hist = np.array([self.vel[0]])
        #self.vy_hist = np.array([self.vel[1]])
        #self.vz_hist = np.array([self.vel[2]])
        self.time_hist = np.array(0)
        self.walk = None
        self.hit_status = ''    # Geeft aan of er een hit heeft plaatsgevonden met iemand anders

    def r_vect(self, other):
        return other.pos - self.pos

    def r_norm(self, other):
        return np.linalg.norm(self.r_vect(other))

    def hit(self, other):
        if self.r_norm(other) <= (self.rad + other.rad):
            self.hit_status = other.name
            return True
        return False

    def save_path_file(self):
        # - seed
        # - path (t,x,y,z)
        # - hit (True, False)
        hit_reg = np.zeros(len(self.z_hist))
        hit_reg[0] = self.hit_status
        log = pd.dataframe([self.time_hist, self.x_hist, self.y_hist, self.z_hist, hit_reg], index=['t', 'x', 'y', 'z', 'Hit'])
        # log.set_axis(self.time_hist, axis='columns', inplace=True)
        # log.to_excel("Halal gaming")
        return log

    def call_path_file(self):
        return 6.9
##################################################################################
                            ### SYSTEM CLASS ###

class Body_System(object):

    def __init__(self, body_list, time_step):
        self.bodies = body_list
        self.n = len(body_list)
        self.dt = time_step
        self.walks = []

    def grav_const(self):
        return 6.67408e-11

    def au(self):
        return 1.4717e11

    def update_force(self): # Kracht berekenen op elke body
        for body in self.bodies:
            body.acc = np.zeros(3)
        for i in range(self.n):
            body_i = self.bodies[i]
            for j in range(i+1, self.n):
                body_j = self.bodies[j]
                force_i_j = self.grav_const() * body_i.mass * body_j.mass * body_i.r_vect(body_j) / (body_i.r_norm(body_j)**3)
                body_i.acc += force_i_j/body_i.mass
                body_j.acc -= force_i_j/body_j.mass

    def initialize_vel(self):  # Snelheid updaten
        self.update_force()
        for body in self.bodies:
            body.vel = body.vel + 0.5 * body.acc * self.dt

    def finalize_vel(self):
        #print("new step")
        for body in self.bodies:
            body.pos += body.vel * self.dt
            #print("{}: {}".format(body.name, body.pos))

        self.update_force()

        for body in self.bodies:
            body.vel += 0.5 * body.acc * self.dt

    def check_hit(self):
        for i in range(self.n):
            for j in range(i+1, self.n):
                if self.bodies[i].hit(self.bodies[j]):
                    print("{} hit {}".format(self.bodies[i].name, self.bodies[j].name))
                    return True
        return False

    def set_newstart(self):
        for body in self.bodies:
            #body.pos = np.array([body.x_hist[-1], body.y_hist[-1], body.z_hist[-1]])
            body.x_hist = np.array(body.pos[0])
            body.y_hist = np.array(body.pos[1])
            body.z_hist = np.array(body.pos[2])
            body.in_vel = body.vel
            body.time_hist = np.array(0)
            body.walk = None
            body.hit_status = ''

    def refresh(self):
        for body in self.bodies:
            body.pos = np.array([body.x_hist[0], body.y_hist[0], body.z_hist[0]])
            body.vel = body.in_vel
            body.acc = np.zeros(3)
            body.x_hist = np.array(body.pos[0])
            body.y_hist = np.array(body.pos[1])
            body.z_hist = np.array(body.pos[2])
            body.time_hist = np.array(0)
            body.walk = None
            body.hit_status = ''

    def initialize_asteroid_pos(self, a, asteroid, earth, sun):
        r = sun.r_vect(earth)
        r_unit = r / sun.r_norm(earth)
        w_unit = np.array([-1*np.sign(earth.vel[0])*r_unit[1], -1*np.sign(earth.vel[1])*r_unit[0], earth.pos[2]])

        asteroid.pos = sun.pos + r + a * w_unit
        asteroid.x_hist = np.array([asteroid.pos[0]])
        asteroid.y_hist = np.array([asteroid.pos[1]])
        asteroid.z_hist = np.array([asteroid.pos[2]])
        asteroid.vel = -1 * np.linalg.norm(asteroid.vel) * w_unit

    def generate_asteroid_hit(self, T, asteroid, earth, sun):

        for i in range(1,40):
            self.initialize_asteroid_pos(i*0.01*self.au(), asteroid, earth, sun)
            hit = False
            t = 0
            self.initialize_vel()

            while not (hit or abs(t) >= T):
                self.Leapfrog_step()
                hit = self.check_hit()
                self.update_pos_hist()
                t += self.dt  # tijdstap variabel maken

            self.set_walks()
            self.plot_path()
            self.distance_plot(t-self.dt, [self.bodies[0], self.bodies[1]])

            print("initial vel: ")
            print("vel: ", self.bodies[0].vel)

            if hit:
                break
            else:
                self.refresh()

    def pos_perturb_sim(self, T, body, r_onz, N):
        self.dt = -self.dt
        self.Leapfrog(T)
        # print("before")
        # print("pos: ", body.pos)
        # print("vel: ", body.vel)
        # print("acc: ", body.acc)
        # print("xhist: ", body.x_hist)
        # print("in vel: ", body.in_vel)
        self.set_newstart()
        hit_counter = 0
        self.dt = -self.dt
        # print("after")
        # print("pos: ", body.pos)
        # print("vel: ", body.vel)
        # print("acc: ", body.acc)
        # print("xhist: ", body.x_hist)
        # print("in vel: ", body.in_vel)

        perturbations = body.pos + np.random.normal(0, r_onz, (1, N, 3))[0]
        for i in range(len(perturbations)):
            print("Start position {}: {}".format(i, perturbations[i]))
            body.pos = perturbations[i]

            self.Leapfrog(T+2*abs(self.dt))
            self.set_walks()
            #self.distance_plot(T+2*abs(self.dt), [self.bodies[0], self.bodies[1]])

            if body.hit_status == 'Earth':
                hit_counter += 1
            print("hit status: ", body.hit_status)
            self.refresh()


        print("amount of hits: ", hit_counter)
        return hit_counter

##################################################################################
                            ### PLOT DEVICES ###

    def update_pos_hist(self):  # Halal positions appenden in arrays
        for body in self.bodies:
            body.x_hist = np.append(body.x_hist, body.pos[0])
            body.y_hist = np.append(body.y_hist, body.pos[1])
            body.z_hist = np.append(body.z_hist, body.pos[2])

    def update_time_hist(self, t):
        for body in self.bodies:
            body.time_hist = np.append(body.time_hist, t)


    def plot_path(self):
        fig = plt.figure()
        ax = plt.axes(projection = "3d")
        for body in self.bodies[:2]:
            ax.plot(body.x_hist, body.y_hist, body.z_hist)

        # a = 1.5
        # ax.set(xlim3d=(-1.496e11*a, a*1.496e11), xlabel='X')
        # ax.set(ylim3d=(-1.496e11*a, a*1.496e11), ylabel='Y')
        # ax.set(zlim3d=(-1.496e11*a, a*1.496e11), zlabel='Z')
        ax.set(xlabel='X')
        ax.set(ylabel='Y')
        ax.set(zlabel='Z')

        plt.show()

    def set_walks(self):
        for body in self.bodies:
            body.walk = np.array(list(zip(body.x_hist, body.y_hist, body.z_hist)))
            self.walks.append(body.walk)

    def animate(self, T, path=True):
        num_steps = int(np.floor(T/self.dt))

        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        lines = [ax.plot([], [], [])[0] for _ in self.walks]

        ax.set(xlim3d=(-2e11, 2e11), xlabel='X')
        ax.set(ylim3d=(-2e11, 2e11), ylabel='Y')
        ax.set(zlim3d=(-2e11, 2e11), zlabel='Z')

        if path:
            ani = animation.FuncAnimation(fig, update_lines_w, num_steps, fargs=(self.walks, lines), interval=1)
        else:
            ani = animation.FuncAnimation(fig, update_lines_wo, num_steps, fargs=(self.walks, lines), interval=1)

        plt.show()

    def distance_plot(self, T, body_names):
        #t_array = np.arange(0, T+2*self.dt, self.dt)
        t_array = body_names[0].time_hist
        r_norm_array = np.array([])
        for i in range(len(body_names[0].walk)):
            r = np.linalg.norm(body_names[0].walk[i]-body_names[1].walk[i])
            r_norm_array = np.append(r_norm_array, r)

        print(min(r_norm_array))
        r_norm_array = r_norm_array # /self.au()

        plt.plot(t_array, r_norm_array)
        plt.show()

    # def save_path_files(self):      "Functie voor halal dataframe"
        # df = pd.concat([i.save_path_file for i in self.bodies])
        # df.to_excel("halalgamingpc.xlsx")


##################################################################################
                            ### INTEGRATORS ###


    def Leapfrog_step(self):

        #print("new step")
        for body in self.bodies:
            body.pos += body.vel * self.dt
            #print("{}: {}".format(body.name, body.pos))

        self.update_force()

        for body in self.bodies:
            body.vel = body.vel + body.acc * self.dt


    def Leapfrog(self, T):
        ''''Function that ...'''
        t = 0
        hit = False

        self.initialize_vel()

        while abs(t) < T-self.dt:
            self.Leapfrog_step()
            hit = self.check_hit()
            self.update_pos_hist()
            t += self.dt
            self.update_time_hist(t)
            if hit:
                #print("time: ", t)
                break

        if not hit:
            self.finalize_vel()
            self.check_hit()
            self.update_pos_hist()
            t += self.dt
            self.update_time_hist(t)

        self.set_walks()
        self.plot_path()
        #self.animate(T)

    def simulation(self, T, integrator):
        self.dt = -self.dt
        integrator(T)
        self.dt = -self.dt
        integrator(T)
        self.dt = 3600
        integrator(86400*31)

import time
def main():
    # get the start time
    st = time.time()
    AU = 1.4717e11
    Nt=876
    T = 86400*20  # seconds
    dt_rich = -T/Nt
    print("dt_rich: ", dt_rich/60, "days")
    # dt_rich_2 = 3600
    Earth = Body('Earth', 5.9722e24, 6.371e6, np.array([0, AU, 0], dtype=np.float64), np.array([-3e4, 0, 0], dtype=np.float64))    # m in kg, al die andere dingen in m(/s)
    # Moon = Body('Moon', 7.34767309e22, 5e6, np.array([0,3.844e8+1.4717e11, 0], dtype=np.float64), np.array([3683/3.6+3e4, 0, 0], dtype=np.float64))
    Sun = Body('Sun', 1.989e30, 0, np.array([0, 0, 0], dtype=np.float64), np.array([0, 0, 0], dtype=np.float64))
    Jupiter = Body('Jupiter', 1.898e27, 0, np.array([7.4076e11,0,0], dtype=np.float64), np.array([0,1.307e4,0], dtype=np.float64))
    Asteroid = Body('Bennu', 7.8e10, 262.5, np.array([6.371e6+262, AU, 0], dtype=np.float64), np.array([-35e3, 0, 0], dtype=np.float64)) # Hoek van 5 graden op baan van de aarde
    Universe = Body_System([Asteroid, Earth, Sun, Jupiter], dt_rich)    # dt in seconden
    # Earth_1 = Body('Earth', 5.9722e24, 6.371e6, np.array([0, AU, 0], dtype=np.float64), np.array([-3e4, 0, 0], dtype=np.float64))  # m in kg, al die andere dingen in m(/s)
    # Sun_1 = Body('Sun', 1.989e30, 0, np.array([0, 0, 0], dtype=np.float64), np.array([0, 0, 0], dtype=np.float64))
    # Jupiter_1 = Body('Jupiter', 1.898e27, 0, np.array([7.4076e11, 0, 0], dtype=np.float64), np.array([0, 1.307e4, 0], dtype=np.float64))
    # Asteroid_1 = Body('Bennu', 7.8e10, 262.5, np.array([6.371e6+262, AU, 0], dtype=np.float64), np.array([-32e3, 0, 0], dtype=np.float64))  # Hoek van 5 graden op baan van de aarde
    # Universe_1 = Body_System([Asteroid_1, Earth_1, Sun_1, Jupiter_1], dt_rich_2)  # dt in seconden


    #Universe.Leapfrog(T)
    #Universe.distance_plot(T, [Earth, Asteroid])
    #Universe.simulation(T, Universe.Leapfrog)
    #Universe.animate(T, path=True)
    #Universe.generate_asteroid_hit(T, Asteroid, Earth, Sun)
    #Universe.pos_perturb_sim(T, Asteroid, 3000, 100)
    Universe.initialize_asteroid_pos(0.1 * AU, Asteroid, Earth, Sun)
    Universe.Leapfrog(T)
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
    # Universe_1.initialize_asteroid_pos(0.1 * AU, Asteroid_1, Earth_1, Sun_1)
    # Universe_1.Leapfrog(T)
    print(len(Asteroid.x_hist))
    # print(len(Asteroid_1.x_hist))
    # boobies = (Asteroid.x_hist[1::12]-Asteroid_1.x_hist[1:])/Asteroid_1.x_hist[1:] * 100
    # print(boobies)
    # print(np.mean(boobies))

if __name__ == '__main__':
    main()
