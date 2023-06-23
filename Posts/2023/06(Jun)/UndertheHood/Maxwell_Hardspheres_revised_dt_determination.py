import numpy as np
import pdb
import pickle

#global values
L       = 1              #meters
Natoms  = 400           
Ratom   = 0.03           #meters
T       = 300            #Kelvin (K)
k_b     = 1.4e-23        #J/K
dt      = 1e-5           #seconds
N_A     = 6.0221408e+23
He_mm   = 4e-3           #kilograms/mole
He_mass = He_mm/N_A      #kilograms

def update_positions(poss,vels):
    poss[...] = poss[...] + vels[...]*dt


def check_for_collisions(poss):
    hitlist = []
    d_sqr   = 4*Ratom*Ratom
    for i in range(len(poss)):
        pos_i = poss[i]
        for j in range(i+1,Natoms):
            pos_j   = poss[j]
            rel_pos = pos_i - pos_j
            if rel_pos.dot(rel_pos) < d_sqr: 
                hitlist.append((i,j))
    
    return hitlist

def norm(a):
    return np.sqrt(a.dot(a))

def determine_dt(r1,r2,v1,v2):
    #construct the relative geometry
    rrel     = r1 - r2
    rrel_mag = norm(rrel)
    rrel_hat = rrel/rrel_mag
    vrel     = v1 - v2
    vrel_mag = norm(vrel)
    vrel_hat = vrel/vrel_mag
    hrel     = np.cross(rrel,vrel)
    z_hat    = hrel/norm(hrel)
    y_hat    = np.cross(z_hat,vrel_hat)
    
    #contruct the collision triangle solution for d
    r_ctheta  = rrel_mag*rrel_hat.dot(vrel_hat)
    radical   = np.sqrt( r_ctheta*r_ctheta - rrel_mag*rrel_mag + 4*Ratom*Ratom )
    dp        = r_ctheta + radical
    dm        = r_ctheta - radical
    
    dtp = dp/vrel_mag
    dtm = dm/vrel_mag
    
    return dtp, dtm

def reflect(pos,vel,n):
    #hand verified by directly manipulating poss
    if np.abs(pos[n]) > L/2:
        ddist  = np.abs(pos[n]) - L/2
        pos[n] = np.sign(pos[n])*(L/2 - ddist)
        vel[n] = -1*vel[n]
        
    return pos,vel

def update_velocities_for_bounce(poss,vels):
    for pos,vel in zip(poss,vels):
        pos,vel = reflect(pos,vel,0)
        pos,vel = reflect(pos,vel,1)
        pos,vel = reflect(pos,vel,2)

def update_velocities_for_collision(poss,vels):
    hitlist = check_for_collisions(poss)
    mass    = He_mass
    for i,j in hitlist:
        #get the individual states
        pos_i = poss[i]
        pos_j = poss[j]
        vel_i = vels[i]
        vel_j = vels[j]
        
        #calculate the relative position and velocity
        rrel     = pos_i - pos_j #the way Sherwood has it ???
        vrel     = vel_j - vel_i #this one is opposite ??
        rrel_mag = np.sqrt(rrel.dot(rrel))        
        vrel_mag = np.sqrt(vrel.dot(vrel))
        rrel_hat = rrel/rrel_mag
        vrel_hat = vrel/vrel_mag        

        #check to see that the overlap persists and was not
        #change by an earlier bounce this time step
        if rrel_mag > 2*Ratom: continue

        #scattering geometry
        dta, dtb = determine_dt(pos_i,pos_j,vel_i,vel_j)
        if max(dta,dtb) > 0:
            deltat = max(dta,dtb)
        else:
            print('deltat is negative - help!')
                
        #back up the particles until the first point of contact
        pos_i = pos_i - vel_i*deltat #back up particle i to contact configuration
        pos_j = pos_j - vel_j*deltat #back up particle j to contact configuration
        
        #check the relative distance at the contact (star) point
        rrel_star      = pos_i - pos_j
        rrel_star_mag  = np.sqrt(rrel_star.dot(rrel_star))        
        if np.abs(2*Ratom - rrel_star_mag) > 1e-7:
            print('still penetrating')
        
        #now adjust the momenta
        mtot   = 2*mass
        p_i    = mass*vel_i
        p_j    = mass*vel_j
        p_tot  = p_i + p_j
        #transform momenta to cm frame
        p_i_cm = p_i - p_tot*mass/mtot 
        p_j_cm = p_j - p_tot*mass/mtot
        #bounce in cm frame
        p_i_cm = p_i_cm - 2*p_i_cm.dot(rrel_hat)*(rrel_hat)
        p_j_cm = p_j_cm - 2*p_j_cm.dot(rrel_hat)*(rrel_hat)
        # transform momenta back to lab frame
        p_i    = p_i_cm + p_tot*mass/mtot 
        p_j    = p_j_cm + p_tot*mass/mtot 
        # move forward deltat in time
        pos_i  = pos_i + (p_i/mass)*deltat
        pos_j  = pos_j + (p_j/mass)*deltat
        
        #update positions and velocities
        poss[i] = pos_i
        poss[j] = pos_j
        vels[i] = p_i/mass
        vels[j] = p_j/mass
        #return hitlist

def update_velocities(poss,vels):
    update_velocities_for_collision(poss,vels)
    update_velocities_for_bounce(poss,vels)

def update_state(poss,vels):
    update_positions(poss,vels)
    update_velocities(poss,vels)

def MB_dist(vs):
    coeff = np.power((He_mass/2.0/np.pi/k_b/T),1.5)*4*np.pi*vs**2
    body  = np.exp(-He_mass*vs**2/2.0/k_b/T)
    
    return coeff*body    

import matplotlib.pyplot as plt
vs = np.arange(0,3000,10)
mb = MB_dist(vs)

num_time_steps = 2001
speed_configs  = np.zeros((num_time_steps,Natoms))
for i in range(num_time_steps):
    if i == 0:
        data_filename = 'c:/Users/byecs/Downloads/temp/sims/initialized_atoms.pkl'
        data_file     = open(data_filename,'rb')
        sim_data      = pickle.load(data_file)
        print(sim_data['name'])
        data_file.close()
        poss = sim_data['positions']
        vels = sim_data['velocities']
        print(poss[0,:])
    else:
        update_state(poss,vels)
    frame_number       = i
    frame_string       = f'{frame_number:05}.png'
    instant_speeds     = np.sqrt(vels[:,0]**2 + vels[:,1]**2 + vels[:,2]**2)
    speed_configs[i,:] = instant_speeds
    if i % 20 == 0:
        print(i,np.sum(instant_speeds/1e8))
        plt.hist(instant_speeds,bins=np.arange(0,3000,100),density=True);
        plt.plot(vs,mb,'r-')
        plt.ylim([0,0.0012])
        plt.xlabel('Atom Speed (m/s)',fontsize=12)
        plt.ylabel('Probability of Occurrence',fontsize=12)
        plt.savefig(frame_string)
        plt.cla()

sim_meta_data = '''
  Generated by Maxwell_Hardspheres_revised_dt_determination.py.  This simulation of Maxwell hard 
  spheres uses the original bounce algorithm of Bruce Sherwood but a revised quadratic deltat
  calculation.
'''
sim_name = 'MHS_revised_dt'

data = {'raw_data':speed_configs,'meta_data':sim_meta_data,'name':sim_name}
fileObj = open('MaxwellHardSphere_revised_dt_data.pkl', 'wb')
pickle.dump(data,fileObj)
fileObj.close()