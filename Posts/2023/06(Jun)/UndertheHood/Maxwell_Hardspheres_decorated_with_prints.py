import numpy as np
import pdb

#global values
L       = 1              #meters
Natoms  = 200           
Ratom   = 0.03           #meters
T       = 300            #Kelvin (K)
k_b     = 1.4e-23        #J/K
dt      = 1e-5           #seconds
N_A     = 6.0221408e+23
He_mm   = 4e-3           #kilograms/mole
He_mass = He_mm/N_A      #kilograms

def initialize_atoms():
    poss     = L*(np.random.rand(Natoms,3)-0.5)
    vel_mag  = np.sqrt(3*k_b*T/He_mass)
    thetas   = np.pi*np.random.rand(Natoms)
    phis     = 2*np.pi*np.random.rand(Natoms)
    vels     = vel_mag*np.array([np.sin(thetas)*np.cos(phis),np.sin(thetas)*np.sin(phis),np.cos(thetas)]).T #verified by hand
    return poss, vels

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
                print(i,j,pos_i,pos_j,rel_pos,rel_pos.dot(rel_pos),d_sqr)
                hitlist.append((i,j))
    
    return hitlist

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
    mass    = He_mass/N_A  #mass of the atom in kg
    mass = 4E-3/6E23 #temp assignment for debugging
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
        print('Inside update_velocities_for_collision(poss,vels) - pre bounce:')
        print('i & j: ', i,j)
        print('poss[i] & poss[j]: \n',poss[i],'\n',poss[j],'\n',rrel)
        print('vels[i] & vels[j]: \n',vels[i],'\n',vels[j],'\n',vrel)
        #scattering geometry
        dx       = rrel.dot(vrel/vrel_mag)
        dy_vec   = np.cross(rrel,vrel/vrel_mag)
        dy       = np.sqrt(dy_vec.dot(dy_vec))
        alpha    = np.arcsin(dy/(2*Ratom))
        if dy > 2*Ratom:
            print('***************************')
            #pdb.set_trace()
            print(hitlist)
            print('i: ',i,'j: ',j)
            print('pos_i: ',pos_i)
            print('pos_j: ',pos_j)
            print('rrel: ',rrel)
            print('vel_i: ',vel_i)
            print('vel_j: ',vel_j)
            print('dy & 2*Ratom: ',dy, 2*Ratom)
        d        = (2*Ratom)*np.cos(alpha)-dx #distance traveled into the atom from first contact
        deltat   = d/vrel_mag                 #time spent moving from first contact to position inside atom
        
        #back up the particles until the first point of contact
        pos_i = pos_i - vel_i*deltat #back up particle i to contact configuration
        pos_j = pos_j - vel_j*deltat #back up particle j to contact configuration
        
        #check the relative distance at the contact (star) point
        rrel_star      = pos_i - pos_j
        rrel_mag_star  = np.sqrt(rrel_star.dot(rrel_star))        
        if np.abs(2*Ratom - rrel_mag_star) > 1e-7:
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
    #hitlist = update_velocities_for_collision(poss,vels)
    update_velocities_for_collision(poss,vels)
    update_velocities_for_bounce(poss,vels)
    #return hitlist

def update_state(poss,vels):
    print('Just inside update_state(poss,vels): \n',poss[0:2],'\n',vels[0:2],'\n')
    update_positions(poss,vels)
    print('Just after update_positions(poss,vels): \n',poss[0:2],'\n',vels[0:2],'\n')    
    #hitlist = update_velocities(poss,vels)
    update_velocities(poss,vels)
    print('Just after update_velocities(poss,vels): \n',poss[0:2],'\n',vels[0:2],'\n')
    #return hitlist

#############################################################################################
poss, vels = initialize_atoms()

#poss = np.array([
#                 [-0.276779, -0.180692, -0.282981],
#                 [-0.31255, -0.144077, -0.284417],
#                ])
#vels = np.array([
#                 [-302.906, -504.161, -1242.61],
#                 [218.142, -676.402, 1176.82]
#                ])

print('Just after initialize_atoms(): \n',poss[0:2],'\n',vels[0:2],'\n')
for i in range(100):
    #hitlist = update_state(poss,vels)
    print(i)
    update_state(poss,vels)
    #i,j = hitlist[0]
    #print('Just after update_state(poss,vels): \n',poss[i],'\n',poss[j],'\n',vels[i],'\n',vels[j],'\n')