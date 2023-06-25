import matplotlib        as mpl
import matplotlib.pyplot as plt
import numpy             as np
import pdb

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

def calculate_post_bounce_conditions_Sherwood(pos_i,pos_j,vel_i,vel_j):
    #calculate the relative position
    rrel_star     = pos_i - pos_j
    rrel_star_hat = rrel_star/norm(rrel_star)

    #now adjust the momenta
    mtot   = 2*mass
    p_i    = mass*vel_i
    p_j    = mass*vel_j
    p_tot  = p_i + p_j
    #transform momenta to cm frame
    p_i_cm = p_i - p_tot*mass/mtot 
    p_j_cm = p_j - p_tot*mass/mtot
    #bounce in cm frame
    p_i_cm = p_i_cm - 2*p_i_cm.dot(rrel_star_hat)*(rrel_star_hat)
    p_j_cm = p_j_cm - 2*p_j_cm.dot(rrel_star_hat)*(rrel_star_hat)
    # transform momenta back to lab frame
    p_i    = p_i_cm + p_tot*mass/mtot 
    p_j    = p_j_cm + p_tot*mass/mtot 
    return p_i/mass, p_j/mass
        

Ratom   = 0.03
dt      = 1.25e-5
N_A     = 6.0221408e+23
He_mm   = 4e-3           #kilograms/mole
He_mass = He_mm/N_A      #kilograms

#case1 - initial conditions at overlap time for headon collision
r1_c1   = np.array([0.02,0.07,0])
r2_c1   = np.array([0.04,0.02,0])
v1_c1   = np.array([800,1100,0])
v2_c1   = np.array([-400,-1200,0])

#case2 - initial conditions at overlap time for sideswipe collision
r1_c2   = np.array([0.02,0.07,0])
r2_c2   = np.array([0.03,0.04,0])
v1_c2   = np.array([200,900,0])
v2_c2   = np.array([300,1350,0])

#case 3 - initial conditions at the overlap time for a tbone collision
r1_c3   = np.array([0.02,0.07,0])
r2_c3   = np.array([0.03,0.04,0])
v1_c3   = np.array([200,900,0])
v2_c3   = np.array([1350,300,0])

dtp, dtm = determine_dt(r1_c1,r2_c1,v1_c1,v2_c1)
r1_c1_star = r1_c1 - v1_c1*dtp
r2_c1_star = r2_c1 - v2_c1*dtp
rrel_c1_star = r1_c1_star - r2_c1_star

dtp, dtm = determine_dt(r1_c2,r2_c2,v1_c2,v2_c2)
r1_c2_star = r1_c2 - v1_c2*dtp
r2_c2_star = r2_c2 - v2_c2*dtp
rrel_c2_star = r1_c2_star - r2_c2_star

dtp, dtm = determine_dt(r1_c3,r2_c3,v1_c3,v2_c3)
r1_c3_star = r1_c3 - v1_c3*dtp
r2_c3_star = r2_c3 - v2_c3*dtp
rrel_c3_star = r1_c3_star - r2_c3_star

#initialize
fig, ax = plt.subplots()
collision_flag = 0
bc_r0          = r1_c1
gc_r0          = r2_c1
bc_v           = v1_c1
gc_v           = v2_c1
mass           = He_mass
dt             = 1.25e-5
#set current time in the past to get ICs
curr_time      = -20*dt
bc             = bc_r0 + curr_time*bc_v
gc             = gc_r0 + curr_time*gc_v
frame_counter  = 0
#reset current time for cosmetic readout
curr_time      = 0

#propagate forward
for i in range(35):
    bc = bc + dt*bc_v
    gc = gc + dt*gc_v 
    #print(i,f'{curr_time*1e6:.2f}',bc,gc)   
    grey   = plt.Circle((gc[0],gc[1]), 1.2*Ratom, color = 'grey', alpha = 0.7 )
    blue   = plt.Circle((bc[0],bc[1]), 0.8*Ratom, color = 'blue', alpha = 0.7 )
    grey_v = plt.arrow(gc[0], gc[1], 2.0*Ratom*gc_v[0]/norm(gc_v), 2.0*Ratom*gc_v[1]/norm(gc_v), color = 'grey')
    blue_v = plt.arrow(bc[0], bc[1], 2.0*Ratom*bc_v[0]/norm(bc_v), 2.0*Ratom*bc_v[1]/norm(bc_v), color = 'blue')
    ax.set_aspect( 1 )
    ax.add_artist( blue )
    ax.add_artist( grey )
    ax.set_xlim([-0.15,0.15])
    ax.set_ylim([-0.15,0.15])
    ax.set_title( 'Case 1: Head on Collision - (time in microseconds)' )
    ax.annotate(f't = {curr_time*1e6:.2f}',xy=(-0.125,0.125))
    frame_name = 'headon_'+f'{frame_counter:02}'+'.png'
    #print(frame_name)
    fig.savefig(frame_name,dpi=500)
    ax.cla()
    if norm(bc-gc) < 2*Ratom:
        for q in range(10):
            grey       = plt.Circle((gc[0],gc[1]), 1.2*Ratom, color = 'grey', alpha = 0.7 )
            blue       = plt.Circle((bc[0],bc[1]), 0.8*Ratom, color = 'blue', alpha = 0.7 )
            grey_vi    = plt.arrow(gc[0], gc[1], 2.0*Ratom*gc_v[0]/norm(gc_v), 2.0*Ratom*gc_v[1]/norm(gc_v), color = 'grey',alpha = 1)
            blue_vi    = plt.arrow(bc[0], bc[1], 2.0*Ratom*bc_v[0]/norm(bc_v), 2.0*Ratom*bc_v[1]/norm(bc_v), color = 'blue',alpha = 1)
            grey_r     = plt.arrow(bc[0],bc[1],gc[0]-bc[0],gc[1]-bc[1],color='black',head_width=0)
            ax.set_aspect( 1 )
            ax.add_artist( blue )
            ax.add_artist( grey )
            ax.set_xlim([-0.15,0.15])
            ax.set_ylim([-0.15,0.15])
            ax.set_title( 'Case 1: Head on Collision - (time in microseconds)' )
            ax.annotate(f't = {curr_time*1e6:.2f}',xy=(-0.125,0.125))
            frame_counter = frame_counter + 1
            frame_name = 'headon_'+f'{frame_counter:02}'+'.png'
            #print(frame_name)
            fig.savefig(frame_name,dpi=500)
            ax.cla()
        dtp, dtm  = determine_dt(bc,gc,bc_v,gc_v)
        deltat    = max(dtp,dtm)
        bc        = bc - deltat*bc_v
        gc        = gc - deltat*gc_v
        curr_time = curr_time - deltat
        collision_flag = 1
        for q in range(10):
            grey       = plt.Circle((gc[0],gc[1]), 1.2*Ratom, color = 'grey', alpha = 0.7 )
            blue       = plt.Circle((bc[0],bc[1]), 0.8*Ratom, color = 'blue', alpha = 0.7 )
            grey_vi    = plt.arrow(gc[0], gc[1], 2.0*Ratom*gc_v[0]/norm(gc_v), 2.0*Ratom*gc_v[1]/norm(gc_v), color = 'grey',alpha = 1)
            blue_vi    = plt.arrow(bc[0], bc[1], 2.0*Ratom*bc_v[0]/norm(bc_v), 2.0*Ratom*bc_v[1]/norm(bc_v), color = 'blue',alpha = 1)
            grey_r     = plt.arrow(bc[0],bc[1],gc[0]-bc[0],gc[1]-bc[1],color='black',head_width=0)
            ax.set_aspect( 1 )
            ax.add_artist( blue )
            ax.add_artist( grey )
            ax.set_xlim([-0.15,0.15])
            ax.set_ylim([-0.15,0.15])
            ax.set_title( 'Case 1: Head on Collision - (time in microseconds)' )
            ax.annotate(f't = {curr_time*1e6:.2f}',xy=(-0.125,0.125))
            frame_counter = frame_counter + 1
            frame_name = 'headon_'+f'{frame_counter:02}'+'.png'
            #print(frame_name)
            fig.savefig(frame_name,dpi=500)
            ax.cla()
    if collision_flag == 1:
        #pdb.set_trace()
        collision_flag = 0
        print('Pre-collision Velocities:\n','  gc_v: ',gc_v,norm(gc_v),'\n','  bc_v: ',bc_v,norm(bc_v))
        #find Sherwood's answer
        new_gc_v, new_bc_v = calculate_post_bounce_conditions_Sherwood(gc,bc,gc_v,bc_v)
        print("Sherwood's computation: \n",' new_gc_v: ',new_gc_v,norm(new_gc_v),'\n','  new_bc_v: ',new_bc_v,norm(new_bc_v))
        #here is where we make the bounce happen
        velocity_scale = 15
        Mtot           = mass + mass
        Rcm            = (mass*gc + mass*bc)/Mtot
        rho_1          = gc - Rcm
        rho_2          = bc - Rcm
        rho_1_hat      = rho_1/norm(rho_1)
        rho_2_hat      = rho_2/norm(rho_2)
        ptot           = mass*gc_v + mass*bc_v
        pi_1           = mass*gc_v - (mass/Mtot)*ptot
        pi_2           = mass*bc_v - (mass/Mtot)*ptot
        rrel_star_hat  = (gc - bc)/norm(gc - bc)
        vrel_star      = (gc_v - bc_v)
        deltaP1        = -2*pi_1.dot(rho_1_hat)*rho_1_hat
        deltaP2        = -2*pi_2.dot(rho_2_hat)*rho_2_hat
        pi_1           = pi_1 + deltaP1
        pi_2           = pi_2 + deltaP2
        p1             = pi_1 + (mass/Mtot)*ptot
        p2             = pi_2 + (mass/Mtot)*ptot
        gc_new_v       = p1/mass
        bc_new_v       = p2/mass
        print("My computation: \n",'  gc_new_v: ',gc_new_v,norm(gc_new_v),'\n','  bc_new_v: ',bc_new_v,norm(bc_new_v))                
        for q in range(10):
            grey       = plt.Circle((gc[0],gc[1]), 1.2*Ratom, color = 'grey', alpha = 0.7 )
            blue       = plt.Circle((bc[0],bc[1]), 0.8*Ratom, color = 'blue', alpha = 0.7 )
            grey_vi    = plt.arrow(gc[0], gc[1], 2.0*Ratom*gc_v[0]/norm(gc_v), 2.0*Ratom*gc_v[1]/norm(gc_v), color = 'grey',alpha = 1 - q/10)
            blue_vi    = plt.arrow(bc[0], bc[1], 2.0*Ratom*bc_v[0]/norm(bc_v), 2.0*Ratom*bc_v[1]/norm(bc_v), color = 'blue',alpha = 1 - q/10)            
            grey_vf    = plt.arrow(gc[0], gc[1], 2.0*Ratom*gc_new_v[0]/norm(gc_new_v), 2.0*Ratom*gc_new_v[1]/norm(gc_new_v), color = 'grey',alpha = q/10)
            blue_vf    = plt.arrow(bc[0], bc[1], 2.0*Ratom*bc_new_v[0]/norm(bc_new_v), 2.0*Ratom*bc_new_v[1]/norm(bc_new_v), color = 'blue',alpha = q/10)
            grey_r     = plt.arrow(bc[0],bc[1],gc[0]-bc[0],gc[1]-bc[1],color='black',head_width=0)
            ax.set_aspect( 1 )
            ax.add_artist( blue )
            ax.add_artist( grey )
            ax.set_xlim([-0.15,0.15])
            ax.set_ylim([-0.15,0.15])
            ax.set_title( 'Case 1: Head on Collision - (time in microseconds)' )
            ax.annotate(f't = {curr_time*1e6:.2f}',xy=(-0.125,0.125))
            frame_counter = frame_counter + 1
            frame_name = 'headon_'+f'{frame_counter:02}'+'.png'
            #print(frame_name)
            fig.savefig(frame_name,dpi=500)
            ax.cla()
        gc_v = new_gc_v
        bc_v = new_bc_v
        bc        = bc + deltat*bc_v
        gc        = gc + deltat*gc_v
        curr_time = curr_time + deltat
        for q in range(10):
            grey       = plt.Circle((gc[0],gc[1]), 1.2*Ratom, color = 'grey', alpha = 0.7 )
            blue       = plt.Circle((bc[0],bc[1]), 0.8*Ratom, color = 'blue', alpha = 0.7 )
            grey_vf    = plt.arrow(gc[0], gc[1], 2.0*Ratom*gc_new_v[0]/norm(gc_new_v), 2.0*Ratom*gc_new_v[1]/norm(gc_new_v), color = 'grey',alpha = 1)
            blue_vf    = plt.arrow(bc[0], bc[1], 2.0*Ratom*bc_new_v[0]/norm(bc_new_v), 2.0*Ratom*bc_new_v[1]/norm(bc_new_v), color = 'blue',alpha = 1)
            #grey_r     = plt.arrow(bc[0],bc[1],gc[0]-bc[0],gc[1]-bc[1],color='black',head_width=0)
            ax.set_aspect( 1 )
            ax.add_artist( blue )
            ax.add_artist( grey )
            ax.set_xlim([-0.15,0.15])
            ax.set_ylim([-0.15,0.15])
            ax.set_title( 'Case 1: Head on Collision - (time in microseconds)' )
            ax.annotate(f't = {curr_time*1e6:.2f}',xy=(-0.125,0.125))
            frame_counter = frame_counter + 1
            frame_name = 'headon_'+f'{frame_counter:02}'+'.png'
            #print(frame_name)
            fig.savefig(frame_name,dpi=500)
            ax.cla()
    curr_time     = curr_time + dt
    frame_counter = frame_counter + 1

import os
os.system('cmd /c "del MB_overlap_and_bounce_on_contact.mp4"')
os.system('cmd /c "ffmpeg -r 5 -start_number 0 -i headon_%02d.png -c:v libx264 -r 80 -pix_fmt yuv420p MB_overlap_and_bounce_on_contact.mp4"');
os.system('cmd /c "del *.png"')