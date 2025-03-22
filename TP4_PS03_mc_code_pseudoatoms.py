import numpy as np
import matplotlib.pyplot as plt

def get_initial_positions(N_atoms,box_size):
  return np.random.uniform(0,box_size,size=(N_atoms,2))

def get_dist(pos1,pos2):
  return ((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)**0.5

def get_energy(pos):
  en=0.
  for i in range(len(pos)):
   for j in range(i+1,len(pos)):
    dist=get_dist(pos[i],pos[j])**4
    en+=4./dist**2-4./dist
  return en

def try_step(pos,iatom,length,box_size):
 radius = np.random.uniform(0,length)
 angle = np.random.uniform(0,2*np.pi)
 pos[iatom][0] += radius*np.cos(angle)
 pos[iatom][1] += radius*np.sin(angle)
 if((pos[iatom][0]>box_size)or(pos[iatom][0]<0.)or(pos[iatom][1]>box_size)or(pos[iatom][1]<0.)): en=1e18
 else: en=get_energy(pos)
 return(pos,en)

# MAIN CODE STARTS HERE

N_atoms = 40
box_size = 10
n_steps = 30001
temperature = 10.
pos_old=get_initial_positions(N_atoms,box_size)
en_old=get_energy(pos_old)

print ("Starting energy:", en_old)

plt.ion()
plt.xlabel('x')
plt.ylabel('y')
plt.xlim([-1, 11])
plt.ylim([-1, 11])

for istep in range(n_steps):
 (pos_new,en_new)=try_step(np.copy(pos_old),np.random.choice(N_atoms),5.,box_size)
 u=np.random.uniform(0,1)
 if (u<np.exp(-(en_new-en_old)/1./temperature)):
   pos_old=np.copy(pos_new)
   en_old=en_new
 if(istep%1000==0):
   plt.scatter(pos_old[:N_atoms,0], pos_old[:N_atoms,1], alpha=.2, s=400)
   plt.xlim([-1, 11])
   plt.ylim([-1, 11])
   plt.draw()
   plt.pause(0.0001)
   plt.clf()
   print("Step",istep,"Energy:",en_old)
