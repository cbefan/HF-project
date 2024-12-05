import numpy as np
import matplotlib.pyplot as plt
import sys

plt.rcParams['figure.dpi'] = 150
plt.rcParams['figure.figsize'] = (6.4*1.2,4.8*1.2)
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['figure.constrained_layout.use'] = True

E_H = -0.466582     #table 3.6
E_He = -2.643876
E_He_plus = -1.975514

# plotting time 😊
f_name = sys.argv[1]
if f_name[8]=='e':   #both atoms are helium
    diff = 2*E_He_plus
    plt.ylabel('E(He2++) - 2E(He+) (a.u.)')
elif f_name[6]=='e':  #atom1 is helium 
    diff = E_He
    plt.ylabel('E(HeH+) - E(He) (a.u.)')
elif f_name[7]=='e': #atom2 is helium
    diff = E_He
    plt.ylabel('E(HeH+) - E(He) (a.u.)')
else:   # both atoms are hydrogen
    diff = 2*E_H   
    plt.ylabel('E(H2) - 2E(H) (a.u.)')

data = np.loadtxt(f_name,delimiter=',')
r = data[:,0]
e = data[:,1]
min_e = np.min(e)
min_r = r[e==min_e]
print(f'minimum energy {min_e} at r = {min_r[0]}')
e -=diff
plt.plot([-1,10],[0,0],'k--')
plt.plot(r,e)

plt.xlabel('R (a.u.)')
plt.ylim([-0.5,1])
plt.xlim([0.5,5])

plt.savefig(f'{f_name[5:-4]}.png')
plt.show()
# Figures 3.5 and 3.8