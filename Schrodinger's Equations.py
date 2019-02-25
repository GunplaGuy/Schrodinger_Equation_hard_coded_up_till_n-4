import scipy.constants as c
import numpy as np

# Angular Wave functions

def angular_wave_func(m,l,theta,phi):
    i = np.complex(0,1)
    if l == 0:
        Y = np.sqrt(1/(4*np.pi))
    elif l == 1:
        if m == 0:
            Y = np.sqrt(3/(4*np.pi))*np.cos(theta)
        elif abs(m) == 1:
            Y = -1*(m/abs(m))*np.sqrt(3/(8*np.pi))*np.sin(theta)*np.exp(i*phi*m)
    elif l == 2:
        if m == 0:
            Y = np.sqrt(5/(16*np.pi))*((3*np.cos(theta)**2)-1)
        elif abs(m) == 1:
            Y = -1*(m/abs(m))*np.sqrt(15/(8*np.pi))*np.cos(theta)*np.sin(theta)*np.exp(m*i*phi)
        elif abs(m) == 2:
            Y = np.sqrt(15/(32*np.pi))*(np.sin(theta)**2)*np.exp(m*i*phi)
    elif l == 3:
        if m == 0:
            Y = np.sqrt(7/(16*np.pi))(5*(np.cos(theta)**3)-3*np.cos(theta))
        elif abs(m) == 1:
            Y = -1*(m/abs(m))*np.sqrt(21/(64*np.pi))*np.sin(theta)*(5*(np.cos(theta)**2)-1)*np.exp(m*i*phi)
        elif abs(m) == 2:
            Y = np.sqrt(105/(32*np.pi))*np.cos(theta)*(np.sin(theta)**2)*np.exp(m*i*phi)
        elif abs(m) == 3:
            Y = -1*(m/abs(m))*np.sqrt(35/(64*np.pi))*(np.sin(theta)**3)*np.exp(i*m*phi)
    Y = round(Y, 5)
    Y = np.complex(Y)
    return Y

# Radial Wave functions

a = c.physical_constants['Bohr radius'][0]

def radial_wave_func(n,l,r):
    if n == 1:
        if l == 0:
            R = 2*np.exp((-1*r)/a)
    elif n == 2:
        if l == 0:
            R = (1/np.sqrt(2))*(1-(r/(2*a)))*np.exp((-1*r)/(2*a))
        elif l == 1:
            R = (1/np.sqrt(24))*(r/a)*np.exp((-1*r)/(2*a))
    elif n == 3:
        if l == 0:
            R = (2/(81*np.sqrt(3)))*(27-18*(r/a)+2*(r/a)**2)*np.exp((-1*r)/(3*a))
        elif l == 1:
            R = (8/(27*np.sqrt(6)))*(1-(r/(6*a)))*(r/a)*np.exp(r/(-3*a))
        elif l == 2:
            R = (4/(81*np.sqrt(30)))*((r/a)**2)*np.exp(r/(-3*a))
    elif n == 4:
        if l == 0:
            R = 0.25*(1-0.75*(r/a)+(1/8)*((r/a)**2)-(1/192)*((r/a)**3))*np.exp(r/(-4*a))
        elif l == 1:
            R = (np.sqrt(5)/(16*np.sqrt(3)))*(r/a)*(1-0.25*(r/a)+(1/80)*(r/a)**2)*np.exp(r/(-4*a))
        elif l == 2:
            R = (1/(64*np.sqrt(5)))*((r/a)**2)*(1-(1/12)*(r/a))*np.exp(r/(-4*a))
        elif l ==3:
            R = (1/(768*np.sqrt(35)))*((r/a)**3)*np.exp(r/(-4*a))
    R = round(R, 5)
    return R
