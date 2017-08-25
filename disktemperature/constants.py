import astropy.constants as co
import astropy.units as u
import numpy as np

pi     = np.pi
mu     = 2.3

AU     = u.au.in_units('cm')

R_sun  = co.R_sun.cgs.value
M_sun  = co.M_sun.cgs.value
L_sun  = co.L_sun.cgs.value
Grav   = co.G.cgs.value
k_b    = co.k_B.cgs.value
m_p    = co.m_p.cgs.value
sig_sb = co.sigma_sb.cgs.value

T_sun  = (L_sun/(4*pi*R_sun**2*sig_sb))**0.25