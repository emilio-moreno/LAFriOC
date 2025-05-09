'''Calculate numerical value of various quantities.'''
from scipy.constants import hbar, epsilon_0, c, m_e, e, eV
import numpy as np

Pi = np.pi

def ground_E(Z=2):
	return -Z**2 * e**4 * m_e / (2 * (4 * epsilon_0 * Pi * hbar)**2)

def ground_E_rep(Z=2):
	return Z * e**4 * m_e / (4 * Pi * epsilon_0 * hbar)**2


print(ground_E(14) / eV)
print(ground_E_rep(14) / eV)
print((ground_E(14) + ground_E_rep(14)) / eV)
print(2 * (2400 - 2285.8) / (2400 + 2285.8))
