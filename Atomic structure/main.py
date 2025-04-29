from atomic_structure import *

n1 = Energylevel(1, m_j=True)
n2 = Energylevel(2, m_j=True)
n3 = Energylevel(3, m_j=True)
n4 = Energylevel(4, m_j=True)
n5 = Energylevel(5, m_j=True)

def graph():
	for n_level in [n4, n5]:
		fig, ax = plt.subplots(figsize=(24, 10), dpi=200)
		n_level.graph_zeeman_effect(fig, ax, E_units='eV')
		plt.savefig(f'Weak Zeeman Effect/Weak Zeeman Effect, n={n_level.n}.png')

def tranistions():
	qnumbers = {"L", "j"}
	E_units = 'nm'

	def get_diff_E(psi1: Eigenstate, psi2: Eigenstate):
	    EFS1 = psi1.HydrogenE0() + psi1.HydrogenEFS1()
	    EFS2 = psi2.HydrogenE0() + psi2.HydrogenEFS1()
	    diff_E = EFS2 - EFS1 # eV
	    diff_E = diff_E / (2 * PI * hbar / eV) # Hz
	    diff_E = c / diff_E * 1E9 # nm
	    return diff_E

	allowed_45 = Energylevel.get_allowed_transitions(n4, n5, get_diff_E,
	                                                 qnumbers, E_units=E_units)
	print(allowed_45, len(allowed_45))
	
	# Energylevel.allowed_transitions(level1, level2)

graph()