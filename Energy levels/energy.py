import sympy as smp
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, epsilon_0, c, m_e, e, eV
import sys
sys.path.insert(0, '../../MiniPys/Formatter')
import minipy_formatter as MF
MF.Format().rcUpdate()

PI = np.pi
mu_B = e * hbar / 2 * m_e

class Hydrogen:
    @staticmethod
    def E0(n, EV=True):
        E = -m_e / (2 * hbar**2) * (e**2 / (4 * PI * epsilon_0))**2 * 1 / n**2
        if EV:
            return E / eV
        return E

    @staticmethod
    def E1(n, j, EV=True):
        E = Hydrogen.E0(n, EV = False)**2 / (2 * m_e * c**2) * (3 - 4 * n / (j + 1 / 2))
        if EV:
            return E / eV
        return E
    
def hide_spines(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.axes.get_xaxis().set_visible(False)
    ax.set_xticklabels([])
    ax.set_xticks([])

    ax.axes.get_yaxis().set_visible(False)
    ax.set_yticklabels([])
    ax.set_yticks([])

class Energylevel:
    
    def __init__(self, n, m_j=True):
        # n, L, j
        eigenstates = []
        for L in np.arange(0, n, 1, dtype=float):
            for j in np.arange(np.abs(L - 1 / 2), L + 1 / 2 + 1, 1,
                               dtype=float):
                for m_j in np.arange(-j, j + 1, 1):
                    eigenstates.append(Eigenstate(n, L, j, 1 / 2, m_j))
        self.n = n
        self.eigenstates = eigenstates
    
    @staticmethod
    def convert_energy(E, unit):
        if unit not in ['eV', 'MHz']:
            raise ValueError(f"Unit: '{unit}' not allowed.")
        if unit == 'eV': return E
        if unit == 'MHz':
            E = E / (2 * PI * (hbar / eV)) # Hz
            return E / 1E6 # MhZ
    
    def graph_structure(self, fig, ax, scale=None, label=True, E_units='eV'):
        if not scale:
            scale: lambda j: 1
        full_line = [0, 4]
        half_line = [1.5, 4]
        E0_line = [0, 1.5]
        E0FS_line = [1.5, 2.5]
        EFS_line = [2.5, 4]
        
        current_j = None
        
        for psi in self.eigenstates:
            if psi.j == current_j:
                current_j = psi.j
                continue
            # Structure
            E0 = self.convert_energy(psi.HydrogenE0(), E_units)
            ax.plot(E0_line, [E0, E0], color='k', linewidth=3)
            E1 = self.convert_energy(psi.HydrogenEFS1(), E_units)
            E1_scaled = E1 * scale(psi.j)
            EFS = E0 + E1_scaled
            ax.plot(half_line, [E0, E0], color='k', linestyle='--')
            ax.plot(E0FS_line, [E0, EFS], color='#444', linestyle=':')
            ax.plot(EFS_line, [EFS, EFS], color='k', linewidth=3)
            ax.set(xlim=(-0.5, 5))
            
            if label:
                # Energies
                num_j = psi.n
                spacing = (EFS_line[1] - EFS_line[0]) / num_j
                energy_pos = EFS_line[0] + ((2 * psi.j - 1) / 2 + 1 / 2) * spacing
                plt.text(np.average(E0_line), E0,
                         f"$n$ = {psi.n}, $E_0$ = {E0:.2E} {E_units}", ha='center',
                         va='bottom')
                ax.arrow(energy_pos, E0, 0, E1_scaled, color='#ccc', linewidth=2,
                         length_includes_head=True, head_width=0, head_length=0)
                plt.text(energy_pos, EFS, f"    j = {int(2 * psi.j)} / 2, " \
                         f"$E_{{FS}}^1$ = {E1:.2E} {E_units}\n", ha='left',
                         va='bottom')
                    
            # Format
            hide_spines(ax)
            ax.set(title=f'Hydrogen Fine Structure for $n$ = {psi.n}')
        
    def zeeman_effect(self, fig, ax, scale=None, label=True, E_units='eV'):
        
        current_L = None
        current_j = None
        for psi in self.eigenstates:
            E0 = psi.HydrogenE0()
            E1 = psi.HydrogenEFS1()
            EFS = E0 + E1
            
            z_length = 0.7E-6
            slope = psi.ZeemanE1(1) / (mu_B / eV) # Normalized
            L_length = 2 * z_length
            L_sep = 4 * z_length
            
            number_L = psi.n
            L_line = [psi.L * (L_sep + L_length),
                      psi.L * L_sep + (psi.L + 1) * L_length]
            E0_line = [0, (psi.n - 1) * L_sep + psi.n * L_length
                       + z_length]
            
            if label:
                # Energies
                plt.text(-1.5 * z_length, E0,
                         f"$E_0$ = {E0:.2E} {E_units}, $n$ = {psi.n}", ha='center',
                         va='center')
                plt.text(-1.5 * z_length, EFS,
                         f"$E_{{FS}}^1$ = {EFS:2f}, $j$ = {int(2 * psi.j)} / 2", ha='center',
                         va='center')
                plt.text(L_line[1] + z_length, EFS + slope * z_length,
                         f"    $m_j$ = {int(2 * psi.m_j)}/2, $m_j g_{{J}}$ = {slope:.2f}", ha='left',
                         va='center')
                plt.text(np.average(L_line), E0, f"\n$l$ = {int(psi.L)}",
                         va='top', ha='center')
            
            E0 = self.convert_energy(psi.HydrogenE0(), E_units)
            ax.plot(E0_line, [E0, E0], color='k', linewidth=1, ls=':')
            ax.plot(E0_line, [EFS, EFS], color='#ccc',
                    linewidth=1, ls='--', zorder=1)
            ax.plot(L_line, [EFS, EFS], color='k', linewidth=2, zorder=2)
            
            zeeman_line = [L_line[1], L_line[1] + z_length]
            ax.plot(zeeman_line, [EFS, EFS + slope * z_length],
                    color='k', linewidth=2)
            
            hide_spines(ax)
            ax.set(title=f'Weak Zeeman Splitting for $n$ = {psi.n}',
                   xlim=(-z_length, E0_line[1] + z_length))
            
            
                            
    def __repr__(self):
        return f"{[str(es) for es in self.eigenstates]}"
        

class Eigenstate:
    
    def __init__(self, n, L, j, s=1/2, m_j=None):
        self.n = n
        self.L = L
        self.j = j
        self.s = 1 / 2
        self.m_j = m_j

    def Lande_g(self):
    	j, L, s = self.j, self.L, self.s
    	g = j* (j + 1) - L * (L + +1) + s * (s + 1)
    	return 1 + g / (2 * j * (j + 1))
    
    def HydrogenE0(self, EV=True):
        n = self.n
        E = -m_e / (2 * hbar**2) * (e**2 / (4 * PI * epsilon_0))**2 * 1 / n**2
        if EV:
            return E / eV
        return E

    def HydrogenEFS1(self, EV=True):
        n, j = self.n, self.j
        E = Hydrogen.E0(n, EV = False)**2 / (2 * m_e * c**2) * (3 - 4 * n / (j + 1 / 2))
        if EV:
            return E / eV
        return E

    def ZeemanE1(self, Bext, EV=True):
        E = self.Lande_g() * self.m_j * mu_B * Bext
        if EV:
            return E / eV
        return E
        
    def __repr__(self):
        return f"|{self.n} {int(self.L)} {int(2 * self.j)}/2 " \
               f"{int(2 * self.s)}/2 {int(2 * self.m_j)}/2⟩"
    
#%%

fig, ax = plt.subplots(figsize=(24, 10), dpi=200)
scale = lambda j: 1
level = Energylevel(3)
#level.graph_structure(fig, ax, scale, E_units='eV')
level.zeeman_effect(fig, ax, scale, E_units='eV')
plt.savefig('Weak Zeeman Effect, n=3.png')
