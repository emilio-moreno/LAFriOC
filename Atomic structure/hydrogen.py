import sympy as smp
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, epsilon_0, c, m_e, e, eV
import sys
from typing import Callable
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

class Eigenstate:
    # TODO: Most of this only works for m_L=None, m_s=None (repr, etc)
    L_dictionary = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H'}
    def __init__(self, n, L, j, s=1/2, m_L=None, m_s=None, m_j=None,
                 term_symbol=True, s_notation=False):
        self.n = n
        self.L = L
        self.j = j
        self.s = s
        self.m_L = m_L
        self.m_s = m_s
        self.m_j = m_j
        self.term_symbol = term_symbol
        self.s_notation = s_notation

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
        if EV: return E / eV
        return E

    def ZeemanE1(self, Bext, EV=True):
        E = self.Lande_g() * self.m_j * mu_B * Bext
        if EV: return E / eV
        return E

    def get_term_symbol(self):
        '''Of the form n ^{s multiplicity} L_symbol_{j}'''
        n = self.n
        s_mult = 2 * self.s + 1
        L_symbol = Eigenstate.L_dictionary[int(self.L)]
        j_str = f"{int(2 * self.j)}/2"
        if not self.s_notation: return f"{n}{L_symbol}_{j_str}"
        return f"{n}^{s_mult}{L_symbol}_{j_str}"

    def get_ket_symbol(self):
        ket_symbol = f"|n={self.n} L={int(self.L)} j={int(2 * self.j)}/2"
        if self.s_notation: ket_symbol + f" s={int(2 * self.s)}/2"
        if self.m_j: ket_symbol + f" m={int(2 * self.m_j)}/2"
        return ket_symbol + "⟩"
        
    def __repr__(self):
        if self.term_symbol: return self.get_term_symbol()
        return self.get_ket_symbol()
        
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

# TODO: Option to graph transitions in wavelength or frequency graphs.
# TODO: Draw arrows in fine structure diagrams for allowed transitions.
class Transition:
    qnumbers = {'L', 'm_L', 'j', 'm_j'}

    def __init__(self, psi1: Eigenstate, psi2: Eigenstate, diff_E: float, E_units='eV'):
        '''
        Here E_units is the unit used for diff_E.
        '''
        self.initial = psi1
        self.final = psi2
        self.diff_E = diff_E
        self.E_units = E_units

    def satisfies_rule(q1: int, q2: int, qnumber: str):
        '''q1, q2 are any pair of quantum numbers.'''  
        if not qnumber in Transition.qnumbers:
            raise ValueError(f"Quantum number '{qnumber}' is not defined.")
        diff_k = int(q1 - q2)
        if qnumber == 'L': return diff_k in [-1, 1]
        if qnumber in ['m_L', 'j', 'm_j']: return diff_k in [-1, 0, 1]

    def is_allowed(psi1: Eigenstate, psi2: Eigenstate, qnumbers: set):
        if not set(qnumbers).issubset(Transition.qnumbers):
            raise ValueError(f"Undefined quantum numbers contained in {qnumbers}.")
        for qnumber in qnumbers:
            q1 = getattr(psi1, qnumber)
            q2 = getattr(psi2, qnumber)

            if q1 == None and q2 == None:
                continue
            if not Transition.satisfies_rule(q1, q2, qnumber):
                return False
        return True

    def __repr__(self):
        return f"{self.initial}->{self.final} - {self.diff_E:.5f} {self.E_units}"


class Energylevel:
    # TODO: I'll have to implement this so that it can accomodate for different

    def get_hydrogen_eigenstates(n, m_j, term_symbol):
        eigenstates = []
        for L in np.arange(0, n, 1, dtype=float):
            for j in np.arange(np.abs(L - 1 / 2), L + 1 / 2 + 1, 1,
                               dtype=float):
                if m_j:
                    for m_j in np.arange(-j, j + 1, 1):
                        eigenstates.append(Eigenstate(n, L, j, 1 / 2, m_j=m_j,
                            term_symbol=term_symbol))
                else:
                    eigenstates.append(Eigenstate(n, L, j, 1 / 2, term_symbol=term_symbol))
        return eigenstates

    def __init__(self, n, m_j=True, term_symbol=True, get_eigenstates):
        # n, L, j
        self.n = n
        self.eigenstates = get_eigenstates(n, m_j, term_symbol)
    
    @staticmethod
    def convert_energy(E, unit):
        if unit not in ['eV', 'MHz', 'GHz', 'THz']:
            raise ValueError(f"Unit: '{unit}' not allowed.")
        if unit == 'eV': return E
        if 'Hz' in unit:
            E = E / (2 * PI * (hbar / eV)) # Hz
            if unit == 'MHz': return E / 1E6 # MHz
            if unit == 'GHz': return E / 1E9 # GHz
            if unit == 'THz': return E / 1E12 # THz

    get_diff_E_type = Callable[[Eigenstate, Eigenstate], float]
    def get_allowed_transitions(level1, level2, get_diff_E: get_diff_E_type,
                                qnumbers: set, E_units='eV'):
        '''Determines allowed transitions from two energy levels (n1, n2).'''
        allowed = []
        for psi1 in level1.eigenstates:
            for psi2 in level2.eigenstates:
                if Transition.is_allowed(psi1, psi2, qnumbers):
                    diff_E = get_diff_E(psi1, psi2)
                    allowed.append(Transition(psi1, psi2, diff_E, E_units))
        allowed.sort(key=lambda transition: transition.diff_E)                    
        return allowed
    
    def graph_FS_energies(self, fig, ax, scale=None, label=True, E_units='eV'):
        # TODO: Make this consistent with graph_fine_strucutre; options for length, sep.
        if not scale:
            scale = lambda j: 1
        full_line = [0, 4]
        half_line = [1.5, 4]
        E0_line = [0, 1.5]
        E0FS_line = [1.5, 2.5]
        EFS_line = [2.5, 4]

        E0 = self.convert_energy(self.eigenstates[0].HydrogenE0(), E_units)
        num_j = self.n
        ax.plot(E0_line, [E0, E0], color='k', lw=3)
        ax.plot(half_line, [E0, E0], color='k', ls='--')

        if label:
            plt.text(np.average(E0_line), E0,
                         f"$n$ = {self.n}, $E_0$ = {E0:.2E} {E_units}",
                         ha='center', va='bottom')
        
        current_j = None
        for psi in self.eigenstates:
            # Avoids extra graphing due to L.
            if psi.j == current_j:
                current_j = psi.j
                continue

            # Structure
            E1 = self.convert_energy(psi.HydrogenEFS1(), E_units)
            E1_scaled = E1 * scale(psi.j)
            EFS = E0 + E1_scaled
            
            ax.plot(E0FS_line, [E0, EFS], color='#444', ls=':')
            ax.plot(EFS_line, [EFS, EFS], color='k', lw=3)
            ax.set(xlim=(-0.5, 5))
            
            if label:
                # Energies
                spacing = (EFS_line[1] - EFS_line[0]) / num_j
                energy_pos = EFS_line[0] + ((2 * psi.j - 1) / 2 + 1 / 2) * spacing
                ax.arrow(energy_pos, E0, 0, E1_scaled, color='#ccc', lw=2,
                         length_includes_head=True, head_width=0, head_length=0)
                plt.text(energy_pos, EFS, f"    j = {int(2 * psi.j)} / 2, " \
                         f"$E_{{FS}}^1$ = {E1:.2E} {E_units}\n", ha='left',
                         va='bottom')
                    
            # Format
            hide_spines(ax)
            ax.set(title=f'Hydrogen Energy Fine Structure for $n$ = {psi.n}')

            current_j = psi.j

    def graph_fine_structure(self, fig, ax, label=True, E_units='eV',
                             L_length: float = 1, L_sep: float = 0.5, max_padding=0):
        E0 = self.convert_energy(self.eigenstates[0].HydrogenE0(), E_units)
        L_length = self.convert_energy(L_length, E_units)
        L_sep = self.convert_energy(L_sep, E_units)
        E0_line = [0, (self.n - 1) * L_sep + self.n * L_length + max_padding]

        ax.plot(E0_line, [E0, E0], color='k', lw=1, ls=':')
        if label:
            plt.text(-0.1 * L_length, E0,
                         f"$E_0$ = {E0:.2E} {E_units}, $n$ = {self.n}", 
                         ha='right', va='center')

        current_j = None
        current_L = None
        for psi in self.eigenstates:
            # Avoids extra graphing due to m_j.
            if psi.j == current_j and psi.L == current_L:
                current_j = psi.j
                current_j = psi.L
                continue

            E1 = self.convert_energy(psi.HydrogenEFS1(), E_units)
            EFS = E0 + E1
            
            L_line = [psi.L * (L_sep + L_length),
                      psi.L * L_sep + (psi.L + 1) * L_length]
                     
            # Graphing
            ax.plot(E0_line, [EFS, EFS], color='#ccc', lw=1, ls='--', zorder=1)
            ax.plot(L_line, [EFS, EFS], color='k', lw=2, zorder=2)

            if label:
                # Energies
                plt.text(np.average(L_line), E0, f"\n$l$ = {int(psi.L)}",
                         va='top', ha='center')
                plt.text(-0.1 * L_length, EFS,
                         f"$E_{{FS}}^1$ = {EFS:.5E} {E_units}, $j$ = {int(2 * psi.j)} / 2", ha='right',
                         va='center')
            
            # Format
            hide_spines(ax)
            ax.set(title=f'Fine Structure for $n$ = {psi.n}',
                   xlim=(-L_length, E0_line[1] + max_padding))

            current_j = psi.j
            current_j = psi.L
            
        
    def graph_zeeman_effect(self, fig, ax, label=True, E_units='eV', Bmax=1):
        # TODO: Extract from here code for FS only.
        E0 = self.convert_energy(self.eigenstates[0].HydrogenE0(), E_units)
        Bmax = Bmax * self.eigenstates[-1].HydrogenEFS1() / 4
        Bmax = self.convert_energy(Bmax, E_units)

        L_length = 2 * Bmax
        L_sep = 4 * Bmax
        for psi in self.eigenstates:
            # Bmax is the maximum value for the external (normalized) electric field.
            # A small number needs to be considered so fine-structure can be discerned
            # alongside Zeeman splitting.
            self.graph_fine_structure(fig, ax, label, E_units, L_length, L_sep,
                max_padding=Bmax)
            
            E1 = self.convert_energy(psi.HydrogenEFS1(), E_units)
            EFS = E0 + E1
            
            L_line_max = psi.L * L_sep + (psi.L + 1) * L_length
            slope = psi.ZeemanE1(1) / (mu_B / eV) # Normalized
            
            # Graphing Zeeman splittings
            zeeman_line = [L_line_max, L_line_max + Bmax]
            ax.plot(zeeman_line, [EFS, EFS + slope * Bmax], color='k', lw=2)

            if label:
                plt.text(L_line_max + Bmax, EFS + slope * Bmax,
                         f"    $m_j$ = {int(2 * psi.m_j)}/2, $m_j g_{{J}}$ = {slope:.2f}",
                         ha='left', va='center')
            
            # Format
            hide_spines(ax)
            ax.set(title=f'Weak Zeeman Splitting for $n$ = {psi.n}',
                xlim=(-Bmax, L_line_max + 2 * Bmax))
            
            
                            
    def __repr__(self):
        return f"{[str(psi) for psi in self.eigenstates]}"
