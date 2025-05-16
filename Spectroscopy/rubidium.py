import numpy as np
from scipy.constants import hbar, epsilon_0, c, m_e, e, eV
import matplotlib.pyplot as plt
from scipy.constants import hbar, epsilon_0, c, m_e, e, eV
import sys
from typing import Callable
sys.path.insert(0, '../../MiniPys/Formatter')
import minipy_formatter as MF
MF.Format().rcUpdate()

#%% Format functions
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


def get_all_allowed(L, S, I):
    J = get_J(L, S)
    print("J: ", J)

    F = []
    for j in J:
        F.extend(get_F(j, I))
    F = list(set(F))
    print("F: ", F)


    K = get_K(F, I, J)
    print("K: ", K)
    return J, F, K


#%%

class HyperfineStructure:
    # Fine structure energy values
    # Same labeling as with Rubidium class dicts.
    
    second_crossover = (229.8518 - 193.7408) / 2 # MHz
    third_crossover = 266.650 / 2 # MHz

class FSLevel:
    L_dictionary = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H'}
    def __init__(self, L, J, I, S=1/2, n=5, A=None, B=None, Rbnumber=None):
        self.n = n
        self.L = L
        self.J = J
        self.I = I
        self.S = S
        self.Rbnumber = Rbnumber

        if not isinstance(self, HFSLevel):
            HFS = []
            for F in np.arange(np.abs(J - I), J + I + 1, 1):
                energy_level = HFSLevel(L, J, F, I, A, B)
                HFS.append(energy_level)
            self.HFS = HFS

    def graph_HFS(self, fig, ax, E0=0, label=True, E_units='eV'):
        # TODO: Make this consistent with graph_fine_strucutre; options for length, sep.
        full_line = [0, 4]
        half_line = [1.5, 4]
        E0_line = [0, 1.5]
        E0HFS_line = [1.5, 2.5]
        HFS_line = [2.5, 4]

        ax.plot(E0_line, [E0, E0], color='k', lw=3)
        ax.plot(half_line, [E0, E0], color='k', ls='--')

        if label:
            name = self.get_term_symbol()
            plt.text(np.average(E0_line), E0,
                         f"${name}$, $E_0$ = {E0:.2E} {E_units}",
                         ha='center', va='bottom')
        
        num_F = len(self.HFS)
        current_F = None
        F_pos = 0
        for level in self.HFS:
            F_pos += 1
            # Avoids extra graphing due to m.
            if level.F == current_F:
                current_F = level.F
                continue

            # Structure
            E1 = level.HFS
            E1_scaled = E1
            HFS = E0 + E1_scaled
            
            ax.plot(E0HFS_line, [E0, HFS], color='#444', ls=':')
            ax.plot(HFS_line, [HFS, HFS], color='k', lw=3)
            ax.set(xlim=(-0.5, 5))
            
            if label:
                # Energies
                spacing = (HFS_line[1] - HFS_line[0]) / num_F
                energy_pos = HFS_line[0] + (F_pos - 1 / 2) * spacing
                ax.arrow(energy_pos, E0, 0, E1_scaled, color='#ccc', lw=2,
                         length_includes_head=True, head_width=0, head_length=0)
                plt.text(energy_pos, HFS, f"    F = {int(level.F)}, " \
                         f"{E1:.1f} {E_units}", ha='left',
                         va='bottom')
                    
            # Format
            hide_spines(ax)
            ax.set(title=f'Rubidium {self.Rbnumber} HFS for ${name}$')

            current_F = level.F

    def get_term_symbol(self):
        '''Of the form n L_symbol_{J} F'''
        n = self.n
        s_mult = 2 * self.S + 1
        L_symbol = FSLevel.L_dictionary[int(self.L)]
        J_str = f"{int(2 * self.J)}/2"
        return f"{n}{L_symbol}_{{{J_str}}}"

    def __repr__(self):
        return f"{self.get_term_symbol()} = {self.HFS}"


class HFSLevel(FSLevel):

    def get_hyperfine_shift(A, B, K, I, J):
        summand1 = A * K / 2
        if B == 0: return summand1
        numerator = 3 * K * (K + 1) - 2 * I * (I + 1) * J * (J + 1)
        denominator = 8 * I * (2 * I - 1) * J * (2 * J - 1)
        shift = summand1 + B * numerator / denominator 

        return shift

    def __init__(self, L, J, F, I, A, B, S=1/2, n=5):
        super().__init__(L, J, I, S=S, n=n)
        self.F = F
        K = F * (F + 1) - I * (I + 1) - J * (J + 1)
        self.HFS = HFSLevel.get_hyperfine_shift(A, B, K, I, J)

    def get_term_symbol(self, F_only=True):
        '''Of the form n L_symbol_{J} F'''
        F = int(self.F)
        if F_only: return f"{F=}"
        term_symbol = super().get_term_symbol()
        return f"{term_symbol} {F=}"
        
    def __repr__(self):
        return self.get_term_symbol()

class Transition:
    def __init__(self, psi1, psi2, E_diff: float, E_units='eV'):
        '''
        Here E_units is the unit used for E_diff.
        '''
        self.initial = psi1
        self.final = psi2
        self.E_diff = E_diff
        self.E_units = E_units

    def __repr__(self):
        return f"{self.initial}->{self.final} - {self.E_diff:.2f} {self.E_units}"        
        
class Rubidium:
    # TODO: Add a way to calculate gross structure and fine structure.
    # Constants
    h = 6.626 * 10**(-34)
    h = 1
    A_S12_87 = h * 3.417341305 * 1E3 # MHz
    A_P12_87 = h * 408.328 # MHz
    A_P32_87 = h * 84.7185 # MHz
    B_P32_87 = h * 12.4965 # MHz
    E_diff_S12_P32_87 = 384.230406373 * 1E9 # MHz

    A_S12_85 = h * 1.011910813 * 1E3 # MHz
    A_P12_85 = h * 120.527 # MHz
    A_P32_85 = h * 25.002 # MHz
    B_P32_85 = h * 25.790 # MHz
    E_diff_S12_P32_85 = 384.2304844685 * 1E9 # MHz

    L_dictionary = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H'}

    # Constants dictionary
    constants_85 = {
            "5S_1/2": A_S12_85,
            "5P_1/2": A_P12_85,
            "5P_3/2": A_P32_85,
            "B5S_1/2": 0,
            "B5P_1/2": 0,
            "B5P_3/2": B_P32_85
    }

    constants_87 = {
            "5S_1/2": A_S12_87,
            "5P_1/2": A_P12_87,
            "5P_3/2": A_P32_87,
            "B5S_1/2": 0,
            "B5P_1/2": 0,
            "B5P_3/2": B_P32_87
    }
    
    # Quantum numbers
    Ls = [0, 1]
    S = 1 / 2
    I_85 = 5 / 2
    I_87 = 3 / 2

    
    def __init__(self, Rbnumber: int):
        if Rbnumber not in [85, 87]:
            raise f'Specified isotope number ({Rbnumber}) not in database.'
        if Rbnumber == 85:
            self.Rbnumber = Rbnumber
            self.I = Rubidium.I_85
            self.constants = Rubidium.constants_85
        else:
            self.I = Rubidium.I_87
            self.constants = Rubidium.constants_87
            
        self.Rbnumber = Rbnumber
        self.S = Rubidium.S
        self.Ls = Rubidium.Ls
        self.HFS = self.get_HFS()
    
    def get_HFS(self):
        I = self.I
        S = self.S
        Rbnumber = self.Rbnumber
        HFS = []
        for L in self.Ls:
            for J in np.arange(np.abs(L - S), L + S + 1, 1):
                key = f"5{self.L_dictionary[int(L)]}_{int(2 * J)}/2"
                B_key = "B" + key
                A, B = self.constants[key], self.constants[B_key]
                HFS.append(FSLevel(L, J, I, S, 5, A, B, Rbnumber))
        return HFS

    def get_allowed_transitions(state1: FSLevel, state2: FSLevel, Rbnumber, E_units='MHz'):
        # As of now, this only calculates diff between ground state and 5P3/2
        transitions = []
        for level1 in state1.HFS:
            for level2 in state2.HFS:
                if not int(level1.F - level2.F) in [-1, 0, 1]:
                    continue
                if Rbnumber == 85: E_diff = Rubidium.E_diff_S12_P32_85
                else: E_diff = Rubidium.E_diff_S12_P32_87
                E_diff += np.abs(level2.HFS - level1.HFS)
                transitions.append(Transition(level1, level2, E_diff, E_units))
        return transitions


#%%
def main():
    Rb85 = Rubidium(85)
    Rb87 = Rubidium(87)

    fig, ax = plt.subplots()
    # rubidium_87.HFS[2].graph_HFS(fig, ax, E_units='MHz')
    trans85 = Rubidium.get_allowed_transitions(Rb85.HFS[0], Rb85.HFS[2], 85)
    trans87 = Rubidium.get_allowed_transitions(Rb87.HFS[0], Rb87.HFS[2], 87)
    trans85 = [trans.E_diff for trans in trans85]
    trans87 = [trans.E_diff for trans in trans87]
    
    print(min(trans87) - max(trans87))
    #ax.scatter(trans85, len(trans85) * [1], color='k')
    #ax.scatter(trans87, len(trans87) * [1])
    #plt.show()


    x = [0, 1.29 * 1E3, 3.03 * 1E3, 2.5 * 1E3]
    ax.scatter(x, len(x) * [1])
    plt.show()
    

if __name__ == '__main__':
    main()