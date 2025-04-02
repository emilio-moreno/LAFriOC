import numpy as np

#%% Quantum number functions

def get_J(L, S):
    return np.arange(np.abs(L - S), L + S + 1, 1)


@np.vectorize
def get_F(J, I):
    return np.arange(np.abs(J - I), J + I + 1, 1)


def get_K(F, I, J):
    return F * (F + 1) - I * (I + 1) - J * (J + 1)

def hyperfine_shifts(A, B, K, I, J):
    summand1 = A * K / 2
    if B == 0: return summand1
    numerator = 3 * K * (K + 1) - 2 * I * (I + 1) * J * (J + 1)
    denominator = 8 * I * (2 * I - 1) * J * (2 * J - 1)
    shift = summand1 + B * numerator / denominator 

    return shift


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
    
    constants_85 = {
            "01": A_S12_85,
            "11": A_P12_85,
            "13": A_P32_85
    }
    pass
    

class Rubidium:
    
    # Constants
    h = 6.626 * 10**(-34)
    h = 1
    A_S12_87 = h * 3.417341305 * 1E3 # MHz
    A_P12_87 = h * 408.328 # MHz
    A_P32_87 = h * 84.7185 # MHz
    B_P32_87 = h * 12.4965 # MHz

    A_S12_85 = h * 1.011910813 * 1E3 # MHz
    A_P12_85 = h * 120.527 # MHz
    A_P32_85 = h * 25.002 # MHz
    B_P32_85 = h * 25.790 # MHz

    # Constants dictionary

    # To access constant I'll use as a key "L 2*J".
    # E.g., L = 1, J = 3/2, then key = "13".
    # I use key = "Bxx" for the quadrupole constant which is
    # non-zero only for P3/2.

    constants_85 = {
            "01": A_S12_85,
            "11": A_P12_85,
            "13": A_P32_85,
            "B01": 0,
            "B11": 0,
            "B13": B_P32_85
    }

    constants_87 = {
            "01": A_S12_87,
            "11": A_P12_87,
            "13": A_P32_87,
            "B01": 0,
            "B11": 0,
            "B13": B_P32_87
    }
    
    # Quantum numbers
    Ls = [0, 1]
    S = 1 / 2
    I_85 = 5 / 2
    I_87 = 3 / 2

    
    def __init__(self, number: int):
        if number not in [85, 87]:
            raise f'Specified isotope number ({number}) not in database.'
        if number == 85:
            self.I = Rubidium.I_85
            self.constants = Rubidium.constants_85
        else:
            self.I = Rubidium.I_87
            self.constants = Rubidium.constants_87
            
        self.number = number
        self.S = Rubidium.S
        self.Ls = Rubidium.Ls
        self.hyperfine_structure = self.get_hyperfine_structure()
    
    def get_hyperfine_structure(self):
        hyperfine_dict = {}
        for L in self.Ls:
            Js = get_J(L, self.S)
            for J in Js:
                key = str(L) + str(int(2 * J))
                B_key = "B" + key
                A, B = self.constants[key], self.constants[B_key]
                Fs = get_F(J, self.I)
                hyperfine_dict[key] = []
                for F in Fs:
                    K = get_K(F, J, self.I)
                    shift = hyperfine_shifts(A, B, K, self.I, J)
                    hyperfine_dict[key].append(shift)
        return hyperfine_dict


#%%
rubidium_85 = Rubidium(85)
rubidium_87 = Rubidium(87)
print(rubidium_87.hyperfine_structure)