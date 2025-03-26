import numpy as np

# %% Constants
h = 6.626 * 10**(-34)
A_S12 = h * 3.417 * 10**3 # MHz
A_P12 = h * 408.328 # MHz
A_P32 = h * 84.7185 # MHz
B_P32 = h * 12.4965 # MHz


def get_J(L, S):
    return np.arange(np.abs(L - S), L + S + 1, 1)


@np.vectorize
def get_F(J, I):
    return np.arange(np.abs(J - I), J + I + 1, 1)


def get_K(F, I, J):
    return F * (F + 1) - I * (I + 1) - J * (J + 1)

def hyperfine_shifts(Ahfs, Bhfs, K, I, J):
    summand1 = Ahfs * K / 2
    numerator = 3 * K * (K + 1) - 2 * I * (I + 1) * J * (J + 1)
    denominator = 2 * I * (2 * I - 1) * 2 * J * (2 * J - 1)
    
    return summand1 + Bhfs * numerator / denominator


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

# %%
L = 1
S = 1 / 2
I = 1 / 2
get_all_allowed(L, S, I)
