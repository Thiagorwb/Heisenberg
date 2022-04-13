from .gates import *


# IBM code, Returns the matrix representation of the XXX Heisenberg model for 3 spin-1/2 particles in a line
def H_heis3():
    # Interactions (I is the identity matrix; X, Y, and Z are Pauli matricies; ^ is a tensor product)
    XXs = (I^X^X) + (X^X^I)
    YYs = (I^Y^Y) + (Y^Y^I)
    ZZs = (I^Z^Z) + (Z^Z^I)
    
    # Sum interactions
    H = XXs + YYs + ZZs
    
    # Return Hamiltonian
    return H

#IBM code

# Returns the matrix representation of U_heis3(t) for a given time t assuming an XXX Heisenberg Hamiltonian for 3 spins-1/2 particles in a line
def U_heis3(t):
    # Compute XXX Hamiltonian for 3 spins in a line
    H = H_heis3()
    
    # Return the exponential of -i multipled by time t multipled by the 3 spin XXX Heisenberg Hamilonian 
    return (t * H).exp_i()

def H_12():

    # Sum interactions
    H = (X^X) + (Y^Y) + (Z^Z)
    # Return Hamiltonian
    return H

def U_12(t):
    # Compute XXX Hamiltonian for 3 spins in a line
    H = H_12()
    
    # Return the exponential of -i multipled by time t multipled by the 3 spin XXX Heisenberg Hamilonian 
    return (t * H).exp_i()

def Heis(n):
    
    dim = 2**n
    H = np.zeros((dim, dim), dtype=np.complex_)
    H12 = H_12().to_matrix()
        
    for i in range(n - 1):
        back = np.identity(2**(n - 2 - i ))
        front = np.identity(2**( i))
        term = np.kron(front, np.kron(H12, back))
        H += term
        
    return H


# classical simulation
def prepare_state(t):
    if t == np.pi:
        return (One^One^Zero)
    # Define initial state |110>
    initial_state = One^One^Zero
    return (U_heis3(float(t)) @ initial_state).eval()

def swap_H(n):
    swap = Operator(swap_gate()).data
    H = np.zeros((2**n, 2**n), dtype=np.complex_)
    for i in range(n-1):
        front = np.eye(2**(i))
        back = np.eye(2**(n - 2 - i))
        interaction = 2*swap - np.identity(4)
        term =  np.kron(front, np.kron(interaction, back)) 
        H += term
    return H