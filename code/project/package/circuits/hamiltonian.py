from .gates import *
from ..helpers import *


###################################################################################
## Helper Functions to generate Hamiltonians and their associated time evolution ##
###################################################################################

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


# Our code starts here

def H_12():

    # Sum interactions
    H = (X^X) + (Y^Y) + (Z^Z)
    # Return Hamiltonian
    return H

def U_12(t):
    # Compute XXX Hamiltonian for 2 spins in a line
    H = H_12()
    
    # Return the exponential of -i multipled by time t multipled by the 2 spin XXX Heisenberg Hamilonian 
    return (t * H).exp_i()

def Heis(n):
    """
    Generates the XXX Hamiltonian for n spins in a line
    Input: n,  int,  number of qubits
    Output: H,  2**n by 2**n numpy array, the Hamiltonian
    """
    
    dim = 2**n
    H = np.zeros((dim, dim), dtype=np.complex_)
    H12 = H_12().to_matrix()
        
    for i in range(n - 1):
        back = np.identity(2**(n - 2 - i ))
        front = np.identity(2**( i))
        term = np.kron(front, np.kron(H12, back))
        H += term
        
    return H

def U_Heis(t, n):
    """
    Generates the time evolution unitary for the n qubit XXX Heisenberg Hamiltonian

    Inputs: 
    - n,  int,  number of qubits
    - t, float, target time

    Output:
    2^n by 2^n numpy array, the time evolution unitary 

    """

    # Compute XXX Hamiltonian for n spins in a line
    H = Heis(n)

    # Return the exponential of -i multipled by time t multipled by the n spin XXX Heisenberg Hamilonian 
    return scipy.linalg.expm(- 1j * t * H)

def swap_H(n):
    """
    Another interesting description of the Heisenberg model
    """
    
    swap = Operator(swap_gate()).data
    H = np.zeros((2**n, 2**n), dtype=np.complex_)
    for i in range(n-1):
        front = np.eye(2**(i))
        back = np.eye(2**(n - 2 - i))
        interaction = 2*swap - np.identity(4)
        term =  np.kron(front, np.kron(interaction, back)) 
        H += term
    return H

#################################################################
## Helper Functions to Classically Simulate the Time Evolution ##
#################################################################

# the input state
def psi_0(num_qubits = 3):
    input_state = np.zeros(2**(num_qubits-2))
    input_state[0] = 1.0
    input_state = np.kron((One^One).to_matrix(), input_state)
    return input_state

# classical simulation of 3 qubit system
def prepare_state(t):

    """
    Function generates the expected quantum state after simulating the 3 qubit XXX Heisenberg model for time t, on input state |011>

    Input: t, float, the time
    Output: numpy array of floats of size 2**3 = 8, the output vector

    """

    # to match IBM code for submission, keep this if statement
    if t == np.pi:
        return (One^One^Zero).to_matrix()

    # Define initial state |110>
    initial_state = One^One^Zero
    return ((U_heis3(float(t)) @ initial_state).eval()).to_matrix()

def prepare_system(t = np.pi, num_qubits = 3, input_state = None):

    """
    Function generates the expected quantum state after simulating the n qubit XXX Heisenberg model for time t, on input state |0^n>

    Input: 
    - t, float, the time
    - num_qubits, int, the number of qubits
    -  input_state, 1D numpy array of length 2**num_qubits, the input state of the time evolution

    Output: numpy array of complex floats of size 2**num_qubits = 8, the output vector

    """

    if input_state is None:
        input_state = psi_0(num_qubits)

    return np.matmul(U_Heis(float(t), num_qubits) , input_state)















