from .libraries import *

"""
Here we place all matrix decomposition and linear algebra helper code
"""

def matrix_decomposition(M, num_qubits = 3):
    
    paulis = [I, X, Y, Z]
    letters = ['I', "X", "Y", "Z"]
    coefficients = []
    
    for A, a in zip(paulis, letters):
        for B, b in zip(paulis, letters):
            for C, c in zip(paulis, letters):
                
                cf = np.trace(np.matmul((A^B^C).to_matrix(), M)) /(2**num_qubits)
                coefficients.append([[a, b, c], cf])
                    
    return coefficients

def coef_to_operator(coefficients,  num_qubits = 3):
    
    pauli_sum = Operator(np.zeros((2**num_qubits, 2**num_qubits)))
    
    for letters, cf in coefficients:
        label = "".join(letters)
        pauli = Operator(Pauli(label = label)) 
        pauli_sum += pauli * cf 
        
    return pauli_sum

def diagonalize(M):
    
    d, V = np.linalg.eig(M)
    
    return V, d

def is_unitary(m):
    return np.allclose(np.eye(len(m)), m.dot(m.T.conj()))

def schur_decomposition(U):
    
    D, Z = linalg.schur(U, output = 'complex')
    return Z, D