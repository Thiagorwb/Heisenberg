from .libraries import *

"""
Here we place all matrix decomposition and linear algebra helper code
"""

def diagonalize(M):
    """
    Diagonalizes matrix M, returns unitary V and eigenvalues
    Input:
    -  M, 2D numpy array,
    Output:
    - V, 2D numpy array, unitary matrix
    - d, 1D numpy array, eigenvalues
    """
    
    d, V = np.linalg.eig(M)
    
    return V, d

def is_unitary(M):
    """
    Checks if matrix M is unitary by computing MM^T and comparing with identity
    Input:
    -  M, 2D numpy array,
    Returns:
    - boolean, True if close to unitary up to machine precision
    """

    return np.allclose(np.eye(len(M)), M.dot(M.T.conj()))

def schur_decomposition(U):
    """
    Computes Schur decomposition of matrix U, often a unitary matrix in our applications

    U = Z D Z^\dagger, where

    Input:
    - U, 2D numpy array,
    Output:
    - U, 2D numpy array, upper triangular complex matrix
    - Z, 2D numpy array, unitary matrix

    """
    
    D, Z = linalg.schur(U, output = 'complex')
    return Z, D


########################################
# NOTE: DEPRECATED - DONT USE BELOW LINE

def matrix_decomposition(M, num_qubits = 3):
    """
    Function decomposes a generic 2**3 by 2**3 matrix into Pauli basis, and returns coefficients
    """
    
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

