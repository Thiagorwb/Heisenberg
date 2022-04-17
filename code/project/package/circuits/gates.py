# import libraries
from ..helpers import *

# Parameterize variable t to be evaluated at t=pi later
t = Parameter('t')
num_qubits = 3

################
# Simple gates #
################

def swap_gate():
    """
    Returns the quantum circuit instructions of a swap gate
    """

    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr)
    qc.cx(qr[0], qr[1])
    qc.cx(qr[1], qr[0])
    qc.cx(qr[0], qr[1])
    return qc.to_instruction()

def QQ_gate(pauli_name, time = t, name = None):
    """
    Function returns a circuit instruction for \exp[i Q \otimes Q t], for a given Pauli operator Q
    Inputs: pauli_name, string, \in 'x', 'y', 'z'
            time, float
            name, string, labels gate
    Returns: QuantumCircuit instructions
    Expected number of 2-qubit gates: 2
    """
    
    qr = QuantumRegister(2)
    qc = QuantumCircuit(qr, name= name)
    
    if pauli_name == 'x':
        qc.ry(np.pi/2,[0,1])
    elif pauli_name =='y':
        qc.rx(np.pi/2,[0,1])
    
    qc.cnot(0,1)
    qc.rz(2 * time, 1)
    qc.cnot(0,1)
    
    if pauli_name == 'x':
        qc.ry(- np.pi/2,[0,1])
    elif pauli_name == 'y':
        qc.rx(- np.pi/2,[0,1])

    # Convert custom quantum circuit into a gate
    return qc.to_instruction()

##########################
# Gates for the Compiler #
##########################

def long_range_cx(target, controle, num_qubits):
    """
    Function generates the "long range CNOT" quantum circuit instructions for a CNOT gate between two qubits
    which are not adjacent in the underlying circuit geometry. 

    An underlying chain arrangement is assumed for the qubits.

    Inputs:
    target, int, target qubit index
    controle, int, controle qubit index
    num_qubits, int, number of qubits

    Returns: QuantumCircuit instructions

    """
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    if np.abs(controle - target) == 1:
        qc.cx(controle, target)
        return qc.to_instruction()
    
    if controle < target:
        
        for i in range(controle, target):
            qc.cx(i, i + 1)
        for i in range(target - 2, controle, -1):
            qc.cx(i, i + 1)
        for i in range(controle, target):
            qc.cx(i, i + 1)
        for i in range(target - 2, controle, -1):
            qc.cx(i, i + 1)
            
    if controle > target:
        
        for i in range(controle, target, -1):
            qc.cx(i, i - 1)
        for i in range(target + 2, controle):
            qc.cx(i, i - 1)
        for i in range(controle, target, -1):
            qc.cx(i, i - 1)
        for i in range(target + 2, controle):
            qc.cx(i, i - 1)
            
    return qc.to_instruction()
        
def cx_on_a_path(path, num_qubits):
    """
    Function generates the "long range CNOT" quantum circuit instructions for a CNOT gate between two qubits
    which are not adjacent in the underlying circuit geometry. 

    Inputs:
    path, list of ints, the list of qubit indices from controle to target
    num_qubits, int, number of qubits on chio

    Returns: QuantumCircuit instructions for effective CNOT supported on the path
    """
    
    controle, target = path[0], path[-1]
    dist = len(path) - 1
    
    # generate circuit
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    if dist == 1:
        qc.cx(controle, target)
        return qc.to_instruction()
    
    for i in range(dist):
        qc.cx(path[i], path[i + 1])
    for i in range(dist - 2, 0, -1):
        qc.cx(path[i], path[i + 1])   
    for i in range(dist):
        qc.cx(path[i], path[i + 1])
    for i in range(dist - 2, 0, -1):
        qc.cx(path[i], path[i + 1])   
            
    return qc.to_instruction()
        
