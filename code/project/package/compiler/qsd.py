from ..circuits import *
from .gray import *

def de_multiplex(a, b):
    '''
    Input: unitary matrices a, b of same shape
    Output: unitary matrices V, W, diagonal matrix d, s.t.
            a = V d W, b = V d^\dagger W
    '''
    
    ab = np.matmul(a, np.matrix(b).getH())
    
    V, T = schur_decomposition(ab)
    d = np.array([np.sqrt(si) for si in np.diagonal(T)])
    W = np.matmul(np.matmul(np.diag(d), np.matrix(V).getH()), b)
        
    return V, np.diag(d), W


def multiplex_diagonal(matrix):
    """
    Input: diagonal, unitary matrix
    Returns: quantum circuit for [matrix, 0], [0, matrix^\dagger]
            Using only 2 qubit gates
    
    """
    
    if not np.all(matrix == np.diag(np.diagonal(matrix))) or not is_unitary(matrix):
        print("Input Matrix is not diagonal!")
        
    dim, dim = np.shape(matrix)
    num_qubits = int(np.log2(dim)) + 1
    
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    # generate angles
    angles = np.array([(np.log(d)/1.j).real for d in np.diagonal(matrix)])
    g_angles = gray_angles(angles)
    
    for i in range(dim):
        
        # rotate last qubit
        angle = g_angles[i]
        if angle != 0.0:
            qc.rz(- 2 * angle, qr[-1])
        
        if i == dim - 1:
            index = num_qubits - 2
        else:
            index = gray_index(i)
               
        """# uncomment here
        if coupling_map:
            lrcx = long_range_cx(num_qubits - 1, index, num_qubits)
            qc.append(lrcx, qr)
        else:
            qc.cx(qr[index], qr[-1])"""
        
        qc.cx(qr[index], qr[-1])
        
    return qc.to_instruction()   

def gray_angles(theta):
    """
    Input: vector theta of angles of controlled rotation gate
    Output: rotated angles after gray code optimization
    
    """
    
    dim = len(theta)
    n_bits = int(np.log2(dim))
    
    M = np.ones((dim, dim))
    for i in range(dim):
        for j in range(dim):
            exp = binary(i, n_bits).dot(binary(gray_code(j), n_bits))
            M[i][j] = (-1)**exp    
    
    M1 = linalg.inv(M)
    new_angles = np.dot(M1, theta)
    
    return new_angles
    

def cs_circuit(theta, num_qubits):
    """
    Input: vector theta of angles of multiplexed ry gate
    Output: 2-qubit gate circuit for multiplexed ry
    """
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    # compute angles after gray code transformation
    g_angles = gray_angles(theta)
    
    for i in range(len(theta)):
        
        # rotate last qubit
        angle = g_angles[i]
        if angle != 0.0:
            qc.ry( 2 * angle, qr[-1])
        
        if i == len(theta) - 1:
            index = num_qubits - 2
        else:
            index = gray_index(i)
                    
        """if coupling_map:
            lrcx = long_range_cx(num_qubits - 1, index, num_qubits)
            qc.append(lrcx, qr)
        else:
            qc.cx(qr[index], qr[-1])"""
        
        qc.cx(qr[index], qr[-1])
                   
    return qc.to_instruction()


def zyz_gate(unitary):
    
    # check if unitary and correct size
    if len(unitary) != 2 :
        print("Cannot decompose larger gates with ZYZ decomposition !")
    if not is_unitary(unitary):
        print("Matrix is not Unitary")
    
    qr = QuantumRegister(1)
    qc = QuantumCircuit(qr)
    
    zyzdecomposer = OneQubitEulerDecomposer("ZYZ")
    a, b, c, d = zyzdecomposer.angles_and_phase(unitary)
    
    if a == 0.0:
        if b + c != 0.0:
            qc.rz(b + c, qr)
    
    else:
        if c != 0.0:
            qc.rz(c, qr)
        qc.ry(a, qr)
        if b != 0.0: 
            qc.rz(b, qr)
    
    #qc = transpile(qc, optimization_level = 3, basis_gates = jbasis)
    
    return qc.decompose().to_instruction()

def physical_zyz(unitary):
    
    # check if unitary and correct size
    if len(unitary) != 2 :
        print("Cannot decompose larger gates with ZYZ decomposition !")
    if not is_unitary(unitary):
        print("Matrix is not Unitary")
        
    qr = QuantumRegister(1)
    qc = QuantumCircuit(qr)
    
    zyzdecomposer = OneQubitEulerDecomposer("ZYZ")
    a, b, c, d = zyzdecomposer.angles_and_phase(unitary)
    
    if a == 0.0:
        if b + c != 0.0:
            qc.rz(b + c, qr)
    
    else:
        if c != 0.0:
            qc.rz(c, qr)
            
        qc.sx(qr) 
        
        if a != -np.pi and a != np.pi:
            qc.rz(a + np.pi, qr)
            
        qc.sx(qr) 
        qc.rz(b + 3*np.pi, qr)
        
    return qc.to_instruction()
    

def kak1_gate(unitary):
    
    dim, dim = np.shape(unitary)
    if dim == 2:
        return physical_zyz(unitary)
    
    # check if unitary and correct size
    if len(unitary) != 4 :
        print("Cannot decompose larger gates with ZYZ decomposition !")
    if not is_unitary(unitary):
        print("Matrix is not Unitary")
        
    kak1_decomposer = TwoQubitBasisDecomposer(CX, euler_basis = 'ZYZ')
    return kak1_decomposer(unitary).decompose().to_instruction()
    
def QSD(U):
    
    dim, dim = np.shape(U)
    num_qubits = int(np.log2(dim))
    
    if dim == 2:
        return physical_zyz(U)
    
    if dim == 4:
        return kak1_gate(U)
    
    # begin with cos sin decomposition, and separate block unitaries
    u_multi, cs, v_multi = linalg.cossin(U, p = dim / 2, q = dim / 2 , separate = True)
    u1, u2 = u_multi
    v1, v2 = v_multi
    
    # next, de-multiplex each side by diagonalization
    uV, uD, uW = de_multiplex(u1, u2)
    vV, vD, vW = de_multiplex(v1, v2)
        
    # recursively generate circuits for multiplexers u, v, via the diagonal multiplexer
    qc_vV = QSD(vV)
    qc_vD = multiplex_diagonal(vD)
    qc_vW = QSD(vW)
    
    qc_uV = QSD(uV)
    qc_uD = multiplex_diagonal(uD)
    qc_uW = QSD(uW)
    
    # generate multiplexed Ry rotation for CS matrix
    qc_cs = cs_circuit(cs, num_qubits)
    
    # Put circuit together
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    # right multiplexer
    qc.append(qc_vW, [qr[i] for i in range(0, num_qubits - 1)])
    qc.append(qc_vD, qr)
    qc.append(qc_vV, [qr[i] for i in range(0, num_qubits - 1)])
        
    # cossine sin circuit
    qc.append(qc_cs, qr)
    
    # Left multiplexer
    qc.append(qc_uW, [qr[i] for i in range(0, num_qubits - 1)])
    qc.append(qc_uD, qr)
    qc.append(qc_uV, [qr[i] for i in range(0, num_qubits - 1)])
    
    #line = [[i, i+1] for i in range(num_qubits - 1)]
    #qc = transpile(qc.decompose(), optimization_level = 3, approximation_degree = 1, basis_gates = jbasis, coupling_map = line)
    
    return qc.decompose().to_instruction()

def gen_qsd(U):
    dim, dim = np.shape(U)
    num_qubits = int(np.log2(dim))
    
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    qc.append(QSD(U), qr)
    
    return qc.decompose().decompose()

