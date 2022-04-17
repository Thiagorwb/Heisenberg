from ..compiler import *
from .gates import *

# Noiseless simulated backend
sim = QasmSimulator()

# Parameterize variable t to be evaluated at t=pi later
num_qubits = 3
basis = ['id', 'rz', 'sx', 'x', 'cx', 'reset']


#########################
### helper functions ####

def reveal(instructions, num_qubits):
    """
    Given circuit instructions and size, prints drawing of circuit and circuit parameters
    """
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    qc.append(instructions, qr)
    qc = qc.decompose().decompose()
    print(qc.decompose().draw())
    print(qc.decompose().count_ops())
    return
    

######################################
### single brick in decomposition ####

def basic_trotter(terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], time = np.pi, num_qubits = 3, qsd = False, name = "Trot"):
    """ 
    Constructs basic trotter decomposition for \exp[i H_heis t], given the decomposition
    Inputs: terms, list of strings, indicating gate sequence
            time, float
            num_qubits, int, number of qubits in H
            qsd, boolean, compiler flag
    Returns: QuantumCircuit instructions
    expected number of cxs is 12 without qsd, 6 with qsd
    """
    
    Trot_qr = QuantumRegister(num_qubits)
    Trot_qc = QuantumCircuit(Trot_qr, name=name)

    for index in terms:
        j, i = index[0], int(index[1]) - 1
        Trot_qc.append(QQ_gate(j, time = time), [Trot_qr[i], Trot_qr[i+1]])
        
    if qsd:
        Trot_qc = gen_qsd(Operator(Trot_qc).data)
        
    # Convert custom quantum circuit into a gate
    Trot_gate = Trot_qc.to_instruction()
    return Trot_qc


def ti_trotter(terms = None, time=np.pi, qsd = False, name = "Brick"):
    """
    Single trotter step using translational invariance - i.e. break apart terms on different qubits
    Inputs:
        time, float
        qsd, boolean, compiler flag
    Returns:
        quantum circuit instruction for single trotter step
    Expected number of 2 qubit gates with compiler flag: 3
    """

    Trot_qr = QuantumRegister(num_qubits)
    Trot_qc = QuantumCircuit(Trot_qr, name= name )

    # note, these terms commute, and we can break them up for free
    """ti_terms = ['x1', 'y1', 'z1']
                ti_evolution = basic_trotter(ti_terms, time, num_qubits = 2, qsd = qsd)
            """

    ti_evolution = gen_qsd(U_12(time).to_matrix())

    Trot_qc.append(ti_evolution, [Trot_qr[0], Trot_qr[1]])
    Trot_qc.append(ti_evolution, [Trot_qr[1], Trot_qr[2]])
    
    # Convert custom quantum circuit into a gate
    Trot_gate = Trot_qc.decompose()
    Trot_gate = Trot_qc.to_instruction()
    return Trot_qc

#######################################
### larger trotterization circuits ####

def basic_trotter_circuit(trotter_gate_function, trotter_steps, time = np.pi, terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], operational_qubits = [0, 1, 2], qsd = False):
    """ 
    Constructs Trotterized Circuit by repeating a single trotter step
    Inputs: trotter_gate_function, function, returns QuantumCircuit instructions 
            trotter_steps, int, number of repetitions of said single step
            time, total simulation time, assumed = np.pi
            terms, list of strings, indicating term decomposition in trotterization
    
    Returns: QuantumRegister, QuantumCircuit 
    """
    q1, q2, q3 = operational_qubits
        
    # Initialize quantum circuit for n qubits
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    # Simulate time evolution under H_heis3 Hamiltonian
    for _ in range(trotter_steps):
        qc.append(trotter_gate_function(terms = terms, time=time/trotter_steps, qsd = qsd), [qr[q1], qr[q2], qr[q3]])
        
    qc = qc.decompose()
    
    return qr, qc 

def second_order_trotter(trotter_steps, num_qubits = 7, time = np.pi, operational_qubits = [0, 1, 2],  qsd = False):
    """
    Hardcoding the second order trotterization for quality controle and optimizing edge cases
    """

    q1, q2, q3 = operational_qubits
        
    # Initialize quantum circuit for n qubits
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    if trotter_steps == 1:
        hot = ti_higher_order_trotter(order = 2, time = t, terms = None)
        instr = hot(time = time, qsd = False)
        qc.append(instr, [qr[q1], qr[q2], qr[q3]])
        return qr, qc
    
    # basic gates
    ttime = time / trotter_steps
    ti_terms = ['x1', 'y1', 'z1']
    ti_half_evolution = basic_trotter(ti_terms, ttime / 2, num_qubits = 2, qsd = qsd)
    ti_evolution = basic_trotter(ti_terms, ttime, num_qubits = 2, qsd = qsd)
     
    # circuit assembly
    qc.append(ti_half_evolution, [qr[q1], qr[q2]])

    for _ in range(trotter_steps - 1):
        qc.append(ti_evolution, [qr[q2], qr[q3]])
        qc.append(ti_evolution, [qr[q1], qr[q2]])
        
    qc.append(ti_evolution, [qr[q2], qr[q3]])
    qc.append(ti_half_evolution, [qr[q1], qr[q2]])
    qc = qc.decompose()
    
    return qr, qc   

def ti_higher_order_trotter(order = 2, time = np.pi, terms = None):
    '''
    Wrapper Function, for single step of the higher order trotter gate function leveraging translational invariance, to be repeated


    Inputs: order, int, order of higher order formula
            time, float, time of simulation of this trotter step
            terms, list of strings, indicating order of decomposition

    Returns: Function, which returns QuantumCircuit instructions
    '''

    if order == 1:
        return ti_trotter
    
    def wrapper(terms = terms, time = time, qsd = True):

        """
        Inner method - function generates single step of higher order trotterization, leveraging translational invariance
        Inputs:
        terms, list of strings, indicating order of decomposition
        time, float,
        qsd, boolean, compiler flag

        """
        
        num_qubits = 3
        Trot_qr = QuantumRegister(num_qubits)
        Trot_qc = QuantumCircuit(Trot_qr, name='Higher_Order_Trot')
        ti_terms = ['x1', 'y1', 'z1']

        coefficients = higher_trotter_coefficients(['1','2'], order, time)

        for term, coefficient in coefficients:
            i = int(term) - 1

            #ti_evolution = gen_qsd(U_12(coefficient).to_matrix())
            ti_evolution = basic_trotter(ti_terms, time = coefficient, num_qubits = 2, qsd = True)

            Trot_qc.append(ti_evolution, [Trot_qr[i], Trot_qr[i+1]])
            
        # Convert custom quantum circuit into a gate
        Trot_gate = Trot_qc.to_instruction()
        return Trot_gate
    
    return wrapper

def higher_order_trotter(order = 2, time = np.pi,  trotter_steps = 4):

    # generate coefficientss and parse out redundancies in composition
    all_coefficients = []
    coefficients = higher_trotter_coefficients(["2", "1"], order, time/trotter_steps)
    for _ in range(trotter_steps):
        all_coefficients += coefficients

    all_coefficients = hot_parse_coefficients(all_coefficients)

    # initialize circuit
    num_qubits = 3
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)

    for pair in all_coefficients:
        
        term, coefficient = pair

        # readout location
        i = int(term) -  1

        # append to circ
        qc.append(gen_qsd(U_12(coefficient).to_matrix()).to_instruction(), [i, i+1])

    return qr, qc

###################################
### Optimized 3 Qubit Circuits ####

def compiled_hot(order = 2, time = np.pi,  trotter_steps = 4, basis_gates = basis, qsd = True, optimization_level = 2):

    # compute higher order trotter
    qr_hot, qc_hot = higher_order_trotter(order = order, time = time, trotter_steps = trotter_steps)

    # compile and parse to geometry
    if qsd:
        nqc = gen_qsd(Operator(qc_hot.decompose()).data)
        gqcinstr =  geometric_compile(nqc.decompose(), 3, line(3))
    else:
        gqcinstr =  geometric_compile(qc_hot.decompose(), 3, line(3))
    geometric_qc = QuantumCircuit(qr_hot)
    geometric_qc.append(gqcinstr, qr_hot)
    geometric_qc = geometric_qc.decompose()

    # efficient transpiltion to basis gates

    if order < 6:
        level = 3
    else:
        level = optimization_level

    nqc = transpile(geometric_qc, basis_gates = basis_gates, optimization_level = level)

    """
    # check if transpilation is accurate - if not, do transpilation at lower optimization_level
    # given unitaries U1, U2, compute U1 U2^{\dagger} and check if multiple of identity
    matrix = np.matmul(Operator(nqc).data, linalg.inv(Operator(geometric_qc).data))
    if not np.allclose(matrix/matrix[0][0] , np.eye(8, 8)):
        nqc = transpile(geometric_qc, basis_gates = basis_gates, optimization_level = 0)

    """

    return qr_hot, nqc

################################################
### Generating Product Formula Coefficients ####

def higher_trotter_coefficients(terms, order, time):
    """ 
    Recursively computes the order and the exponents of the higher order trotter decomposition
    Implementation based on the paper The Theory of Trotter Error by Childs et al.
    Inputs: terms, list of strings, indicating term decomposition of trotterization
            order, int, order of product formula
            time, float
    Returns: list of tuples, (string, float)
    """
    
    coefficients = []
    
    if order == 1:
        for term in terms:
            coefficients.append((term, time))
        return coefficients
        
    if order == 2:
        for term in terms[:-1]:
            coefficients.append((term, time/2))
        
        coefficients.append((terms[-1], time))
        
        for term in reversed(terms[:-1]):
            coefficients.append((term, time/2))
            
        return coefficients
    
    u_order = 1/(4 - 4**(1/(order-1)))
    
    previous = higher_trotter_coefficients(terms, order - 2, u_order*time)
    middle = higher_trotter_coefficients(terms, order - 2, (1 - 4*u_order)*time)
    
    coefficients = previous + previous + middle + previous + previous
    
    return hot_parse_coefficients(coefficients)


def hot_parse_coefficients(coefficients):
    """
    Function parses through trotterization coefficients, removing redundancies

    """
    
    new_c = []
    current_term, current_coefficient = coefficients[0] 
     
    for term, coef in coefficients[1:]:
        if term != current_term:
            new_c.append([current_term, current_coefficient])
            current_term = term
            current_coefficient = coef
        else:
            current_coefficient += coef
            
    new_c.append([current_term, current_coefficient])
    return new_c

###############################
### BrickWork Architecture ####

def single_layer(time, num_qubits):
    """
    Function generates single layer of a Brickwork Architecture for the 1D heisenberg model
    Inputs:
    time, float, target time
    num_qubits, int,  number of qubits in system
    qsd, boolean, comppiler flag 

    Output:
    QuantumCircuit object

    """
    
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    # The Brick
    ti_terms = ['x1', 'y1', 'z1']
    ti_evolution = gen_qsd(U_12(time).to_matrix())
    
    # even terms
    for i in range(0, num_qubits - 1, 2):
        qc.append(ti_evolution, [i, i+1])
        
    # odd terms
    for i in range(1, num_qubits - 1, 2):
        qc.append(ti_evolution, [i, i+1])
    
    return qc

def brickwork(time, num_qubits, trotter_steps = 1, qsd = True):

    """
    Function generates Brickwork Architecture circuit for time evolution of the 1D heisenberg model

    Inputs:
    time, float, target time
    num_qubits, int,  number of qubits in system
    trotter_steps, int, number of trotter steps in decomposition
    qsd, boolean, comppiler flag 

    Output:
    QuantumCircuit object
    """
    
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)

    if num_qubits == 2:
        qc.append(gen_qsd(U_12(time).to_matrix()).to_instruction(), qr)
        return qr, qc
    
    layer = single_layer(time/trotter_steps, num_qubits)
    
    for _ in range(trotter_steps):
        qc.append(layer, qr)
        
    qc = qc.decompose()
        
    return qr, qc



    
        
