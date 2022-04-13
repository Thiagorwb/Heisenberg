from .gates import *


# Noiseless simulated backend
sim = QasmSimulator()

# Parameterize variable t to be evaluated at t=pi later
t = Parameter('t')
num_qubits = 3

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
    

def basic_trotter(terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], time = t, num_qubits = 3, backend = sim, compiler = False, qsd = False, coupling_map = None):
    """ 
    Constructs basic trotter decomposition for \exp[i H_heis t]
    Inputs: terms, list of strings, indicating gate sequence
            time, float
    Returns: QuantumCircuit instructions
    expected number of cxs is 12
    """
    
    Trot_qr = QuantumRegister(num_qubits)
    Trot_qc = QuantumCircuit(Trot_qr, name='Trot')
    

    for index in terms:
        j, i = index[0], int(index[1]) - 1
        Trot_qc.append(QQ_gate(j, time = time), [Trot_qr[i], Trot_qr[i+1]])
        
    if compiler:
        Trot_qc = transpile(Trot_qc, basis_gates = backend.configuration().basis_gates, optimization_level = 3)
        
    if qsd:
        Trot_qc = gen_qsd(Operator(Trot_qc).data)
        
    # Convert custom quantum circuit into a gate
    Trot_gate = Trot_qc.to_instruction()
    return Trot_qc


def ti_trotter(terms = None, time=t, backend = sim, compiler = False, qsd = False, coupling_map = None):
    Trot_qr = QuantumRegister(num_qubits)
    Trot_qc = QuantumCircuit(Trot_qr, name='Trot')
    
    ti_terms = ['x1', 'y1', 'z1']
    ti_evolution = basic_trotter(ti_terms, time, num_qubits = 2, backend=backend, compiler = compiler, qsd = qsd, coupling_map = coupling_map)
    Trot_qc.append(ti_evolution, [Trot_qr[0], Trot_qr[1]])
    Trot_qc.append(ti_evolution, [Trot_qr[1], Trot_qr[2]])
    
    if compiler:
        Trot_qc = transpile(Trot_qc, basis_gates = backend.configuration().basis_gates, optimization_level=3)

    # Convert custom quantum circuit into a gate
    Trot_gate = Trot_qc.decompose()
    Trot_gate = Trot_qc.to_instruction()
    return Trot_qc

def basic_trotter_circuit(trotter_gate_function, trotter_steps, num_qubits = 7, time = t, terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], operational_qubits = [0, 1, 2], backend = sim, compiler = False, qsd = False, coupling_map = None):
    """ 
    Constructs Trotter Circuit for Jakarta machine by Repeating a given gate
    Inputs: trotter_gate_function, function, returns QuantumCircuit instructions 
            trotter_steps, int
            target_time, final simulation time, assumed = np.pi
            terms, list of strings, indicating gate sequence
    
    Returns: QuantumRegister, QuantumCircuit 
    """
    q1, q2, q3 = operational_qubits
        
    # Initialize quantum circuit for n qubits
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    # Simulate time evolution under H_heis3 Hamiltonian
    for _ in range(trotter_steps):
        qc.append(trotter_gate_function(terms = terms,time=time/trotter_steps, backend = backend, compiler = compiler, qsd = qsd, coupling_map = coupling_map), [qr[q1], qr[q2], qr[q3]])
        
    # Evaluate simulation at target_time (t=pi) meaning each trotter step evolves pi/trotter_steps in time
    qc = qc.decompose()
    if compiler:
        qc = transpile(qc, optimization_level = 3, basis_gates  = backend.configuration().basis_gates)
    
    #qc = qc.bind_parameters({t: target_time/trotter_steps})
    
    return qr, qc 

def second_order_trotter(trotter_steps, num_qubits = 7, time = t, operational_qubits = [0, 1, 2], backend = sim, compiler = False, qsd = False, coupling_map = None):
    q1, q2, q3 = operational_qubits
        
    # Initialize quantum circuit for n qubits
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    if trotter_steps == 1:
        hot = ti_higher_order_trotter(order = 2, time = t, terms = None)
        instr = hot(time = time, backend = backend, compiler = False, qsd = False, coupling_map = None)
        qc.append(instr, [qr[q1], qr[q2], qr[q3]])
        return qr, qc
    
    # basic gates
    ttime = time / trotter_steps
    ti_terms = ['x1', 'y1', 'z1']
    ti_half_evolution = basic_trotter(ti_terms, ttime / 2, num_qubits = 2, backend=backend, compiler = compiler, qsd = qsd)
    ti_evolution = basic_trotter(ti_terms, ttime, num_qubits = 2, backend=backend, compiler = compiler, qsd = qsd)
     
    # circuit assembly
    qc.append(ti_half_evolution, [qr[q1], qr[q2]])

    for _ in range(trotter_steps - 1):
        qc.append(ti_evolution, [qr[q2], qr[q3]])
        qc.append(ti_evolution, [qr[q1], qr[q2]])
        
    qc.append(ti_evolution, [qr[q2], qr[q3]])
    qc.append(ti_half_evolution, [qr[q1], qr[q2]])
    qc = qc.decompose()
    
    return qr, qc   

def ti_higher_order_trotter(order = 2, time = t, terms = None):
    def wrapper(terms = terms, time = time, backend = sim, compiler = False, qsd = False, coupling_map = None):
        
        num_qubits = 3
        Trot_qr = QuantumRegister(num_qubits)
        Trot_qc = QuantumCircuit(Trot_qr, name='Higher_Order_Trot')
        ti_terms = ['x1', 'y1', 'z1']

        coefficients = higher_trotter_coefficients(['1','2'], order, time)

        for term, coefficient in coefficients:
            i = int(term) - 1
            ti_evolution = basic_trotter(ti_terms, time = coefficient, num_qubits = 2, qsd = qsd)
            Trot_qc.append(ti_evolution, [Trot_qr[i], Trot_qr[i+1]])
            
        if compiler:
            Trot_qc = transpile(Trot_qc, basis_gates = backend.configuration().basis_gates, optimization_level=3)

        # Convert custom quantum circuit into a gate
        Trot_gate = Trot_qc.to_instruction()
        return Trot_gate
    
    return wrapper

def higher_order_trotter(order = 2, time = t, random = False, terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2']):
    '''
    Generates higher order trotter gate function.
    Inputs: order, int, order of higher order formula
            time, float, time of simulation of this trotter step
            random, boolean, random choice of term ordering?
            terms, list of strings, indicating order of decomposition
    Returns: Function, which returns QuantumCircuit instructions
    '''
    
    def wrapper(terms = terms, time = time,  backend = sim, compiler = False, qsd = False, coupling_map = None):

        num_qubits = 3
        Trot_qr = QuantumRegister(num_qubits)
        Trot_qc = QuantumCircuit(Trot_qr, name='TI Higher_Order_Trot')

        if random:
            terms = np.random.permutation(['x1', 'y1', 'z1', 'x2', 'y2', 'z2'])

        coefficients = higher_trotter_coefficients(terms, order, time)

        for term, coefficient in coefficients:
            j, i = term[0], int(term[1]) - 1
            Trot_qc.append(QQ_gate(j, coefficient), [Trot_qr[i], Trot_qr[i+1]])
            
        if compiler:
            Trot_qc = transpile(Trot_qc, basis_gates = backend.configuration().basis_gates, optimization_level=3)

        # Convert custom quantum circuit into a gate
        Trot_gate = Trot_qc.to_instruction()
        return Trot_gate
    
    
    return wrapper

def higher_trotter_coefficients(terms, order, time):
    """ 
    Recursively computes the order and the exponents of the higher order trotter decomposition
    Implementation based on the paper The Theory of Trotter Error [need citation]
    Inputs: terms, list of strings, indicating gate sequence
            order, int, order of formula
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
    
        
