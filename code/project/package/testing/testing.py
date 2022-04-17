from .tomography import *



######################
## Testing Routines ##
######################

def test_formula(order = 4, trotter_steps = 10, time = np.pi, backend = sim, repetitions = 3, shuttle = False, verbose = False, compiler = True, return_prob = False):


    ####################################
    # initialize backend and qubits
    is_jakarta = False
    
    if backend.configuration().backend_name == 'qasm_simulator':
        system_size = 3
        graph = clique(system_size)
        basis = None
    else:
        system_size = backend.configuration().num_qubits
        basis = backend.configuration().basis_gates

    # pick operational qubits
    q1, q2, q3 = 0, 1, 2
    if backend.configuration().backend_name == 'ibmq_jakarta' or backend.configuration().backend_name == "qasm_simulator(ibmq_jakarta)":
        q1, q2, q3 = 1, 3, 5
        is_jakarta = True
        
    if backend.configuration().backend_name == 'ibmq_quito' or backend.configuration().backend_name == "qasm_simulator(ibmq_quito)":
        q1, q2, q3 = 0, 1, 2
        
    if backend.configuration().backend_name == 'ibmq_manila' or backend.configuration().backend_name == "qasm_simulator(ibmq_manila)":
        q1, q2, q3 = 2, 3, 4
        

    #################################
    # generate formula circuit

    # initialize register
    qr = QuantumRegister(system_size)

    # compute higher order trotter
    if compiler:
        qr_hot, qc_hot = compiled_hot(order = order, time = time,  trotter_steps = trotter_steps, basis_gates = basis)

    else:
        qr_hot, qc_hot = higher_order_trotter(order = order, time = time, trotter_steps = trotter_steps)

    instructions = qc_hot.decompose().to_instruction()
    
    #################################
    # construct tomography circuit
    qc = QuantumCircuit(qr)

    if shuttle and is_jakarta:
        # initialize 011 state
        qc.x([qr[5], qr[4]])
        qc.append(instructions, [qr[3], qr[5], qr[4]])
        qc.append(swap_gate(), [qr[1], qr[3]])
        qc.append(swap_gate(), [qr[3], qr[5]])
        qc.append(swap_gate(), [qr[5], qr[4]])
    else:
        # initialize 011 state
        qc.x([qr[q2], qr[q3]])
        qc.append(instructions, [qr[q1], qr[q2], qr[q3]])

    # compile into native basis gates for simulation and consistency
    qc2 =  QuantumCircuit(qr)
    qc2.append(transpile(qc, backend, optimization_level = 0), qr)
        
    if verbose:
        print("Circuit Description")
        print(qc2.decompose().count_ops())
        print(qc2.decompose().draw())
        
    # tomography circuit
    st_qcs  =  state_tomography_circuits(qc2, [qr[q1], qr[q2], qr[q3]])
        
    # assign jobs and evaluate fidelity
    jobs = evaluate_circuit(st_qcs, reps = repetitions, backend = backend, verbose = verbose)
    
    result_state = prepare_state(time)
    fid = average_fidelity(jobs, st_qcs, result_state)

    if return_prob:

        prob = average_prob(jobs, st_qcs, psi_0(3))
        return fid, prob
    
    
    return fid


################# deprecated
def test_3_qubits(instructions, time = np.pi, backend = sim, repetitions = 3, shuttle = False, verbose = False, return_prob = False):
    
    """
    
    Evaluates a higher order product formula circuit on 3 Qubits and computes the fidelity to the time evolution of the 3-Qubit XXX Heisenberg model

    Inputs:
    - instructions, QuantumCircuit instructions, 
    - time, float,
    - backend, Backend obj
    - repetitions, 
    - verbose,  bool
    - return_prob, bool, flag if true returns time-experiment information
    - shuttle, bool, heurstic flag

    Ouputs:
    - float, the fidelity to the time evolution of the 3-Qubit XXX Heisenberg model

    """

    is_jakarta = False
    
    if backend.configuration().backend_name == 'qasm_simulator':
        num_qubits = 3
        graph = clique(num_qubits)
    else:
        num_qubits = backend.configuration().num_qubits

    # pick operational qubits
    q1, q2, q3 = 0, 1, 2
    if backend.configuration().backend_name == 'ibmq_jakarta' or backend.configuration().backend_name == "qasm_simulator(ibmq_jakarta)":
        q1, q2, q3 = 1, 3, 5
        is_jakarta = True
        
    if backend.configuration().backend_name == 'ibmq_quito' or backend.configuration().backend_name == "qasm_simulator(ibmq_quito)":
        q1, q2, q3 = 0, 1, 2
        
    if backend.configuration().backend_name == 'ibmq_manila' or backend.configuration().backend_name == "qasm_simulator(ibmq_manila)":
        q1, q2, q3 = 2, 3, 4
        
    # initialize circuit
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    if shuttle and is_jakarta:
        # initialize 011 state
        qc.x([qr[5], qr[4]])
        qc.append(instructions, [qr[3], qr[5], qr[4]])
        qc.append(swap_gate(), [qr[1], qr[3]])
        qc.append(swap_gate(), [qr[3], qr[5]])
        qc.append(swap_gate(), [qr[5], qr[4]])
    else:
        # initialize 011 state
        qc.x([qr[q2], qr[q3]])
        qc.append(instructions, [qr[q1], qr[q2], qr[q3]])
    
    # compile into native basis gates
    qc2 =  QuantumCircuit(qr)
    qc2.append(transpile(qc, backend, optimization_level = 0), qr)
        
    if verbose:
        print("Circuit Description")
        print(qc2.decompose().count_ops())
        print(qc2.decompose().draw())
        
    # tomography circuit
    st_qcs  =  state_tomography_circuits(qc2, [qr[q1], qr[q2], qr[q3]])
        
    # assign jobs and evaluate fidelity
    jobs = evaluate_circuit(st_qcs, reps = repetitions, backend = backend, verbose = verbose)
    
    result_state = prepare_state(time)
    fid = average_fidelity(jobs, st_qcs, result_state)

    if return_prob:

        prob = average_prob(jobs, st_qcs, psi_0(3))
        return fid, prob
    
    return fid

def test_brickwork(num_qubits, time = np.pi, trotter_steps = 5, backend = sim, repetitions = 3, qsd = False, verbose = True):

    """
    Evaluates the brick-work architecture and computes the fidelity to the time evolution of the n-Qubit XXX Heisenberg model

    Inputs:
    - num_qubits, int, the intended system size
    - time, float,
    - trotter_steps, int, number of formula steps
    - backend, Backend obj
    - repetitions, 
    - verbose,  bool
    - qsd, compiler flag (keep off for large system sizes)

    Ouputs:
    - float, the fidelity to the time evolution of the 3-Qubit XXX Heisenberg model

    """
    
    # initialize backend circuit geometry 
    if backend.configuration().backend_name == 'qasm_simulator':
        system_size = num_qubits
        graph = clique(system_size)
    else:
        system_size = backend.configuration().num_qubits
        graph = backend.configuration().coupling_map
    
    # compute brickwork circuit 
    qr1, qc_brickwork = brickwork(time, num_qubits, trotter_steps = trotter_steps, qsd = True)

    # compile and parse to circuit geometry
    if qsd:
        geometric_instr = geometric_compile(gen_qsd(Operator(qc_brickwork).data).decompose(), system_size, graph)
    else:
        geometric_instr = geometric_compile(qc_brickwork.decompose(), system_size, graph)
    
    # initialize circuit
    qr = QuantumRegister(system_size)
    qc = QuantumCircuit(qr)
    qc.x([qr[num_qubits-2], qr[num_qubits-1]])

    # add instructions
    qc.append(geometric_instr, qr)
    nqc = QuantumCircuit(qr)

    # transpile to circuit for consistency
    nqc.append(transpile(qc, backend, optimization_level = 2), qr)

    if verbose:
        print("The Initial Brickwork Circuit", qc_brickwork.decompose().count_ops())

        # to draw, uncomment below:
        #print(qc_brickwork.draw())

        print("The Compiled Circuit", nqc.decompose().count_ops())

        # to draw, uncomment below:
        print(nqc.decompose().draw())
    
    # tomography circuit
    st_qcs  =  state_tomography_circuits(nqc, [qr[i] for i in range(num_qubits)])
    
    # assign jobs and evaluate fidelity
    jobs = evaluate_circuit(st_qcs, reps = repetitions, backend = backend, verbose = verbose)
    
    result_state = prepare_system(time, num_qubits)
    fid = average_fidelity(jobs, st_qcs, result_state)
    
    return fid


def test_gate_count(order, time, trotter_steps, backend = sim, compiler = False):

    if backend.configuration().backend_name == 'qasm_simulator':
        basis = None
    else:
        basis = backend.configuration().basis_gates

    if compiler:
        qr_hot, qc_hot = compiled_hot(order = order, time = time,  trotter_steps = trotter_steps, basis_gates = basis)

    else:
        qr_hot, qc_hot = higher_order_trotter(order = order, time = time, trotter_steps = trotter_steps)

    return qc_hot.decompose().count_ops().get('cx', 0)



    
   


