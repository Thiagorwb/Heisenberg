from .testing import *


###################
#### Compiler Exp #
###################

def gate_count_experiment(max_order, min_steps, max_steps, time = np.pi, backend = sim):

    orders = [1] + [i for i in range(2, max_order + 1, 2)]

    compiled_counts = []
    uncompiled_counts = []

    for order in orders:

        order_ccounts = []
        order_unccounts = []
        print("Starting kth order formulas, 2k = ", order)

        for L in range(min_steps, max_steps + 1):

            compiled_count = test_gate_count(order, time, L, backend = backend, compiler = True)
            uncompiled_count = test_gate_count(order, time, L, backend = backend, compiler = False)

            order_ccounts.append(compiled_count)
            order_unccounts.append(uncompiled_count)

            print("at number of steps ", L, "(compiled, uncompiled) CNOT counts are:", compiled_count, uncompiled_count)

        compiled_counts.append(order_ccounts)
        uncompiled_counts.append(order_unccounts)

    return compiled_counts, uncompiled_counts


#################################
#### 3 Qubit Trotterization Exp #
#################################

def product_experiment(time, order = 2, min_steps = 4, max_steps = 10, backend = backend, repetitions = 3, compiler = True, verbose = False):

    """
    
    Runs the kth order product formula circuit for a range of trotter steps. Outputs the average fidelity at each step.

    Inputs: 
        time, float, 
        order, int, the order of the product formula = 2k
        min_steps, int, min number of trotter_steps
        max_steps, int, max number of trotter_steps
        backend, backend obj, the IBMQ backend
        repetitions, # different job queues
        compiler, bool, compiler flag
        verbose, bool
    Outputs:
        list of floats, the fidelities
    
    """

    fids = []

    for L in range(min_steps, max_steps+1):

        print(L , "Trotterization steps")
        fid = test_formula(order = order, trotter_steps = L, time = time, backend = backend, repetitions = repetitions, verbose = verbose, compiler = compiler)
        fids.append(fid)

    return fids

def time_experiment(order, trotter_steps ,  max_time = np.pi, time_steps = 10, backend = sim,  repetitions = 5, compiler = True):

    ts = [ i * max_time/time_steps for i in range(1, time_steps + 1)]
    fids = []
    probs  = []
    
    for t in ts:
        print("Simulating until time ", t)
        fid, prob = test_formula(order = order, trotter_steps = trotter_steps, time = t, backend = backend, repetitions = repetitions, verbose = False, return_prob = True)

        fids.append(fid)
        probs.append(prob)
        
    return ts, fids, probs

###############################
#### Experiments for larger n #
###############################

def brick_steps_experiment(num_qubits, time = np.pi, min_steps = 1, max_steps = 10, backend = backend, repetitions = 3, qsd = False, verbose = False):

    fids = []
    for L in range(min_steps, max_steps + 1):

        print(L , "Trotterization steps")

        f = test_brickwork(num_qubits, time, trotter_steps = L, backend = backend, repetitions = repetitions, qsd = qsd, verbose = verbose)

        fids.append(f)

    print("\n")

    return fids

def size_experiment(max_size, min_size = 3, time = np.pi, min_steps = 4, max_steps = 10, backend = backend, repetitions = 3, qsd = False, verbose = False):

    size_fids = []
    for size in range(min_size, max_size + 1):

        fids = brick_steps_experiment(size, time = np.pi, min_steps = min_steps, max_steps = max_steps, backend = backend, repetitions = repetitions, qsd = qsd, verbose = verbose)

        nsteps = np.argmax(fids)
        max_fid = np.max(fids)

        print("The fidelity for the ", size, "qubit experiment was ", max_fid, " at ", min_steps + nsteps, "trotter steps")
        size_fids.append(max_fid)

    return size_fids


################
## Unfinished ##
################

"""def depth_experiment(num_qubits, min_size = 3, time = np.pi, min_steps = 4, max_steps = 10, backend = sim, repetitions = 3, qsd = False, verbose = False):
        
    # initialize backend circuit geometry 
    if backend.configuration().backend_name == 'qasm_simulator':
        system_size = num_qubits
        graph = clique(system_size)
    else:
        system_size = backend.configuration().num_qubits
        graph = backend.configuration().coupling_map

    counts = []

    for nsteps in range(min_steps, max_steps):
    
        # compute brickwork circuit 
        qr1, qc_brickwork = brickwork(time, num_qubits, trotter_steps = nsteps, qsd = True)

        # compile and parse to circuit geometry
        if qsd:
            geometric_instr = geometric_compile(gen_qsd(Operator(qc_brickwork).data).decompose(), system_size, graph)
        else:
            geometric_instr = geometric_compile(qc_brickwork.decompose(), system_size, graph)

        qr = QuantumRegister(system_size)
        qc.append(geometric_instr, [range(num_qubits)])

        counts.append(qc.count_ops().get("cx", 0))

    return counts
"""




































