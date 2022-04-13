from ..helpers import *
from ..ciruits import *
from ..compiler import *


###### Evaluation and Tomography ###########

def evaluate_circuit( st_qcs, backend = sim, shots = 5000, reps = 8, verbose = False):
    """
    Repeated evaluation and tomography of quantum circuit
    Inputs: st_qcs, state tomography circuit
            backend, 
            shots, number of repetitions for state tomography
            reps, number of overall evaluations
    """

    jobs = []
    for _ in range(reps):
        # execute
        #if noise_model:
        #    job = execute(st_qcs, backend, shots=shots, noise_model = noise_model)
        #else:
        job = execute(st_qcs, backend, shots=shots)
        if verbose:
            print('Job ID', job.job_id())
        jobs.append(job)
        
    return jobs

# Compute the state tomography based on the st_qcs quantum circuits and the results from those circuits
def state_tomo(result, st_qcs):
    # The expected final state; necessary to determine state tomography fidelity
    target_state = (One^One^Zero).to_matrix()  # DO NOT MODIFY (|q_5,q_3,q_1> = |110>)
    # Fit state tomography results
    tomo_fitter = StateTomographyFitter(result, st_qcs)
    rho_fit = tomo_fitter.fit(method='lstsq')
    # Compute fidelity
    fid = state_fidelity(rho_fit, target_state)
    return fid

def tomography(result, st_qcs, target_state):
    # Fit state tomography results
    tomo_fitter = StateTomographyFitter(result, st_qcs)
    rho_fit = tomo_fitter.fit(method='lstsq')
    # Compute fidelity
    fid = state_fidelity(rho_fit, target_state.to_matrix())
    return fid

def average_fidelity(jobs, st_qcs, result_state,  verbose = True):
    # Compute tomography fidelities for each repetition
    fids = []
    for job in jobs:
        
        # uncomment here for submission
        #fid = state_tomo(job.result(), st_qcs)
        fid = tomography(job.result(), st_qcs,  result_state)
        fids.append(fid)

    if verbose:
        print('state tomography fidelity = {:.4f} \u00B1 {:.4f}'.format(np.mean(fids), np.std(fids)))
    
    return np.mean(fids)


def gen_testing_circuit(trotter_function, trotter_steps, num_qubits = 7, target_time = np.pi,  terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'],  backend = sim, operational_qubits = [1, 3, 5], compiler = False, pauli_compiler = False, coupling_map = None, shuttle = False, qsd = False):
    q1, q2, q3 = operational_qubits
    
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
        
    # Prepare initial state (remember we are only evolving 3 of the 7 qubits on jakarta qubits (q_5, q_3, q_1) corresponding to the state |110>)
    qc.x([operational_qubits[1], operational_qubits[2]])  # DO NOT MODIFY (|q_5,q_3,q_1> = |110>)
    
    # generate circuit, compile, parse to instructions
    qr2, time_ev_circ = basic_trotter_circuit(trotter_function, trotter_steps, terms = terms, target_time = target_time, num_qubits = len(operational_qubits), compiler = compiler )
    
    
    if pauli_compiler:
        time_ev_circ = pauli_transpilation(time_ev_circ, 3, backend = backend, coupling_map = coupling_map)
        
    if qsd:
        qsd_circ = gen_qsd(Operator(time_ev_circ).data)
        time_ev_circ = transpile(qsd_circ, optimization_level = 3, backend = backend)
        
    circuit_instructions = time_ev_circ.to_instruction()
    
    if shuttle:
        q1, q2, q3 = 3, 5, 6
        qc.x([qr[3], qr[6]])
        qc.append(circuit_instructions, [qr[q1], qr[q2], qr[q3]])
        qc.append(swap_gate(), [qr[1], qr[3]])
        qc.append(swap_gate(), [qr[3], qr[5]])
        qc.append(swap_gate(), [qr[5], qr[6]])
    else:
        qc.append(circuit_instructions, [qr[q1], qr[q2], qr[q3]])
        
    print("Circuit Description ", qc.decompose().count_ops())
    
    return qr, qc


def test_circuit(trotter_function, trotter_steps, target_time = np.pi, num_qubits = 7, repetitions = 5, backend_name = None, backend = sim, terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], noise_model = None, compiler = False, pauli_compiler =  False, coupling_map = None, active_qubits = 3, shuttle = False, qsd = False):
    """
    Runs trotterized circuit and returns average fidelity
    Inputs: trotter_function, function which returns QC instructions for trotter gate
            trotter_steps, int, number of trotter repetitions
            repetitions, int, number of experiments averaged over
            backend
            terms, list of strings, indicating order of decomposition
    """
    
    # need to fix thiss.....
    if num_qubits <= 3 or backend_name == 'quito':
        q1, q2, q3 = 0, 1, 2
    elif (not backend_name) or backend_name == 'jakarta':
        q1, q2, q3 = 1, 3, 5
        
    operational_qubits = [q1, q2, q3]
        
    # generate circuit, compile, parse to instructions
    qr, qc = gen_testing_circuit(trotter_function, trotter_steps, target_time = target_time, num_qubits = num_qubits, terms = terms, operational_qubits = operational_qubits, compiler = compiler,  coupling_map = coupling_map, shuttle=shuttle, pauli_compiler = pauli_compiler, qsd = qsd )
    
    # tomography circuit
    st_qcs  =  state_tomography_circuits(qc, [qr[q1], qr[q2], qr[q3]])
        
    # assign jobs and evaluate
    jobs = evaluate_circuit(st_qcs, reps = repetitions, backend = backend)
    
    result_state = prepare_state(target_time)
    fid = average_fidelity(jobs, st_qcs, result_state)
    
    return fid

def test(instructions, target_time = np.pi, num_qubits = 7, backend = sim, shuttle = False, verbose = False):
    repetitions = 3
    
    # need to fix thiss.....
    if num_qubits <= 3:
        q1, q2, q3 = 0, 1, 2
    elif num_qubits >= 6:
        q1, q2, q3 = 1, 3, 5
    
    qr = QuantumRegister(num_qubits)
    qc = QuantumCircuit(qr)
    
    if shuttle:
        qc.x([qr[5], qr[6]])
        qc.append(instructions, [qr[3], qr[5], qr[6]])
        qc.append(swap_gate(), [qr[1], qr[3]])
        qc.append(swap_gate(), [qr[3], qr[5]])
        qc.append(swap_gate(), [qr[5], qr[6]])
    else:
        qc.x([qr[q2], qr[q3]])
        qc.append(instructions, [qr[q1], qr[q2], qr[q3]])
    
    """if verbose:
        print(qc.decompose().decompose().count_ops())
        print(qc.decompose().draw())"""
    
    qc2 =  QuantumCircuit(qr)
    qc2.append(transpile(qc, backend, optimization_level = 0), qr)
        
    if verbose:
        print("Circuit Description")
        print(qc2.decompose().count_ops())
        print(qc2.decompose().draw())
        
    # tomography circuit
    st_qcs  =  state_tomography_circuits(qc2, [qr[q1], qr[q2], qr[q3]])
        
    # assign jobs and evaluate
    jobs = evaluate_circuit(st_qcs, reps = repetitions, backend = backend, verbose = verbose)
    
    result_state = prepare_state(target_time)
    fid = average_fidelity(jobs, st_qcs, result_state)
    return fid
   
###### Testing Sistematically ###########

def steps_experiment(trotter_function, n_max , n_min = 4, target_time = np.pi, step=1, repetitions = 3, backend_name=None, backend = sim, terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], noise_model = None, compiler = False, shuttle=False, pauli_compiler = False, qsd = False):
    """
    Repeatedly runs test_circuit to evaluate a given trotter decomposition as we increase the number of steps
    Inputs: trotter_function, function which returns QC instructions for trotter gate
            n_max, int, max number of trotter_steps
            n_min, int, min number of trotter_steps
            step, int, spacing between trotter_steps
            repetitions, int, number of experiments averaged over
            backend
            terms, list of strings, indicating order of decomposition
    
    """
    
    fids = []
    n_steps = []
    
    for i in range(n_min, n_max, step):
        n_steps.append(i)
        fids.append(test_circuit(trotter_function, i, target_time = target_time, repetitions = repetitions, backend_name = backend_name, backend = backend, noise_model = noise_model, compiler = compiler, pauli_compiler=pauli_compiler, shuttle=shuttle, qsd = qsd))
        
    return n_steps, fids        

def time_experiment(trotter_function, trotter_steps , time_steps = 10, max_time = np.pi, step=1, repetitions = 3, backend_name=None, backend = sim, terms = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2'], noise_model = None, compiler = False, shuttle=False, pauli_compiler = False):
    
    ts = [ i * max_time/time_steps for i in range(1, time_steps + 1)]
    fids = []
    
    for t in ts:
        print("Simulating until time ", t)
        fids.append(test_circuit(trotter_function, trotter_steps, target_time = t, repetitions = repetitions, backend_name = backend_name, backend = backend, noise_model = noise_model, compiler = compiler, pauli_compiler=pauli_compiler, shuttle=shuttle))
        
    return ts, fids

def plot_results(x, y, plot_name = None, x_axis_name = None, y_axis_name = "Average State Fidelity"):
    
    plt.scatter(x, y)
    
    if plot_name:
        plt.title(plot_name)
        
    if x_axis_name:
        plt.xlabel(x_axis_name)
    if y_axis_name:
        plt.ylabel(y_axis_name)
        
    plt.show()



