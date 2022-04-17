from ..helpers import *
from ..circuits import *
from ..compiler import *


###### Evaluation and Tomography ###########

def evaluate_circuit( st_qcs, backend = sim, shots = 5000, reps = 8, verbose = False):
    """
    Repeated evaluation and tomography of quantum circuit
    Inputs: st_qcs, state tomography circuit
            backend, IBMQ backend to run evaluation on
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
# THIS IS IBM CODE

def state_tomo(result, st_qcs):
    # The expected final state; necessary to determine state tomography fidelity
    target_state = (One^One^Zero).to_matrix()  # DO NOT MODIFY (|q_5,q_3,q_1> = |110>)
    # Fit state tomography results
    tomo_fitter = StateTomographyFitter(result, st_qcs)
    rho_fit = tomo_fitter.fit(method='lstsq')
    # Compute fidelity
    fid = state_fidelity(rho_fit, target_state)
    return fid


# Compute the state tomography based on the st_qcs quantum circuits and the results from those circuits

def tomography(result, st_qcs, target_state):
    """
    Generic tomography procedure which compares output state to noiseless simulated result, and computes fidelity. 
    This function lets us test our circuits for all target times, as opposed to just np.pi

    Inputs:
    result, job.result output
    st_qcs, state_tomography_circuits 
    target_state, numpy array of size 2**3, the ideal output quantum state

    Returns: fidelity, float
    """

    # Fit state tomography results
    tomo_fitter = StateTomographyFitter(result, st_qcs)
    rho_fit = tomo_fitter.fit(method='lstsq')
    # Compute fidelity
    fid = state_fidelity(rho_fit, target_state)
    return fid

def average_fidelity(jobs, st_qcs, target_state,  verbose = True):
    """
    Computes average fidelity of jobs
    Inputs:
    jobs, 
    st_qcs, state_tomography_circuits
    target_state, numpy array of size 2**3, the ideal output quantum state 
    verbose, boolean, if True prints all individual fidelities

    Returns: average fidelity, float
    """

    # Compute tomography fidelities for each repetition
    fids = []
    for job in jobs:
        
        # uncomment here for submission
        #fid = state_tomo(job.result(), st_qcs)
        fid = tomography(job.result(), st_qcs,  target_state)
        fids.append(fid)

    if verbose:
        print('state tomography fidelity = {:.4f} \u00B1 {:.4f}'.format(np.mean(fids), np.std(fids)))
    
    return np.mean(fids)

def output_prob(result, st_qcs, vector):
    """
    Generic tomography procedure which computes the squared inner product of the
     circuit output and a fixed computational basis state 

    Inputs:
    result, job.result output
    st_qcs, state_tomography_circuits 
    vector, numpy array, a fixed computational basis state 

    Returns: fidelity, float
    """

    # Fit state tomography results
    tomo_fitter = StateTomographyFitter(result, st_qcs)
    rho_fit = tomo_fitter.fit(method='lstsq')
    # Compute inner product
    prob = np.abs(vector @ rho_fit @ vector)
    return prob


def average_prob(jobs, st_qcs, vector = psi_0(3),  verbose = True):
    """
    Computes average fidelity of jobs
    Inputs:
    jobs, 
    st_qcs, state_tomography_circuits
    vector, numpy array, a fixed computational basis state 
    verbose, boolean, if True prints all individual fidelities

    Returns: average fidelity, float
    """

    # Compute tomography fidelities for each repetition
    probs = []
    for job in jobs:
        
        # uncomment here for submission
        #fid = state_tomo(job.result(), st_qcs)
        p = output_prob(job.result(), st_qcs,  vector)
        probs.append(p)

    if verbose:
        print('probability of measuring psi = {:.4f} \u00B1 {:.4f}'.format(np.mean(p), np.std(p)))
    
    return np.mean(probs)


