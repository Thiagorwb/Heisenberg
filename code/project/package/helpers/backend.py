from .libraries import *
"""
Here we place all access codes to backend
"""

token = '2559ad320e7ab14dd9ccff0a9b0bdbc9218b8ca0a62e063cfe2fd3de05b9fbb194861a202c3d7c213b456e4cc8446cf2bd19560114b48be79cb8dc751b2e1a49'

IBMQ.save_account(token, overwrite=True)

def fetch_backend(name = 'sim', simulated = True, noisy = False, verbose = False):
    
    if (not name) or name == 'sim':
        if verbose:
            print("Fetching Simulator")
        return QasmSimulator()
    
    provider = IBMQ.load_account()
    provider = IBMQ.get_provider(hub='ibm-q-community', group='ibmquantumawards', project='open-science-22')
    backend = provider.get_backend(name)
    
    if (not simulated):
        if verbose:
            print("Fetching the Quantum Computer")
        return backend
    if simulated or noisy:
        if verbose:
            print("Fetching the Simulator for ", name)
        return QasmSimulator.from_backend(backend)

def get_basis(backend):
    return backend.configuration().basis_gates

def get_map(backend):
    return backend.configuration().coupling_map