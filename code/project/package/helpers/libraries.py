# plotting and numpy
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})  # enlarge matplotlib fonts

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

# qiskit stuff

# Import qubit states Zero (|0>) and One (|1>), and Pauli operators (X, Y, Z)
from qiskit.opflow import Zero, One, I, X, Y, Z, CX
from qiskit.quantum_info.operators import Operator, Pauli

# Importing standard Qiskit modules
from qiskit import QuantumCircuit, QuantumRegister, IBMQ, execute, transpile
from qiskit.converters import dag_to_circuit, circuit_to_dag

from qiskit.providers.aer import QasmSimulator
from qiskit.tools.monitor import job_monitor
from qiskit.circuit import Parameter
from qiskit.providers.aer.noise import NoiseModel
from qiskit.quantum_info import OneQubitEulerDecomposer, TwoQubitBasisDecomposer, Statevector

# Import state tomography modules
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter
from qiskit.quantum_info import state_fidelity


# optimizers
import scipy
from scipy import linalg
from scipy.linalg import cossin
from scipy.stats import unitary_group
import openql as ql
import quantum_decomp