from .qsd import *

##################
## Basic Graphs ##
##################

def clique(n):
    """ 
    Computes adjacency matrix for a Fully Connected Graph
    Input: n, int, size
    Output: List of 2-tuples, the edges
    """
    a = []
    for i in range(n-1):
        for j in range(i+1, n):
            a.append([i,j])
    return a

def line(n):
    """ 
    Computes adjacency matrix for a Line Graph
    Input: n, int, size
    Output: List of 2-tuples, the edges
    """
    a = []
    for i in range(n-1):
        a.append([i,i+1])
    return a

##################
## Graph Parser ##
##################

class Graph:
    """
    Graph Class - underlying map that the circuit lives on
    
    Inputs 
    - nvertices, int, number of vertices, 
    - coupling_map, list of 2-tuples, list of undirected edges

    Methods 

    - construct-adjacency, constructs adjacency list of graph
    - APSP, computes all pairs shortest paths in G
    - Query_path, returns path between query qubits i and j
    
    """
    
    def __init__(self, nvertices,  coupling_map = None):
        """
        Inputs 
        - nvertices, int, number of vertices, 
        - coupling_map, list of 2-tuples, list of undirected edges

        """

        # coupling_map must be list of edges
        if not coupling_map:
            print("Where is your graph?")
            
        # initialize
        self.n = nvertices
        self.map = coupling_map
        
        # precompute APSP
        self.adjacency = self.Construct_adjacency()
        self.distances, self.preds = self.APSP()
        
    def Construct_adjacency(self):
        """
        Constructs adjacency matrix
        """
        
        adjacency = {}
        
        for edge in self.map:
            u, v = edge
            
            adjacency[u] = adjacency.get(u, [])
            adjacency[u].append(v)
            
            adjacency[v] = adjacency.get(v, [])
            adjacency[v].append(u)
                    
        return adjacency
        
    def APSP(self):
        """
        Computes all pairs shortest paths using Floyd-Warshall's algorithm
        """
        dist = {}
        pred = {}
        
        # initialize distances and paths for all vertices
        for u in range(self.n):
            dist[u] = {}
            pred[u] = {}
            for v in range(self.n):
                dist[u][v] = self.n
                pred[u][v] = None
            
            # in place
            dist[u][u] = 0
            pred[u][u] = None

            for v in self.adjacency[u]:
                dist[u][v] = 1
                pred[u][v] = u

        for mid in range(self.n):
            for u in range(self.n):
                for v in range(self.n):
                    pathlength = dist[u][mid] + dist[mid][v]
                    if pathlength < dist[u][v]:
                        dist[u][v] = pathlength
                        pred[u][v] = pred[mid][v]

        return (dist, pred)
    
    
    def Query_path(self, i, j):
        """ queries path in couping map from i to j"""
        
        path = [j]
        dist = self.distances[i][j]
        
        current = j
        for d in range(dist):
            prev = self.preds[i][current]
            path.append(prev)
            current = prev
        path.reverse()
        
        return path



class Graph_Parser:
    
    # Question --- can one compress the circuit too? Occasional repeated cnots...
    
    """
    
    Graph Parser class
    
    -  initializes underlying circuit topology,
    - parses fully connected circuit to a circuit of given topology

    Inputs 
    - nvertices, int, number of vertices, 
    - coupling_map, list of 2-tuples, list of undirected edges
    
    Methods
    - parse: given an input quantum circuit, returns another quantum circuit defined on the topology    
    
    """
    
    def __init__(self, nvertices,  coupling_map = None):
        """
        Inputs 
        - nvertices, int, number of vertices, 
        - coupling_map, list of 2-tuples, list of undirected edges
        """
        
        if not coupling_map:
            print("Where is your graph?")
            
        # initialize
        self.n = nvertices
        self.map = coupling_map
        
        # precompute APSP
        self.graph = Graph(self.n,  self.map)
        
    def parse(self, circuit):
        """
        Inputs 
        - circuit, Quantum Circuitt 

        Returns
        - quantum circuit instructions, with 2-qubit gates only supported on adjacenct qubits in the coupling_map

        Note: Input quantum circuit must only have CNOTs as 2 qubit gates - consider compiling or transpiling first
        """
        
        # parse circuit into DAG
        dag = circuit_to_dag(circuit)
        
        # initialize new qc
        qr = QuantumRegister(self.n)
        new_qc = QuantumCircuit(qr)

        # parse each operation at a time, expanding two qubit cxs into long range gates
        for node in dag.topological_op_nodes():
            
            inst = node.op.copy()
        
            if inst.name != 'cx':
                q1 = node.qargs[0]
                p1 = q1.index
                new_qc.append(inst, [qr[p1]])

            if inst.name == 'cx':
                
                # index controle and target qubits
                q1, q2 = node.qargs[0], node.qargs[1]
                p1, p2 = q1.index, q2.index
                
                # query path in map
                path = self.graph.Query_path(p1, p2)
                
                # perform long range cx
                longcx = cx_on_a_path(path, self.n)
                new_qc.append(longcx, qr)
                
        return new_qc.decompose().to_instruction()
        

def geometric_compile(circuit, nqubits, coupling_map):
    """
    Wrapper function for geometric (graph) compilation scheme

    Inputs:
    - circuit, QuantumCircuit, circuit to be compiled into the geometry
    - nqubits, int, number of qubits on the underlying quantum chip
    - coupling_map, list of 2-tuples, list of undirected edges

    Function initializes the Graph_Parser class on the coupling map, and compiles the circuit

    Returns:
    - Quantum Circuit instructions, with 2-qubit interactions all adjacent on the underlying coupling map

    Note: Input quantum circuit must only have CNOTs as 2 qubit gates - consider compiling or transpiling first

    """
    
    parser = Graph_Parser(nqubits, coupling_map)
    instr = parser.parse(circuit)
    
    return instr





