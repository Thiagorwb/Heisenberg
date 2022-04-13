from .qsd import *

class Graph:
    
    """
    Graph Class - underlying map that the circuit lives on
    
    Inputs - int, number of vertices, and list coupling_map, list of directed edges
    Methods - construct-adjacency, 
            - apsp
            - query path, returns path between query qubits i to j
    
    """
    
    def __init__(self, nvertices,  coupling_map = None):
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
        
        adjacency = {}
        
        for edge in self.map:
            u, v = edge
            
            adjacency[u] = adjacency.get(u, [])
            adjacency[u].append(v)
            
            adjacency[v] = adjacency.get(v, [])
            adjacency[v].append(u)
                    
        return adjacency
        
    def APSP(self):
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
    
    Methods
    
    - parse 
    Inputs: Qiskit Quantum circuit, returns qiskit quantum circuit defined on topology
    
    
    """
    
    def __init__(self, nvertices,  coupling_map = None):
        # coupling_map must be list of edges
        
        if not coupling_map:
            print("Where is your graph?")
            
        # initialize
        self.n = nvertices
        self.map = coupling_map
        
        # precompute APSP
        self.graph = Graph(self.n,  self.map)
        
    def parse(self, circuit):
        
        # parse circuit into DAG
        dag = circuit_to_dag(circuit)
        num_qubits = dag.num_qubits()
        
        # initialize new qc
        qr = QuantumRegister(num_qubits)
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
                longcx = cx_on_a_path(path, num_qubits)
                new_qc.append(longcx, qr)
                
        return new_qc.decompose().to_instruction()
        

def geometric_compile(circuit, nqubits, coupling_map):
    
    parser = Graph_Parser(nqubits, coupling_map)
    instr = parser.parse(circuit)
    
    return instr