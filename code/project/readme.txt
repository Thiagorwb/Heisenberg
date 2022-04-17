

Hi!

here you can find my submission to the IBMQ open science prize. 

In main.ipynb, you can find my main submission. Please run it for me! I only managed the time to run it on noisy simulations of Jakarta.

In this folder you should be able to find my writeup, where I go over the main methods used. Here one can find their implementation:

In package, you will find a python package for the tools I implemented. 
- in package/helpers, i import libraries and define other simple linear algebra tasks
- in package/circuits, one can find the implementation of 
	-- higher order product formulas and other circuits, in addition to basic gates and hamiltonian helper functions
- in package/compiler, one can find the implementation of
	-- the QSD compiler, see package/compiler/qsd.py
	-- the gray code optimization, used in qsd.py
	-- the geometric circuit parser, see 
- in package/testing, 
	-- tomography and plotting helper functions
	-- basic testing scripts to run the circuits, see package/testing/testing.py
	-- main experimental scripts, see package/testing/experiments.py

Finally, in this folder you will also find a series of other ipynbs, which I placed here for reproducibility of my plots. Namely,

- In Compiler Experiments, the experiments on compiling Haar random unitaries, and the gate count of the time evolution circcuits

- In 3 Qubit Experiments, the main experiments on 3 Qubit systems on the IBM Jakarta, Manila and Quito quantum chips

- In Brickwork Experiments, the proof of concept experiments for the brickwork architecture for larger number of qubits
