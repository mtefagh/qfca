
Quantitative Flux Coupling Analysis
====

*Quantitative flux coupling analysis* (QFCA) is a quantitative approach to identify and remove the redundancies of the steady-state flux cone. The full details of our approach are discussed in the associated papers.

Example
-------
The following code uses QFCA to compute the table of flux coupling relations and the list of blocked reactions for the E. coli core model and also returns the reduced metabolic network.
```Matlab
load('ecoli_core_model.mat');
[S_reduced, rev_reduced, fctable, blocked] = QFCA(model.S, model.rev, 'linprog', 1e-6);
```
The output of the above code is as follows.
```Matlab
Identifying the blocked reactions and removing them from the network: 1.109
Reduced number of:
	metabolites = 72;	reactions = 95;	nonzero elements = 360
Finding the full coupling relations: 0.031
Reduced number of:
	metabolites = 39;	reactions = 60;	nonzero elements = 279
Correcting the reversibility types: 0.063
Finding the directional and partial couplings: 0.797
Metabolic network reductions postprocessing: 0.000
Inferring by the transitivity of full couplings: 0.141
Reduced number of:
	metabolites = 39;	reactions = 58;	nonzero elements = 276
```

License
----

QFCA is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).