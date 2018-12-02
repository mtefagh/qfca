Quantitative Flux Coupling Analysis
====

Introduction
----
*Quantitative flux coupling analysis* (QFCA) is a quantitative approach to identify and remove the redundancies of the steady-state flux cone. The full description and theory of the algorithm implemented by QFCA is discussed in the associated papers. When referencing QFCA, please cite the following:
- 1st paper
- 2nd paper

Quick Start
----
The following example uses the `QFCA` function to compute the table of **flux coupling relations** and the list of **blocked reactions** for the [core *E. coli* model](http://systemsbiology.ucsd.edu/Downloads/EcoliCore) and also returns the **reduced metabolic network**.
```matlab
load('ecoli_core_model.mat');
[S_reduced, rev_reduced, fctable, blocked] = QFCA(model.S, model.rev, true, 'linprog');
```
The output of the above code is as follows.
```matlab
Original number of:
	metabolites = 72;	reactions = 95;	nonzero elements = 360
Original number of:
	reversible reactions = 59;	irreversible reactions = 36
Identifying the blocked reactions and removing them from the network: 1.531
Reduced number of:
	metabolites = 72;	reactions = 95;	nonzero elements = 360
Finding the full coupling relations: 0.016
Reduced number of:
	metabolites = 39;	reactions = 60;	nonzero elements = 279
Correcting the reversibility types: 0.047
Finding the directional and partial coupling relations: 1.094
Inferring by the transitivity of full coupling relations: 0.000
Metabolic network reductions postprocessing: 0.000
Reduced number of:
	metabolites = 37;	reactions = 58;	nonzero elements = 274
The number of solved:
	linear programs = 54;	systems of linear equations = 6
```

License
----
QFCA is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).