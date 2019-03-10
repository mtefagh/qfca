Quantitative Flux Coupling Analysis
====

Introduction
----
*Quantitative flux coupling analysis* (QFCA) is a quantitative approach to identify and remove the redundancies of the steady-state flux cone. The full description and theory of the algorithm implemented by QFCA is discussed in the associated papers. When referencing QFCA, please cite the following:
- Tefagh, M. & Boyd, S.P. J. Math. Biol. (2018). [https://doi.org/10.1007/s00285-018-1316-9](https://doi.org/10.1007/s00285-018-1316-9)
- Tefagh, M. & Boyd, S.P. bioRxiv (2018). [https://doi.org/10.1101/499251](https://doi.org/10.1101/499251)

![Metabolic Network Reductions](https://connect.biorxiv.org/qr/qr_img.php?id=499251)

Quick Start
----
The following example uses the `QFCA` function to compute the table of **flux coupling relations** and the list of **blocked reactions** for the [core *E. coli* model](http://systemsbiology.ucsd.edu/Downloads/EcoliCore) and also returns the **reduced metabolic network**.
```matlab
load('ecoli_core_model.mat');
[reduced_net, fctable, blocked] = QFCA(model, true, 'linprog');
```
The output of the above code is as follows.
```matlab
Original number of:
	metabolites = 72;	reactions = 95;	nonzero elements = 360
Original number of:
	reversible reactions = 59;	irreversible reactions = 36
Identifying the blocked reactions and removing them from the network: 0.034
Reduced number of:
	metabolites = 72;	reactions = 95;	nonzero elements = 360
Finding the full coupling relations: 0.008
Reduced number of:
	metabolites = 34;	reactions = 60;	nonzero elements = 236
Correcting the reversibility types: 0.037
Finding the directional and partial coupling relations: 0.411
Inferring by the transitivity of full coupling relations: 0.000
Metabolic network reductions postprocessing: 0.001
Reduced number of:
	metabolites = 30;	reactions = 58;	nonzero elements = 224
The number of solved:
	linear programs = 54;	systems of linear equations = 8
```
Furthermore, the `directionallyCoupled` function can be utilized as a stand-alone function to find **fictitious metabolite** certificates.

License
----
*QFCA* is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html).
