# Zebrafish-stripe-model
Model for zebrafish pattern formation - as in elife paper (accepted): A quantitative model for zebrafish pattern formation by Jennifer P. Owen, Robert N Kelsh and Christian Yates.

This folder contains the following files suitable for use with matlab:

The user is reccommended to use the 'reader' file for simulations, followed by the plotting functions for visualisation.

Reader function
1. reader.m - run this file to generate plots for zebrafish mutants; WT, nac, shd, pfe, rse, cho, sbr, seurat, leo and save for different time points (subject to answers from prompts)
Full list of mutant types:
WT, nac, pfe, shd, leo, shd_pfe, nac_pfe, nac_shd, leo_shd, leo_nac, leo_pfe, sbr, seurat, ablation, vertical_stripe, ablate_iridophore, move_up_stripe, initially_stripey, tall_domain, small_domain

Main function
1. main.m - This file is the main document which generates the zebrafish simulations. For best use of this file, run reader.m

Inputs 
mutant_type (WT, nac, shd etc), 
time_str- 1 if want to save the domain at every dev stage, 2 if every half day or 3 if only at the end point of the simulation (juvenile stage)
name_dir - name of the directory that you want the domain at different time points to be saved in (if using in reader this is generated automatically)

Outputs:
domain_matrix,domain_matrix_ir,domain_matrix_X
Matrices corresponding to the M, I and X layers respectively.
For the matrices listed above if a space (i,j) corresponds to a 0, then it is empty.
m=1; %melanocytes
xb=2; %xanthoblasts
xd=4; %xanthophores
il=5; %loose iridophore
id=6; %dense iridophore

Subfunctions (in order of appearance)

1. determine_mutant- for mutant that affect cell number updates the appropriate parameters relating to the initial conditions and end points
Input
mutant_type (WT, nac, shd, pfe, rse, cho, sbr, seurat, leo etc)
Outputs
final_rec - number of important timepoints for the mutant
N_xb - initial number of xanthoblasts
N_m - initial number of melanophores
size_rec - list of important timepoints for the mutant
xan_on - whether or not xanthophores are present in mutant
mel_on - whether or not melanophores are present in mutant
irr_on - whether or not iridophores are present in mutant

2. fixed_parameters - creates a list of parameters that are fixed across time (not influenced by number of cells).
Input - mutant_type
Outputs - all parameters that are fixed across time.

3. Initial_conditions - generates the initial condition of the given simulation.
Inputs
mutant_type (WT, nac, shd, pfe, rse, cho, sbr, seurat, leo etc)
N_xb - initial number of xanthoblasts (as determined in determine_mutant)
N_m - initial number of melanophores (as determined in determine_mutant)
mel_on - whether or not melanophores are present in mutant
irr_on - whether or not iridophores are present in mutant
Outputs - The initial conditions for the simulation

4. check_timed_events - updates model based on events that occur during timeline.
Inputs - domain_matrix at time  t (with other relevant parameters)
Outputs - updated domain_matrix at time  t given timed events occur (with other relevant parameters)

5. parameters - generates parameters that change with time.
Input: numbers of cells etc.
Outputs - parameter list for time t (for Gillespie algorithm)

6. FIND - finds position on domain where a cell of a certain type is located.

7. calc_AR_p_speedup_diag_m - Function checks the surroundings of R,C when cell type is located on X or I and outputs matrices of the numbers of cells in each position - to be used with the weighting vector (based on attraction/repulsion to cells to determine probability of moving in each given direction

8. calc_AR_p_speedup_diag_m - Function checks the surroundings of R,C when cell type is located on M and outputs matrices of the numbers of cells in each position - to be used with the weighting vector (based on attraction/repulsion to cells to determine probability of moving in each given direction

9. FIND_any - finds position on domain where a cell of any type is located.

10. prolif_p_speedup - updates the model based on a proliferation event

11. calc_neighbours - Uses the uniform norm to calculate number of neighbours of type x on the SAME domain matrix at distance max_distance from space R,C.

12. calc_neighbours_m -Uses the uniform norm to calculate number of neighbours to an m of type x on the x, i domain matrix at distance max_distance from space R,C.

13. calc_neighbours_xi_count_m -Uses the uniform norm to calculate number of neighbours to an x or an i of type m at distance max_distance from space R,C.

14. check_empty - checks a position is empty of all cell types in all layers

15. check_iridophore_loose_or_dense - determines whether an iridophore should be loose or dense based on its surroundings.

16. p_movement_diag - updates the model based on a movement event

17. pull_xb_p - updates the model based on a 'pull' event

Saving files
1. PROF_SAVE_NOW - saves the domain matrices by calling ...
2. parsave_2 - save the domain_matrices

Plotting files
1. plot_and_save_zebrafish - plots and saves the cell positions determined from the domain matrices as a visual pattern
2. plot_zebrafish - plots the cell positions determined from the domain matrices as a visual pattern without saving



