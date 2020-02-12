# state-space-diffusion-approximation

Code to obtain results from Johnston, Simpson and Crampin "Predicting population extinction in lattice-based birth-death-movement models"

To run code, run "Comparison_Script.m" in MATLAB. All relevant parameters are defined at the top of "Comparison_Script.m"

"Comparison_Script.m" will perform the random walk ("Random_Walk.m"), solve the state space diffusion approximation PDE ("SSDA_PDE_Solver.m") and solve the mean field ODE ("Mean_Field_ODE.m"), and finally compare results obtained from all three methods ("Plot_Results.m")

The remaining files are functions called  during the aforementioned script to: to map coordinates to a vector index ("Mapping.m"); map a vector index  to coordinates ("Reverse_Mapping.m"), and; solve a tridiagonal matrix system ("Thomas_Algorithm.m")
