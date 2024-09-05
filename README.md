# Code for "Dynamics of plant phenome can be accurately predicted from genetic markers".

This code was developed and run on operating systems: `Ubuntu 22.04 LTS and 24.04 LTS`

This code was developed and run using: `R version 4.3.3`

R packages used: 
- `Matrix (1.6-5)`
- `pracma (2.4.4)`
- `reshape2 (1.4.4)`
- `ggplot2 (3.5.1)`

Code was developed and run on a 12th Gen Intel® Core™ i7-1265U.

Package can be installed using: 
`devtools::install_github("dobby978/dynamicGP")`

Install time is within 1 to 10 minutes on a "normal" desktop computer.

Load data:
`data <- read.csv("example_data.csv")[, -1]`

A list containing a tensor of A matrices directly calculated with algorithm 1 along with a tensor of X matrices is returned by: 
`DMD_components <- find_DMD_components(data, method = "calculated.A")`

A corresponding list containing tensors of A_r along with all intermediate component matrices for algorithm 2 is returned by: 
`DMD_components <- find_DMD_components(data, r = 2, method = "schur.DMD")`

Longitudinal prediction accuracies across the complete time series can be determined with:
`accuracies <- accuracy_aggregate_and_plot_training(DMD_components[[1]], DMD_components[[2]])`

Tensors can be converted into line by matrix entry matrices with `tensor_unroller()`.

- `U_r_transpose_all <- tensor_unroller(DMD_components$U_r_transpose_all)`
- `V_r_all           <- tensor_unroller(DMD_components$V_r_all)`
- `sigma_r_all       <- tensor_unroller(DMD_components$sigma_r_all)`
- `A_tilde_all       <- tensor_unroller(DMD_components$A_tilde_all)`
- `Q_all             <- tensor_unroller(DMD_components$Q_all)`
- `R_all             <- tensor_unroller(DMD_components$R_all)`
- `phi_all           <- tensor_unroller(DMD_components$phi_all)`

DMD_components can be calculated and reshaped into matrices within a few seconds. 
