# Cardiorespiratory control computational model: Simulation and fitting

This folder is composed of a series of files which try to execute the simulation of a cardio respiratory control model, also its fitting with data. So there are mixed files, ones dedicated just for the simulation, others just for fitting, others for data processing, and others for analysis of results.

### Simulation

For the simulation the following files are needed:

- `ForwardModel`: it's like the `main` for the simulation, it contains:
  - `run_ode`: it's the solver's driver:
    - `model_basic`: it's the model, is the code which contains all the equations that are used within the solver
*IMPORTANT:
To run ForwardModel, is very important the hyperparameters assigned in the set_up functions. To run the simulation with the fitted parameters the dates of the fitted parameters must be assigned.
Also, if respiratory times will be given, it has to be specified. The default for the moment is this:
[setup] = set_up("simulation", patient_idx, state, "mix", "dt", 0.1,'pars_from_fitting', 1, 'fitting_mat_file', {'07-09-2025', '30-09-2025', '20-10-2025', '24-10-2025', '02-11-2025'}, 'time_from_data', 1);
"pars_from_fitting" allows loading parameters from fitted files
"fitting_mat_file" are all the fitting dates from which the parameters will be loaded
"time_from_data" allows to take respiratory times from data.
  
 
### Sensitivity Analysis

To run sensitivity analysis use the function  `parameter_analysis_fun` which recieves: 
patient_number: char (1,2,3,4);
pindex: char; the parameter index of the parameter to be perturbed
load_rb: char (1,0); is 1 if results_base will be loaded and 0 if it will be computed and saved as a new results_base.
 This function calls:
     -  `sens_functions`: the driver for all the senstivities computations: run base, compute perturbed simulations and sensitivities, build STensor. It calls:
          -  `data_processing`
          -  `run_ode`

### Fitting

`ForwardModelFitting`: it runs the fitting operation, contains:

* `data_preprocessing`: Gives the data preprocessed for fitting
* `obj_fun`: it is the objective function which computes the error of simulation vs experimental data, runs:

  * `run_ode`


### Data processing

`data_preprocessing`: recieves the original data and tries to read it and process it, to generate a decent dataset to fit. It also computes the polynomial estimations for VO2, VCO2 and fiO2:

* `bestPolynomialFit`: it computes the polynomial estimation. It allows to plot the data and the polynomial function


### Others

#### Data related

    `estimate_newton_ohm` :Reads parameters file, and using ohm newton equation based on simple parameters of each subject, estimate a bunch of parameters of vasculature.

    `readBeatscopeData20`: Its an special file to read data from FINAPRES results.

#### Parameter pre-estimation

    `sensitivity_analysis`: Its a code similar to `ForwardModel` with the difference that runs multiple instances with parallel computing to perfome the sensitivity analysis of the model

    `sensitivity_lab`: a tool to test the variation of parameters in the model

    `identificability_analysis`: as it says, computes the FIM from the sensitivity matrix, to reveal the amount of correlation between parameters.

#### Setting up

    `vectorize_dicts`: it transforms both solver's driver and the model from dictionary based structures to vector based structures, using:

`transformation_imported_dict`: using regex, it goes all into `run_ode` and `model_basic` to know which lines to update into vector type

`Optimize_percentages:` its a code to know which percentages of different vascular conditions have to be defined in order to estimate correctly the parameters using newton ohm equations. It's not needed anymore, after is already used.

#### Plots

    `plot_for_document:` Different codes for plotting figures that will go into the document

    `deploy_papers_results`: For plotting the steady state simulation
