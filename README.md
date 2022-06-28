# BI-GMRF

This repo provides data and code used in the paper 

`A High-dimensional Mediation Model for a Neuroimaging Mediator: Integrating Clinical, Neuroimaging, and Neurocognitive Data to Mitigate Late Effects in Pediatric Cancer`

With the data and code provided, the 2D simulation -- butterfly scenario in Section 4.1 can be reproduced.

The other simulation scenarios can be reproduced by slightly modifying the code.



### Step 1.  Generate data
Generate simulation data with the folder `gene` by running `gene.m`

##### Setup:
* Install  [MATLAB Tensor Toolbox, Version 3.1](http://www.tensortoolbox.org/)
* Install the MATLAB toolbox `TensorReg` in `gene/`

##### Script: 
`gene.m` provided at `gene/`.
The *inputs* are `.png` files specifying shapes of true imaging coefficients; the *outputs* are `.txt` files including generated data and parameter true values. 
Details are provided in folder `gene`.



### Step 2. Inference
Apply the BI-GMRF algorithm (Section 3.1) with the folder `BI_GMRF` (C++ project).

##### Setup:

- Prerequisites:
    - [Eigen 3.3.7](https://eigen.tuxfamily.org/index.php?title=Main_Page)      
    - cmake 3.3.1
    - C++ 14
- Set eigen path in `CmakeLists.txt` by revising `include_directories(<your-path>/eigen3 .)`
- Compile and run in the terminal:
``` shell
mkdir build
cd build
cmake ..
make
./BI-GMRF
```
The BI-GRMF binary program should be stored in a folder, which is in the same folder with the data folder storing the input files and the output files storing the results. 

##### Code:
Main function: `BI_GMRF/main.cpp`.
The *inputs* are the neighbourhood information and (simulated) data; the *outputs* are the MCMC samples and estimated parameters.
Details are provided in folder `BI_GMRF`



### Step 3. Evaluate and visualize the estimates

One are encouraged to use any preferred software (e.g. `R`) to evaluate and visualize the estimates provided by `BI-GMRF`. Some sample codes are provided in `evaluation/`.

