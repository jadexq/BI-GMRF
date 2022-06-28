## Inference

Apply the BI-GMRF algorithm (Section 3) with the folder `BI_GMRF`.
Taking 2D simulation -- butterfly scenario as an example.

##### Setup:

- Prerequisites:
    - Eigen 3.3.7
    Details at: https://eigen.tuxfamily.org/index.php?title=Main_Page
    - cmake 3.3.1
    - C++ 14
- Set eigen path in `CmakeLists.txt` by revising `include_directories(<your-path>/eigen3 .)`
- Compile and run:
``` shell
mkdir build
cd build
cmake ..
make
./BI-GMRF
```

##### Inputs:

- neighbourhood information provided at `BI-GMRF/data/`: 

  `nbh.txt` -- the $v$th row contains the index of the neighbourhoods ($\partial_v$) of the $v$th pixel

  `nNbh.txt` -- the $v$th row is the number of neighbourhoods ($|\partial_v|$) of the $v$th pixel

- outputs of the data step 1 (genderate data)
   `x_gen.txt` -- exposure $x_i$
   `m_gen.txt` -- imaging mediator $m_i(v)$
   `y_gen.txt` -- outcome $y_i$
   `c1_gen.txt` -- confounder $c_{1i}$

  (put these under `BI-GMRF/data/` together with `nbh.txt` and `nNbh.txt`)

##### C++ code documentation:

- `main.cpp` the main function

- `Constant` (C++ class) for specify model settings and constants

  see comments in `Constant.h` for details
  
  symbol `<<` marked out variables can be specified e.g. sample size
  
- `Initial` (C++ class) set initial values

- `Mcmc` (C++ class) gibbs sampling

- `Generate` (C++ class, not in use) generate data

- functions

  `rnorm` -- for generating standard normal random numbers

  `runif` -- for generating uniform(0,1) random numbers

  `runifVec` -- for generating uniform(0,1) random vectors

  `strToVec` -- convert a string to vector (double)

- `CMakeLists.txt` the cmake list I use, just for your reference

##### Outputs:

Compile the source code with cmake and run the executive file gives .txt files save in `BI-GMRF/output/` including:

- `rep_alpha` , `rep_zetaA` -- $\alpha$ and $\zeta^{\alpha}$ estimates (each row for one replication)
- `rep_beta` , `rep_zetaB` -- $\beta$ and $\zeta^{\beta}$ estimates
- `rep_lambda1` -- $\lambda_1$ estimates
- `rep_gamma` -- $\gamma$ estimates
- `rep_sigE` -- $\sigma^2_{\epsilon}$ estimates

One can costomize the outputs by modifying `main.cpp`

##### Other scenarios:

For acquiring the estimates corresponding to other 2D scenarios, following the steps:
1. Replace the simulated data (`x_gen.txt`, `m_gen.txt`, `y_gen.txt`, `c1_gen.txt`) at `BI-GMRF/data/`, keep the neighbourhood information unchanged (`nbh.txt`, `nNbh.txt`). 

2. Update `Constant.h`. 

   For example, for `square` scenario, instead of using the original `Constant.h` (the one for `butterfly` scenario), please copy and paste the one at `BI_GRMF/const/square/` to `BI-GRMF/` (replace the original one). 
   The `Constant.h` corresponding to `square`, `cross`, `triangle` and `circle` scenarios are provided in `BI_GRMF/const/` under the folder with the same name.
   
   