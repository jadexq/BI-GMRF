## Generate data

Generate simulation data with the folder `gene` by running `gene.m`.
Taking 2D simulation -- butterfly scenario as an example.

##### Setup:

* Install  [MATLAB Tensor Toolbox, Version 3.1](http://www.tensortoolbox.org/)
* Install TensorReg in `gene/`

##### Inputs:

- `alphaTrue.png` -- shape of $\alpha(v)$

- `betaTrue.png`  -- shape of $\beta(v)$

provided at `gene/input/`

##### Script: 

`gene.m` provided at `gene/`

##### Outputs:

Run `gene.m` to generate outputs as .txt files saved in `gene/output/` including

- data: 
   `x_gen.txt` -- exposure $x_i$
   
   `m_gen.txt` -- imaging mediator $m_i(v)$
   
   `y_gen.txt` -- outcome $y_i$
   
   `c1_gen.txt` -- confounder $c_{1i}$

- true values: 
   `alpha_true.txt` (matrix form), `alphaVec_true.txt` (vectorized form) -- $\alpha(v)$ true value
   
   `beta_true.txt`, `betaVec_true.txt` -- $\beta(v)$ true value
   
   `lambda1_true.txt`, `lambda1Vec_true.txt` -- $\lambda_1$ true value
   
   `gamma_true.txt` -- $\gamma$ true value

##### Other scenarios:

For generating simulation data corresponding to other 2D scenarios, simply replace the `gene/input/alphaTrue.png` and `gene/input/betaTrue.png`.

The true imaging-related coefficients shapes (`.png` files) used for `square`, `cross`, `triangle` and `circle`  are provided in folders with the same name under `gene/input/`.
