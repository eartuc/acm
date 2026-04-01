# Simple solution code for Artuc, Chaudhuri and McLaren (2010)

I realized that the replication package provided with the paper was unclear and inefficient. When we wrote the code around 2006, Matlab did not support broadcasting, Julia and Python 3 did not even exist. I rewrote the simulation code in Python, Julia and Matlab using a much simpler and shorter solution algorithm that I have been using for many years (since 2010, I think). Each script is self-contained, with references to equations in the paper. In this version, I am considering the simplest case in our paper, where all sectors are traded and prices are exogenous.

See Artuc and Ortega (2026) for a more comprehensive dicussion on solving static trade models, but with many more moving parts.  

## Solution Steps 

### Step 0: 

First, guess $V^i_t$ and $L^i_t$. Note that $L^i_0$ is the initial labor allocation before the shock. 

### Step 1: 

Calculate consumer price index: 

$\phi_t = \prod \left( p^i_t \right)^{\theta^i}$

then calculate real wages using MPL: 

$w^i_t = \frac{p^i_t}{\phi_t} \alpha^i \psi^i \left( (L^i_t)^{\rho^i-1} \right) \left(\alpha^i (L^i_t)^{\rho^i-1} + (1-\alpha^i) \right)^{\frac{1}{\rho^i}-1}$

### Step 2: 

Calculate moving probabilities: 

$m^{ij}_t = \frac{\exp \left( \left[ \beta V^j_t - C^{ij}\right]/\nu \right)}{\sum_k \exp \left( \left[ \beta V^k_t - C^{ik}\right]/\nu \right)}$

### Step 3: 

Calculate option values:

$\Omega^i_t = -\nu \log \left( m^{ii}_t \right)$

### Step 4: 

Implied values:

$V_t^{i*} = w^i_t + \beta V^i_{t+1} + \Omega^i_t$

### Step 5: 

Implied labor allocation:

$L_{t+1}^{i*} = \sum_k L^k_t m^{ki}_t$

### Step 6: 

Check errors:

$\xi_1  = \sum_i \sum_t |V^{i*}_t - V^{i}_t |$

$\xi_2  = \sum_i \sum_t |L^{i*}_t - L^{i}_t |$

Update:

$V^{i}_t = \vartheta V^{i*}_t + (1-\vartheta) V^{i}_t$

$L^{i}_t = \vartheta L^{i*}_t + (1-\vartheta) L^{i}_t$

If errors are small then stop, otherwise go to Step 1.






