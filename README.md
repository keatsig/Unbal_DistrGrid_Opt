# ThreePhasePowerModels
Solve optimal power flow for unbalanced three-phase distribution systems 

## Nonlinear three-phase OPF formulations:
1. `TPOPF_pol()`: Conventional OPF in polar coordinates (bus injection) 
2. `TPOPF_rect()`: Conventional OPF in rectangular coordinates (bus injection) 
3. `TPOPF_ivr()`: Conventional OPF in polar coordinates (branch flow) 

## Linearized three-phase OPF formulations:
1. `FOT_opf_pol()`: Iterative first-order Taylor (FOT) approximation OPF in polar coordinates (bus injection)
2. `FP_opf_pol()`: Iterative fixed-point (FP) OPF in polar coordinates (bus injection)
3. `FP_opf_rect()`: Iterative FP-OPF in rectangular coordinates (bus injection)
4. `BFS_opf_rect()`: Iterative backward-forward sweep (BFS) OPF in rectangular coordinates (branch flow)
4. `D3F_opf_rect()`: LinDist3Flow OPF (branch flow)

## Three-phase power flow formulations:
1. `FP_pf()`: Fixed-point method
1. `NR_pf()`: Newton-Raphson method
1. `BFS_pf()`: Backward-forward sweep method

## Running code:
To run the code, use the following command:
```@example overview
include("main.jl"); TPOPF("case_name",obj_op)
``` 
   - *case_name*: subfolder name in data folder to choose test feeder
   - *obj_op*: choose objective function 
     - 1: no objective
     - 2: minimize VUF
	 - 3: minimize substation power injection
	 - 4: minimize PVUR
	 - 5: minimize LVUR
	 - 6: minimize network losses with unbalance constraints
	 - 7: minimize solar PV reactive power injections