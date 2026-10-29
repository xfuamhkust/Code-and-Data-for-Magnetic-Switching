# Magnetic Switching

This repository contains Python codes used for the FM and AFM switching simulations in our paper *A Unified Theory of Deterministic Magnetic Switching*.

There are two main files:
- `FM_LLG.py` for ferromagnets (FM)
- `AFM_LLG.py` for antiferromagnets (AFM)

## How to run

The codes use Python and NumPy (`pip install numpy`).

To run a simulation, open `FM_LLG.py` or `AFM_LLG.py`, change the parameters at the bottom of the file using the tables below, and run:

```bash
python FM_LLG.py
# or
python AFM_LLG.py
```

Each run saves text files in the working folder. Please move or rename them before running another case, since the next run may overwrite them.

## Parameters

The following parameters are used in all eight cases:

| Parameter | Value |
|---|---|
| `alpha` | `0.1` |
| `m` | `1.72` (Bohr magnetons) |
| `K` | `1e-23` J |
| `k` | `[0,0,1]` |
| `dt` | `1e-15` s |

The vectors for spin directions and torque polarizations are normalized in the code. `ttF` and `ttD` give the start time, duration, and time after the pulse, in time steps.

### FM

| Case | `tF` | `tD` | `pF` | `pD` | `m0` | `ttF = ttD` |
|---|---:|---:|---|---|---|---|
| FLT, up to down | -1 | 0 | `[0,.02,1]` | `[0,.02,1]` | `[0,0,1]` | `[0,8000,0]` |
| FLT, down to up | 1 | 0 | `[0,.02,1]` | `[0,.02,1]` | `[0,0,-1]` | `[0,8000,0]` |
| DLT, up to down | 0 | -0.2 | `[0,.02,1]` | `[.01,-.01,1]` | `[0,0,1]` | `[0,4000,0]` |
| DLT, down to up | 0 | 0.2 | `[0,.02,1]` | `[.01,-.01,1]` | `[0,0,-1]` | `[0,4000,0]` |

Here, FLT and DLT mean field-like and damping-like torque. The code converts `tF` and `tD` by multiplying them by `1e-21/(m*muB)`, with `muB=9.27400999e-24` J/T.

### AFM

For all AFM cases, use `J=-1e-21` J, `d=[1e-22,0,0]` J, `pF=[0,0,1]`, `tF2=tF1`, and `tD2=tD1`.

| Case | `pD` | `tF1 = tF2` | `tD1 = tD2` | `hN` | `ttF` | `ttD` |
|---|---|---:|---:|---|---|---|
| FFS, up to down | `[0,1,0]` | 0 | 0.01 | `0.5*tD1` | `[0,0,6000]` | `[0,3000,3000]` |
| FFS, down to up | `[0,1,0]` | 0 | -0.01 | `0.5*tD1` | `[0,0,6000]` | `[0,3000,3000]` |
| FAS, up to down | `[1,0,0]` | 0.02 | -0.01 | 0 | `[0,13000,0]` | `[0,10000,3000]` |
| FAS, down to up | `[1,0,0]` | 0.02 | 0.01 | 0 | `[0,13000,0]` | `[0,10000,3000]` |

The initial spin directions are:

| Case | `m10` | `m20` |
|---|---|---|
| FFS, up to down | `[0,.0493241086,.998782825]` | `[0,.0493241086,-.998782825]` |
| FFS, down to up | `[0,-.0493241086,-.998782825]` | `[0,-.0493241086,.998782825]` |
| FAS, up to down | `[0,-.000248833713,.999999969]` | `[0,.0997611693,-.995011412]` |
| FAS, down to up | `[0,.000248833713,-.999999969]` | `[0,-.0997611693,.995011412]` |

FFS means field-free switching, and FAS means field-assisted switching. The code converts the torque parameters `tF1`, `tF2`, `tD1`, and `tD2` by multiplying them by `abs(J)/(m*muB)`. For FFS, set `hN=0.5*tD1` **after** this conversion.

## Results

The codes save `time_data.txt` and `magnetization_data.txt`. The AFM code also saves `neelvector_data.txt`. The time file contains step numbers; multiply them by `0.001` to get time in ps.
