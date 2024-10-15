## UQ for the PB equation in molecular electrostatics

This software "shakes" pqr files according to thermal fluctuations, to generate many conformations that run with the [PBJ](https://github.com/bem4solvation/pbj/tree/main/pbj).

To run:

```
python sampler.py -nt N_TESTS -nw N_WORKERS -pqr BASE_PQR_FILE -f FOLDER
python MC_Solver.py -f FOLDER_WITH_PQRS
python get_moments.py -f FOLDER
```
* `sampler.py` generates the "shaken" pqrs, and stores them in `FOLDER/job_X/file.pqr`. Inside `FOLDER` there are `N_WORKERS` `job_X` folders.
* `MC_Solver.py` runs `PBJ` on pqr structures that are in `FOLDER_WITH_PQRS`, (*ie.*, `FOLDER_WITH_PQRS=FOLDER/job_0/`) and stores the results in a `.csv` file.
* `get_moments.py` finds all the `.csv` files in `FOLDER/job_X` and computes mean and variace. This is stored in an `output_summary.csv` file in `FOLDER`.

There are a few more flags available in each case. To explore them, use the `--help` flag. 
