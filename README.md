EUTER Example Program
=======

This is a test setup, which shows a minimal example how to use the BrainScaleS C++ interface. 
The example follows the sample Python-based script found in the [guidebook](https://electronicvisions.github.io/hbp-sp9-guidebook/pm/using_pm_newflow.html) to a big extend.


Compiling
----------
In contrast to the BS software stack, we make use of cmake for compiling the code. 
After pulling this repository, please execute the following inside the project main folder:

```bash
mkdir build
cd build 
run_nmpm_software cmake ..
run_nmpm_software make 
```
 
In this case, we assume that the repository is situated somewhere at the BrainScaleS front-end server, and that the respective modules are loaded, e.g. execute 
```bash
module load nmpm_software/current
```

Compiling will create an executable called euter_test. Execute via 
```bash
srun -p experiment --wmod 33 run_nmpm_software ./euter_test
```
