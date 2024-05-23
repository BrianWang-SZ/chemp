*chemp* is a C++ program developed by Ruohe (Brian) Wang based on the *Programming Projects* by the Crawford Group. This program requires *Psi4* as a python module to generate the molecular integrals and utilizes the Eigen library, a C++ template library, as its linear algebra solver.

# File Preparation (`psi4chemp.py`)
You can find a python script named *psi4chemp.py* in the `src` directory to generate the one-electron (kinetic and electron-attraction) and two-electron integrals and the nuclear repulsion energy from *Psi4*.

## Environment Setup
*Psi4* is used as a python module for this script. You can create a conda environment for *Psi4* using the following command:
```
conda create -n p4env python=3.8 psi4 pydantic<2.0 -c psi4/label/dev
```

## Usage
After activating the environment, you can use the following command to generate the required files:
```
python psi4chemp.py <geometry_file> <method> <basis_set> [directory]
```
The argument `geometry_file` is the path to the input geometry, which supports the same input file types as *Psi4* (e.g., `xyz`, `zmat`, ...). The argument `method` and `basis_set` specifies the level of theory to be used. The `method` needs to be specified to obtain the hessian matrix for frequency analysis. The argument `directory` is an optional path you can specify to put the generated files. The default directory would be the current directory. **Note that the directory needs to be created before running the command.**

The *Psi4* output file `output.dat` will also be saved in the specified directory. This will contain results of energy and frequency analysis at the specified level of theory. 

You should see "Files prepared! Ready to chemp." if the files are generated successfully.

# Calculation (`chemp`)
## Compilation
The files required for the compilation of `chemp` are located in the `src` directory. With a `gcc` compiler, you should be able to run the following command to compile the executable file in the `src` directory based on the provided `Makefile`:
```
make
```
After compilation, you should see the executable named `chemp` in the `src` directory. You can run the following command to undo the compilation:
```
make clean
```
## Usage
After compiling `chemp`, you can use the executable in the following manner in the `src` directory:
```
./chemp [-i path] [-m vib|geom|hf|mp2|ccsd] [-s] [-n cycle] [-e decimal] [-d decimal] [-h] [-v]
```
See below for the explanations for the options:


-i

The directory of the prepared files (default: current directory)


-m

The operation to perform, including vibrational analysis (vib), geometry analysis (geom), restricted Hartree-Fock (hf), second-order Møller-Plesset perturbation (mp2) and coupled-cluster singles and doubles (ccsd). **Note that CCSD has error.** (default: "hf")


-s

Whether to use the Direct Inversion of the Iterative Subspace (DIIS) technique for the Hartree-Fock calculation (default: false)


-n  

Max number of iteractions for SCF methods (default: 100)


-e

The decimal place for the tolerance of energy convergence (default: 12)


-d

The decimal place for the tolerance of density matrix convergence (default: 11)


-v

Turn off detail printing for the Hartree-Fock method (default: on)


-h   

Print usage information

# Test Files
In the `geom` directory, there are some files to be used for `psi4chemp.py` for H_2, CH_4, and H_2O. In the `int` directory, the molecular integrals for all three molecules with the `STO-3G`, `DZ`, and `DZP` basis set. 


# Example
Here is an example of a sample calculation with diis calculation in the `src` directory with DIIS with a maximum of 1000 cycles and a convergence of 6 decimal places for both energy and density matrix.
```
> mkdir ../test
> source activate p4env
> python psi4chemp.py ../geom/h2o.zmat hf sto-3g ../test
> ./chemp -i ../test/ -m hf -s -n 1000 -e 6 -d 6
```
