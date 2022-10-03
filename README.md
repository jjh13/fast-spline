
# Fast-Spline
This is the companion code repository to the paper "Generating Fast Interpolants for Lattice Samplings". This code 
performs both the analysis of spline spaces, as well as the LLVM code 
generation. 

## Who is this repository for?
If you have an awkward multivariate piece-wise polynomial spline, and want to generate fast interpolation code,
then this repository is for you.

## Getting Started
- SageMath 9.2: The original code was written for SageMath 6.8, however SageMath
  is under rapid development, and breaks compatability often. Moreover, older versions of
  Sage are difficult to find, and no longer hosted by SageMath. The code has been updated
  to work with Sage 9.2, but hasn't been as thourougly tested as the code for 6.8. If you
  notice any bugs, please open an issue. Sage 9.2 can be installed by following the
  instructions found here: https://doc.sagemath.org/html/en/installation/
  
- llvmlite: This is needed for code generation, specifically it's a lightweight
  LLVM Python binding. Keep in mind you'll also need LLVM.

- LLVM: This is only neccesary if you want to use the code you've generated. I use the LLVM toolchain that comes with Ubuntu, installed simply by using
  ```sudo apt-get install llvm``` this may change based on your distribution.


## Running the Code
You can choose to get started in two ways. Using a Jupyter notebook is more user friendly, but can become somewhat 
cumbersome when you're trying to automate tasks. To solve that issue, we've also provided a simple command line utility.

#### Jupyter
If you want to experiment with Jupyter, first start a Jupyter notebook server within Sage via 
```sage --notebook=jupyter```. Make sure to extract this repository to a location that Sage will be able 
to find. Then, open the ```Open Me First.ipynb```. This
notebook will install llvmlite, and unpack the Voronoi splines. Next, to get a feel for
the code, you'll want to check out the ```2D Examples```, ```3D Examples``` and 
```4D Examples``` notebooks --- I suggest looking at them in that order. 

#### Command Line
To automate testing, we've also provided a simple command line utility ``fastspline.sage``. 
Before running this command, first run ``sage -m pip install llvmlite``. This will ensure that
Sage has access to the llvmlite library. Next, run the command
``sage fastspline.sage -h`` to see details on how to generate code for splines.

##### Example


## Reproducing Results
To reproduce the results in the paper, first run the command ``make codegen -j[n]``. This will generate all the code for 
all the test cases in the paper.

## Roadmap
While this codebase is relatively general, it has a number of issues that still need to be addressed before 
- Modularize: The code needs to be broken down into smaller modules. Currently there are a handful of monolithic files that perform the analyis and code generation. I'd like to break these down into smaller, testable components.
- Remove dependencies: SageMath is a huge software suite with many dependencies (that changes quite often). Eventually I want to remove this dependency (and rewrite everything in a more performant language).
- BB Form: We need to eventually provide and evaluate generating code for the BB-form. This requires a bit of theoretical work.
- Remove Cruft: There's a lot of left over code from earlier stages of prototyping, which needs to be removed.

## License and Citation 
This project is licensed under the MIT licence (see LICENSE). In that sense, it is completely free (as in beer). However if you use this in your paper/project please cite

```
@InProceedings{Horacsek2022,
title = {FastSpline: Automatic Generation of Interpolants for Lattice Samplings},
author = {Joshua Horacsek, Usman Alim},
}
```
