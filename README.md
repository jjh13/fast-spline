# Readme
----
This is the companion code to the papers "Generating Fast Interpolants for Lattice Samplings" 
Part I & II. This code performs both the analysis of spline spaces, as well as the LLVM code 
generation. This code has been tested under Linux and may have issues running under other 
platforms. 



## Dependencies
- SageMath 9.2 --- The original code was written for SageMath 6.8, however SageMath
  is under rapid development, and breaks compatability often. Moreover, older versions of
  Sage are difficult to find, and no longer hosted by SageMath. The code has been updated
  to work with Sage 9.2, but hasn't been as thourougly tested as the code for 6.8. If you
  notice any bugs, please open an issue. Sage 9.2 can be installed by following the
  instructions found here: https://doc.sagemath.org/html/en/installation/
  
- llvmlite --- This is needed for code generation, specifically it's a lightweight
  LLVM Python binding. Keep in mind you'll also need LLVM.

- LLVM --- I use the LLVM toolchain that comes with Ubuntu, installed simply by using
  ```sudo apt-get install llvm``` this may change based on your distribution.


## Getting Started
First, start a Jupyter notebook server within Sage via ```sage --notebook=jupyter``` and 
make sure this repository is accessible. Then, open the ```Open Me First.ipynb```. This
notebook will install llvmlite, and unpack the Voronoi splines. Next, to get a feel for
the code, you'll want to check out the ```2D Examples```, ```3D Examples``` and 
```4D Examples``` notebooks --- I suggest looking at them in that order. To get a better
feel for the code that can be generated, take a look at the ```Codegen Details``` notebook.
