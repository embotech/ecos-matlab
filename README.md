ecos-matlab
===========

Matlab interface for ECOS, the Embedded Conic Solver. For detailed documentation please visit [this page](https://www.embotech.com/ECOS).

Repository Content
====

You will find the following directories in this repository:

* `bin`: directory with the script for making a MEX binary (makemex.m) and other helper files for ecos-matlab.
* `conelp`: Matlab implementation of ECOS with different linear system solver options.
* `src`: mex interface C file
* `test`: testing code - run batchtest.m to run different tests

Using ECOS in MATLAB
====

ECOS can be used in three different ways within MATLAB:

- native MEX interface: [learn more](https://www.embotech.com/ECOS/Matlab-Interface/Matlab-Native)
- via CVX [learn more](https://www.embotech.com/ECOS/Matlab-Interface/CVX)
- via YALMIP [read more](https://www.embotech.com/ECOS/Matlab-Interface/Yalmip)

In either case, you need a compiled mex binary for your platform, which you can download from [embotech](https://embotech.com/ECOS/Download). Alternatively, you can also build it from source (the bottom part of the [embotech download page](https://embotech.com/ECOS/Download) explains how).
