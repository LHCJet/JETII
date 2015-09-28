JETJet II
=================================================================
An implementation of the JET II algorithm for finding 2-prong jets, described in arXiv:1509.07522.

Requirements:
--------------------------------------------
```
GCC 
FastJet
```

Quick Start:
--------------------------------------------
The main code is distributed as a single header file, ```double_cone.hh```. It requries ```fastjet``` to compile and run. An example is given in ```example.cc```, which reads a sample datafile with two events and find the leading jets. To run the example, download all files and put them in the same directory. Make sure ```fastjet-config``` is callable or you may need to edit the ```Makefile``` file. Then do
```
make
./example
```

If the code executes correctly, the output should be the same as ```output.txt```. 
Usage:
--------------------------------------------
See the description at the beginning of ```double_cone.hh```.
