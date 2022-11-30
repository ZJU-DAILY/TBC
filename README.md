# Exact and Approximate Temporal Betweenness Computation

The source code for "Efficent Exact and Approximate Betweenness Centrality Computation for Temporal graphs"

## Running Enviroment

a 64-bit Linux-based OS;

## Exact Temporal Betweenness Computation

### How to run

step1: Download the whole data and code from "data-set" and "exact tbc computation" directory.

step2: Put the data in the corresponding code folder.

step3: Enter "exact tbc computation" directory  and Select parameters in source code.

step4: Execute command line statement:

```shell
g++ -O3 -fopenmp -o run main.cpp
./run
```

## Approximate Temporal Betweenness Computation

### Files

All of the C++ source code can be found in "exact tbc computation" directory

### How to run

step1: To compile and run ATBC, you need the source of exact tbc computation and ABRA, the NetworKit sources and the NLopt library.

step2: To obtain the Network sourcesplease see the [Networkit](https://github.com/networkit/networkit.git)

step3: To install NLopt, please see the [NLopt homepage](http://ab-initio.mit.edu/wiki/index.php/NLopt)

step4: To obtain the Network sources, please download the sources from [ABRA sources](http://matteo.rionda.to/software/ABRA-radebetw.tbz2) and then copy the `.cpp` and `.h` files in file `Abra/src/abra` into the directory `NetworKit/networkit/cpp/centrality/` and `NetworKit/include/networkit/centrality/` , where `Abra` is the directory where you have the ABRA sources and   `NetworKit` is the directory where you have the NetworKit sources

step5: Then, create a file `atbc` in `NetworKit` and copy sources in both `TBC/approximate tbc computation` and `TBC/exact tbc computation` into it

step6: Select parameters in source code

step7: Compile and run

```shell
g++ sample.cpp ../networkit/cpp/centrality/RadeAux.cpp -fopenmp -L /usr/local/lib/ -lnetworkit -lnlopt -O3 -o sp
./sp
```

