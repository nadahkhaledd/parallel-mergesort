### Execution options:
            1. Sequentially: run only original sequential merge sort. by sending run argument (e.g. s).
            2. Parallel: run only added parallel implementation for merge sort. by sending run argument (e.g. p).
            3. run sequential and parallel together and calculate speedup for comparison. by sending run argument "both".


#### Commands: 
        > make
        > ./mergesort-co <number-elements> <mode: s, p, both> <number-threads>


#### Output examples: 
##### 1
   > ./mergesort-co 10000 both 4

        Initialization...
        Sorting 10000 elements of type int (0.038147 MiB)...
        Sequential execution completed in 0.002502 sec.
        Start parallel implementation with  4 threads.
        Parallel execution completed in 0.010697 sec.
        Speedup: 0.233897
        Verification... successful.

##### 2
   > ./mergesort-co 200000 p 5

        Initialization...
        Sorting 2000000 elements of type int (7.629395 MiB)...
        Start parallel implementation with  5 threads.
        Parallel execution completed in 0.198699 sec.
        Verification... successful.

##### 3
   > ./mergesort-co 2000 s 

        Initialization...
        Sorting 2000 elements of type int (0.007629 MiB)...
        Sequential execution completed in 0.000400 sec.
        Verification... successful.