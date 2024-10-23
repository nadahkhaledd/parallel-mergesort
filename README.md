### Execution options:
            1. Sequentially: run only original sequential merge sort. by sending run argument (e.g. s).
            2. Parallel: run only added parallel implementation for merge sort. by sending run argument (e.g. p).
            3. run sequential and parallel together and calculate speedup for comparison. by sending run argument "both".


#### Commands: 
        > make
        > ./mergesort-co <number> <mode: s, p, both>


#### Output examples: 
##### 1
   > ./mergesort-co 1000 both

        Initialization...
        Sorting 1000 elements of type int (0.003815 MiB)...
        Sequential execution completed in 0.000154 sec.
        Parallel execution completed in 0.000425 sec.
        Speedup: 0.362353
        Verification... successful.

##### 2
   > ./mergesort-co 200000 par

        Initialization...
        Sorting 200000 elements of type int (0.762939 MiB)...
        Parallel execution completed in 0.049664 sec.
        Verification... successful.

##### 3
   > ./mergesort-co 2000 s 

        Initialization...
        Sorting 2000 elements of type int (0.007629 MiB)...
        Sequential execution completed in 0.000400 sec.
        Verification... successful.