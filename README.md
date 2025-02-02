# Fast Maximal Independent Sets on Dynamic Graphs

This is the computational artefact for the paper **Fast Maximal Independent Sets on Dynamic Graphs**. It includes the source
code for *experiments* in three separate `c++` files: `tvbl.cpp`, `tvb.cpp`, and `inc_Bat.cpp`.

#### tvbl.cpp
1. **TVBL Class**: 
   - Represents the Tuple-Valued Bitmask with Levelling graph structure
   - Main data structure: `tvbl` (Stores TVBL-specific information for each node)
2. **PDMI_TVBL Function**: Implements Parallel Dynamic MIS Insertion for TVBL
3. **PDMD_TVBL Function**: Implements Parallel Dynamic MIS Deletion for TVBL
4. **readCSR Function**: Reads CSR format and MIS from files, initializes the TVBL object
5. **Main Function**: Orchestrates the experiment workflow

#### tvb.cpp
- Similar structure to `tvbl.cpp`
- Implements the TVB (Tuple-Valued Bitmask) data structure
- Contains PDMI_TVB and PDMD_TVB functions for insertion and deletion

#### inc_Bat.cpp
- Uses a graph with CSR data structure
- Contains PDMI_INC function for incremental updates
- PDMI_BAT is implemented directly in the main function
- No separate deletion function

All three implementations use OpenMP for parallel processing.

## Dependencies and requirements

We run all experiments a server that has two AMD EPYC
7742 64-Core Processors running at 2.25 GHz, built on the
x86 64 architecture. Each processor core encompasses a 32-
bit and 64-bit CPU op-mode and operates with 2 threads per
core. The system boasts 128 CPU cores distributed across
2 sockets, enabling extensive parallel processing capabilities.
The processor ensures efficient data access and management
with 512 MB of L3 cache spread across 32 instances. Ac-
companied by 512 GB of DDR4 memory, the server runs on
Ubuntu 22.04.

We use `13` graphs in *Matrix Market (.mtx)* file format from the *SuiteSparse Matrix Collection* and converted them in the CSR format (.txt) as our input dataset. These must be placed in the `~/csr`
directory **before running** the experiments. Thereafter, we create the batches of different sizes, which must be placed in the `~/batches` directory, and scratch MIS is calculated for the graph using Luby's algorithm, which is stored as `.txt` files in the `~/mis` directory.

## Compilation

To compile the source files, use the following commands:

```bash
g++ -O3 -fopenmp tvbl.cpp -o tvbl
g++ -O3 -fopenmp tvb.cpp -o tvb
g++ -O3 -fopenmp inc_Bat.cpp -o inc_Bat
```

These commands compile the source files with optimization level 3 (-O3) and OpenMP support (-fopenmp) for parallel processing.

## Running the Experiments

To run an individual algorithm, use the following command template:

```bash
./<executable> <CSR_file> <MIS_file> <Deletion_folder> <Insertion_folder> <num_threads> <Num_batch_for_deletion> <Num_batch_for_insertion>
```

Example:
```bash
./tvbl csr/kron.txt mis/kron.txt ./batches/kron/3 ./batches/kron/3 64 10 10
```

This command runs the PDMI-TVBL and PDMD-TVBL algorithms on the `Kronecker` graph, using 64 threads, with 10 batches each for deletion and insertion operations.

## Citation  
If you find our work useful, please consider citing:  

**Fast Maximal Independent Sets on Dynamic Graphs**, Prajjwal Nijhara, Aditya Trivedi, Dip Sankar Banerjee, *Proceedings of Euromicro PDP 2025*


```bibtex
@inproceedings{Nijhara2025FastMIS,
  author    = {Prajjwal Nijhara and Aditya Trivedi and Dip Sankar Banerjee},
  title     = {Fast Maximal Independent Sets on Dynamic Graphs},
  booktitle = {Proceedings of Euromicro PDP},
  year      = {2025}
}
