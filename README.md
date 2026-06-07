# Epidemics Engine

A high-performance epidemic simulation engine written in Fortran 90, designed to model disease spread on complex networks using SIS and SIR models.

**Author:** Adrià Rojo  
**Year:** 2026

---

## 📋 Overview

This project implements a computationally efficient simulator for epidemic models on network topologies. It supports:

- **Multiple epidemic models:** SIS (Susceptible-Infected-Susceptible) and SIR (Susceptible-Infected-Recovered)
- **Network types:** Weighted and unweighted networks
- **High-performance computing:** OpenMP parallelization for multi-threaded execution
- **Flexible input:** Support for various network formats
- **Detailed output:** Statistics and event tracking for analysis

The core simulation engine is written in Fortran 90 for maximum performance, with analysis and visualization tools provided via Jupyter notebooks.

---

## 🔧 Requirements

### System Requirements
- **Operating System:** Linux (tested on Linux with aarch64 architecture)
- **Compiler:** `gfortran` (GNU Fortran compiler)

### Build & Development
- **IDE:** Visual Studio Code (recommended)
- **Conda:** Anaconda/Miniconda (for dependency management)

---

## 🚀 Quick Start

### 1. Compile the code

Choose one of the following compilation scripts based on your use case:

```bash
# Standard single-threaded build
./build.sh

# Parallel build with OpenMP support
./build-parallel.sh

# Alternative main program (main2.f90)
./build-main2.sh
```

### 2. Run a simulation

After compilation, the executable will be ready to run. The simulation parameters are configured in the main program files.

### 3. Clean build artifacts

```bash
./clear-builds.sh
```

---

## 📁 Project Structure

```
epidemics-tfg/
├── main.f90                    # Primary simulation program
├── main2.f90                   # Alternative simulation program
├── skiplist_test.f90           # Skip list data structure tests
├── include/                    # Fortran module headers
├── post-process/               # Post-processing scripts
├── files/                      # Network input files
├── output/                     # Simulation output directory
├── build.sh                    # Build script (standard)
├── build-main2.sh              # Build script (main2 variant)
├── build-parallel.sh           # Build script (OpenMP)
├── clear-builds.sh             # Clean build artifacts
├── compile_modules.sh          # Module compilation helper
└── fortls.config.json          # Fortran LS configuration

```

---

## 📊 Key Features

### Simulation Models

#### SIS Model
The Susceptible-Infected-Susceptible model where infected individuals recover and can be reinfected.

#### SIR Model
The Susceptible-Infected-Recovered model where recovered individuals gain immunity.

### Network Handling

- **Network Loading:** Reads edge lists from standard network format files
- **Weighted Networks:** Full support for weighted network topologies
- **Sparse Representation:** Efficient memory usage for large networks
- **Hash Map Management:** Custom hash map implementation (`fhash`) for non-consecutive node indices

### Parallelization

The code includes OpenMP directives for:
- Parallel simulation batches
- Thread-safe file I/O
- Dynamic scheduling for load balancing

---

## 📦 Dependencies

### Fortran Dependencies

#### fhash
- **Description:** Dictionary/hash map implementation for Fortran 90
- **Purpose:** Handles networks with non-consecutive node numbering
- **Source:** [jl2922/fhash](https://github.com/jl2922/fhash/blob/master/fhash.f90)
- **License:** MIT

#### mt19937
- **Description:** Mersenne Twister pseudo-random number generator
- **Purpose:** Provides deterministic random number generation with seed control
- **Type:** Built-in module

---

## 💻 Usage

### Configuration

Edit `main.f90` or `main2.f90` to configure:

1. **Network file path**
   ```fortran
   open(unit=11, file='./files/test-file-1.txt', action='read')
   ```

2. **Simulation parameters**
   ```fortran
   call execute_simulation(net, &
       infection_rate, recovery_rate, seed, &
       time_limit, model_type, &
       stats_unit=21, events_unit=22, net_name=name)
   ```

3. **Output locations**
   ```fortran
   open(unit=stats_unit, file='./output/stats-'// filename //'.dat', action='write')
   ```

### Running Simulations

```fortran
! Single simulation
call execute_simulation(net, 2.0_dp, 1.0_dp, 42072, 50.0_dp, SIR_MODEL, &
    stats_unit=21, events_unit=22, net_name='mynetwork')

! Batch parallel execution
!$omp parallel do private(i_lambdas) schedule(dynamic)
do i_lambdas = 100, 1, -1
    call execute_simulation(net, real(i_lambdas, dp)/100, 1.0_dp, &
        42069+i_lambdas, 100.0_dp, SIR_MODEL, stats_unit=(i_lambdas+50))
end do
!$omp end parallel do
```

### Output Files

- **Stats file:** `./output/stats-<name>.dat` - Time series of infected/recovered densities
- **Events file:** `./output/events-<name>.dat` - Timestamped event log (infection/recovery events)
- **Metadata:** Header comments in output files with simulation parameters

---

## 🔨 Advanced Build Options

### Compile Specific Modules

```bash
./compile_modules.sh
```

### Manual Compilation

```bash
gfortran -fopenmp -o epidemics_sim *.f90 -Iinclude
```

### Release Build

Add optimization flags:

```bash
gfortran -O3 -fopenmp -o epidemics_sim *.f90
```

---

## 📝 Configuration Files

### `fortls.config.json`
Configuration for Fortran Language Server (VS Code integration)

### `.envrc`
Environment setup for direnv (optional, for automatic environment activation)

---

## 🐛 Troubleshooting

### Compilation Errors

**Issue:** `gfortran: command not found`
```bash
# Install gfortran
sudo apt-get install gfortran
```

**Issue:** Module not found
```bash
# Ensure include directory is in the build path
./compile_modules.sh
```

### Runtime Issues

**Issue:** Network file not found
- Check that network files are in `./nets/` directory
- Verify file path in main program

**Issue:** Output directory missing
```bash
mkdir -p output
```

---

## 📊 Sample Workflow

```bash
# 1. Compile
./build-parallel.sh

# 3. Run simulation (output to ./output/)
./main

# 4. Analyze results
jupyter lab test.ipynb
```

---

## 🎯 Project Goals

This project was developed as a **TFG** - Final Degree Project.

**Main objectives:**
- Implement efficient epidemic simulation on complex networks
- Compare SIS and SIR models on various network topologies
- Provide tools for analyzing disease spreading dynamics
- Demonstrate high-performance computing techniques (Fortran + OpenMP)

---

## 📄 License

---

## 👤 Author

**Adrià Rojo**
- GitHub: [@AdriRed](https://github.com/AdriRed)
- Repository: [epidemics-tfg](https://github.com/AdriRed/epidemics-tfg)

---

## 📚 References

- **Epidemic Models:** [SIS and SIR models](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology)
- **Complex Networks:** Network science and graph theory
- **Fortran Performance:** [Modern Fortran Guide](https://fortran-lang.org/)

---

## ✨ Key Technologies

| Component | Technology |
|-----------|-----------|
| **Core Engine** | Fortran 90 |
| **Parallelization** | OpenMP |
| **Data Structures** | Hash Maps, Skip Lists |
| **Editor** | VS Code |
