# **Optimal and Fast Wavelength Selection for Spectral Unmixing in Photoacoustic Imaging**

This repository contains MATLAB tools for generating synthetic absorption spectra, organizing real chemical spectra, and implementing algorithms for optimal wavelength selection in spectral unmixing. 

For additional details, see in-code documentation and comments within each file.

---

##  Algorithms 

We include three core algorithms for identifying optimal wavelength subsets. Each algorithm aims to find a set of wavelengths that improves spectral unmixing performance (e.g., minimizes the Frobenius norm of the pseudoinverse or maximizes the smallest singular value).

### **1. Brute Force Algorithm**
Evaluates **all possible combinations** of wavelengths and selects the subset that minimizes  the Frobenius norm of the pseudoinverse.  
This provides the *globally optimal* solution but becomes expensive for large `k` or `n`.

---

### **2. LUKE Algorithm**
A greedy algorithm that iteratively removes the wavelength (column) whose removal **maximizes the smallest singular value** of the remaining matrix.  
This provides a fast and effective approximation to the optimal set.

---

### **3. NK Algorithm**
Re-parameterizes the wavelength-selection problem into a **binary vector NK-landscape** and performs a greedy search over the space, keeping the highest-performing configurations.  
Useful for large search spaces where brute force is not feasible.

---

##  Generators

We include code for both synthetic spectrum generation and real spectral dataset construction.

### **Synthetic Spectra**
Used to test algorithm behavior under controlled conditions.

- **`build_curve.m`**  
  Generates a synthetic absorption spectrum using a sum of Gaussian peaks.  
  Configurable parameters include:
  - number of peaks  
  - wavelength range  
  - number of species  
  - peak amplitudes and widths

- **`generate_spectrum_curve.m`**  
  Calls `build_curve` repeatedly to form a full spectral matrix  
  (each row = one synthetic species).

---

### **Real Spectra**

Real spectral datasets are stored in the **`Spectrum_Data/`** directory as `.csv` or `.mat` files.

Key utilities:

- **`pick_bandwidth.m`**  
  Crops a spectrum to a desired wavelength range and interpolates to a fixed number of sampling points.  
  Useful for normalizing heterogeneous datasets.

- **`build_absorption_matrix.m`**  
  Constructs a 5-species absorption matrix from individual CSV files.  
  *(Paths may need to be updated based on your environment.)*

---

##  Examples

We include several example scripts demonstrating how to use the generators and algorithms:

- **`example_build_spectrum.m`**  
  Demonstrates how to crop and interpolate real spectral curves using `pick_bandwidth`.

- **`example_generate_synthetic_curve.m`**  
  Shows how to generate synthetic spectra using `build_curve` and `generate_spectrum_curve`.


