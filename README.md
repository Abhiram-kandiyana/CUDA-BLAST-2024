# CUDA-BLAST-2024
### This repository is to compile [NCBI BLAST](https://github.com/TadasBaltrusaitis/OpenFace) with CUDA backend. The repo implements both `k-mer filtering` and `Smith-Waterman Alignment` processes. 

## **CUDA-BLAST-2024 Installation**
* Follow the below steps to clone this repo:
```
git clone https://github.com/Abhiram-kandiyana/CUDA-BLAST-2024.git
cd CUDA-BLAST-2024
```
## **Dependency Installation**
* Install CUDA to your machine so you can compile the .cu files. Follow this guide for [Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/) and [Windows](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/).

## **Running CUDA-BLAST-2024**
1. The SRC directory has all the source code (C, CUDA). Move to the SRC directory
```
cd SRC
```
2. It has both the CPU and GPU implementation of BLAST. The CPU implementation is not as fast the GPU implementation. If you are unfamiliar with CUDA, go through the CPU implementation first and then try to read through the GPU implementation.

* To run the CPU (C) CODE:
  ```
  cd CPU_codes
  
  g++ -std=c++17 -o kmer_cpu_final kmer_cpu_final.cpp -pthread
  
  ./kmer_cpu_final
  ```
* To run the GPU (CUDA) CODE:
  ```
  cd GPU_codes
  
  nvcc kmer_gpu_final.cu -o kmer_gpu_final
  
  ./kmer_gpu_final
  ```

## **Results**

| Method                | K-mer Length | Duration (in seconds) | Throughput (queries/sec) |
|-----------------------|--------------|-----------------------|--------------------------|
| Multithreaded CPU     | 15           | 10.014                | 49.930                   |
| CUDA-BLAST-2024 (Ours)| 15           | **2.682**             | **186.428**              |
| NCBI BLAST            | 15           | 11.000                | 43.478                   |
| Multithreaded CPU     | 10           | 10.458                | 47.811                   |
| CUDA-BLAST-2024 (Ours)| 10           | **2.909**             | **171.880**              |
| NCBI BLAST            | 10           | 14.000                | 43.478                   |

* The best results for each K-mer length are bolded. 
* CUDA-BLAST-2024 achieves a **5x** improvment in speed over NCBI-BLAST.



