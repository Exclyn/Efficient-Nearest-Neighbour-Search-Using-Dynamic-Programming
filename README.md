# Efficient Nearest Neighbour Search Using Dynamic Programming

This project explores different algorithms for **nearest neighbour search** in multi-dimensional spaces, with an emphasis on efficiency and experimental comparison.  
It provides implementations of several approaches and prepares the ground for a **Dynamic Programming (DP) based query-table method**.

---

##  Features
- **Brute Force Search** – simple baseline for nearest neighbour queries.
- **Farthest Point Sampling (FPS)** – subsampling heuristic to reduce search space.
- **KD-Tree** – recursive space-partitioning structure for efficient queries.
- **R-Tree (stub)** – hierarchical bounding-box search (basic version provided).
- **Delaunay / Query Table (experimental)** – placeholder for DP-based approach.

---

##  Project Structure
├── CustomMethods.h # Core data structures and algorithms (Point, KDTree, etc.)
├── newupdate.cpp # Interactive demo (build structures, query a point)
├── results.cpp # Benchmarking & timing experiments

---

##  Build Instructions

### Prerequisites
- C++17 or later  
- CMake (optional, can also use g++ directly)

### Compile & Run
# Compile the benchmarking program
g++ -std=c++17 results.cpp -o results

# Compile the interactive demo
g++ -std=c++17 newupdate.cpp -o demo
Run with:

./results   # Runs benchmarks (KD-tree, brute force, etc.)
./demo      # Interactive demo with user input

## Usage Examples
Benchmark (results.cpp)
Generates random points, builds data structures, and compares:

  Build times
  
  Nearest neighbour query times
  Outputs results in the terminal.

3Example output:


Nearest neighbour search with 10000 points:
Brute force time: 1.23 ms
KD-tree build + query time: 0.09 ms
