#include <iostream>
#include <vector>
#include <chrono>  // For measuring time
#include "CustomMethods.h"

using namespace std;
using namespace chrono;

int main() {
    // Generate random points
    size_t numPoints;
    cout << "Enter the number of points to generate: ";
    cin >> numPoints;

    vector<Point> points = generateRandomPoints(numPoints);

    // Ask for the query point after displaying prompt
    double qx, qy, qz = 0;
    cout << "Enter the query point (x, y, z for 3D or x, y for 2D): ";
    cin >> qx >> qy;
    if (cin.peek() != '\n') {
        cin >> qz;
    }

    Point query(qx, qy, qz);

    // 1. Brute-force Nearest Neighbor Search
    NNSAlgorithm nns(points);
    auto start = high_resolution_clock::now();
    Point nearestBruteForce = nns.findNearestNeighbor(query);
    auto end = high_resolution_clock::now();
    auto durationBF = duration_cast<microseconds>(end - start);
    cout << "\nNearest Point using NNS: (" 
         << nearestBruteForce.x << ", " 
         << nearestBruteForce.y << ", " 
         << nearestBruteForce.z << ")\n";
    cout << "Time taken by NSS: " 
         << durationBF.count() << " microseconds.\n";

    // 2. KD-Tree Nearest Neighbor Search
    auto startKD = high_resolution_clock::now();
    KDTree kdTree(points);
    Point nearestKDTree = kdTree.nearestNeighbor(query);
    auto endKD = high_resolution_clock::now();
    auto durationKD = duration_cast<microseconds>(endKD - startKD);
    cout << "\nNearest Point using KD-Tree: (" 
         << nearestKDTree.x << ", " 
         << nearestKDTree.y << ", " 
         << nearestKDTree.z << ")\n";
    cout << "Time taken by KD-Tree method: " 
         << durationKD.count() << " microseconds.\n";

    // 3. R-Tree Nearest Neighbor Search
    auto startRTree = high_resolution_clock::now();
    RTree rTree(points);
    Point nearestRTree = rTree.nearestNeighbor(query);
    auto endRTree = high_resolution_clock::now();
    auto durationRTree = duration_cast<microseconds>(endRTree - startRTree);
    cout << "\nNearest Point using R-Tree: (" 
         << nearestRTree.x << ", " 
         << nearestRTree.y << ", " 
         << nearestRTree.z << ")\n";
    cout << "Time taken by R-Tree method: " 
         << durationRTree.count() << " microseconds.\n";

    return 0;
}
