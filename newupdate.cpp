#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <algorithm>
#include <queue>
using namespace std;

// Point structure
struct Point {
    double x, y, z;

    Point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    double distanceTo(const Point &other) const {
        return sqrt((x - other.x) * (x - other.x) +
                    (y - other.y) * (y - other.y) +
                    (z - other.z) * (z - other.z));
    }

    bool operator<(const Point &other) const {
        if (x != other.x)
            return x < other.x;
        if (y != other.y)
            return y < other.y;
        return z < other.z;
    }

    bool operator==(const Point &other) const {
        return (x == other.x && y == other.y && z == other.z);
    }
};

// Delaunay structure with incremental construction
struct Delaunay {
public:
    Delaunay(const vector<Point> &points) : points(points) {
        triangulate();
    }

    vector<Point> points;
    vector<vector<Point>> neighbors;
    vector<vector<Point>> queryTable;

    void triangulate() {
        size_t n = points.size();
        neighbors.resize(n);
        queryTable.resize(n);

        // Simulating adjacency for neighbors (simplified)
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    neighbors[i].push_back(points[j]);
                }
            }
        }
    }

    void updateQueryTable(int m) {
        set<int> adjacentPoints;
        for (size_t i = 0; i < m; ++i) {
            if (find(neighbors[m].begin(), neighbors[m].end(), points[i]) != neighbors[m].end()) {
                adjacentPoints.insert(i);
            }
        }

        for (int adj : adjacentPoints) {
            queryTable[adj].push_back(points[m]);
        }
    }

    void insertAndUpdate(int m) {
        points.push_back(points[m]);
        updateQueryTable(m);
    }

    // Find the nearest point to the new point pk using the query table
    int findNearestPointIndex(const Point &pk) {
        double minDist = numeric_limits<double>::max();
        int nearestPointIndex = -1;

        for (size_t i = 0; i < points.size(); ++i) {
            double dist = pk.distanceTo(points[i]);
            if (dist < minDist) {
                minDist = dist;
                nearestPointIndex = i;
            }
        }

        return nearestPointIndex;
    }
};

// Farthest Point Sampling Class
class FPS {
public:
    FPS(const vector<Point>& points) : points(points) {}

    vector<Point> farthestPointSample(size_t sampleSize) {
        vector<Point> sampledPoints;
        vector<bool> selected(points.size(), false);
        
        // Randomly select the first point
        sampledPoints.push_back(points[0]);
        selected[0] = true;

        // Perform Farthest Point Sampling
        for (size_t i = 1; i < sampleSize; ++i) {
            double maxDist = -1;
            int farthestPointIndex = -1;
            for (size_t j = 0; j < points.size(); ++j) {
                if (!selected[j]) {
                    double dist = closestDistance(points[j], sampledPoints);
                    if (dist > maxDist) {
                        maxDist = dist;
                        farthestPointIndex = j;
                    }
                }
            }
            sampledPoints.push_back(points[farthestPointIndex]);
            selected[farthestPointIndex] = true;
        }

        return sampledPoints;
    }

private:
    vector<Point> points;

    double closestDistance(const Point& point, const vector<Point>& sampledPoints) {
        double minDist = numeric_limits<double>::max();
        for (const Point& p : sampledPoints) {
            minDist = min(minDist, point.distanceTo(p));
        }
        return minDist;
    }
};

// Nearest Neighbor Search Algorithm
class NNSAlgorithm {
public:
    NNSAlgorithm(const vector<Point>& points) : points(points) {}

    Point findNearestNeighbor(const Point& queryPoint) {
        double minDist = numeric_limits<double>::max();
        Point nearestPoint;

        for (const Point& p : points) {
            double dist = queryPoint.distanceTo(p);
            if (dist < minDist) {
                minDist = dist;
                nearestPoint = p;
            }
        }

        return nearestPoint;
    }

private:
    vector<Point> points;
};

// Main function to integrate all
int main() {
    size_t numPoints;
    cout << "Enter the number of points: ";
    cin >> numPoints;

    vector<Point> points;
    cout << "Enter the points (x, y, z for 3D or x, y for 2D):\n";

    for (size_t i = 0; i < numPoints; ++i) {
        double x, y, z = 0;
        cout << "Point " << i + 1 << ": ";
        cin >> x >> y;
        if (cin.peek() != '\n') {
            cin >> z;
        }
        points.push_back(Point(x, y, z));
    }

    // Farthest Point Sampling
    FPS fps(points);
    size_t sampleSize = numPoints / 2;  // Example: sample half the points
    vector<Point> sampledPoints = fps.farthestPointSample(sampleSize);

    cout << "Sampled Points using Farthest Point Sampling: \n";
    for (const Point& p : sampledPoints) {
        cout << "(" << p.x << ", " << p.y << ", " << p.z << ")\n";
    }

    // Incremental Delaunay Construction
    Delaunay delaunay(points);
    
    // Insert points one by one in the farthest point sampled order
    for (size_t i = 0; i < sampledPoints.size(); ++i) {
        // Find the nearest point to the current sampled point
        int nearestPointIndex = delaunay.findNearestPointIndex(sampledPoints[i]);
        cout << "Inserting Point: (" << sampledPoints[i].x << ", " << sampledPoints[i].y << ", " << sampledPoints[i].z << ")\n";
        cout << "Nearest Point: (" << delaunay.points[nearestPointIndex].x << ", " << delaunay.points[nearestPointIndex].y << ", " << delaunay.points[nearestPointIndex].z << ")\n";
        
        // Insert and update the Delaunay structure
        delaunay.insertAndUpdate(i);
    }

    // Nearest Neighbor Search
    NNSAlgorithm nns(points);
    double qx, qy, qz = 0;
    cout << "Enter the query point (x, y, z for 3D or x, y for 2D): ";
    cin >> qx >> qy;
    if (cin.peek() != '\n') {
        cin >> qz;
    }

    Point query(qx, qy, qz);
    Point nearest = nns.findNearestNeighbor(query);

    cout << "Nearest Point: (" << nearest.x << ", " << nearest.y << ", " << nearest.z << ")\n";

    return 0;
}
