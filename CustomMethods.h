#ifndef CUSTOMMETHODS_H
#define CUSTOMMETHODS_H

#include <vector>
#include <cmath>
#include <limits>
#include <set>
#include <algorithm>
#include <random>  // Added to fix random generation
#include <iostream>  // For printing results

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
    Delaunay(const std::vector<Point> &points) : points(points) {
        triangulate();
    }

    std::vector<Point> points;
    std::vector<std::vector<Point>> neighbors;
    std::vector<std::vector<Point>> queryTable;

    void triangulate() {
        size_t n = points.size();
        neighbors.resize(n);
        queryTable.resize(n);

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    neighbors[i].push_back(points[j]);
                }
            }
        }
    }

    void updateQueryTable(int m) {
        std::set<int> adjacentPoints;
        for (size_t i = 0; i < m; ++i) {
            if (std::find(neighbors[m].begin(), neighbors[m].end(), points[i]) != neighbors[m].end()) {
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

    int findNearestPointIndex(const Point &pk) {
        double minDist = std::numeric_limits<double>::max();
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
    FPS(const std::vector<Point>& points) : points(points) {}

    std::vector<Point> farthestPointSample(size_t sampleSize) {
        std::vector<Point> sampledPoints;
        std::vector<bool> selected(points.size(), false);

        sampledPoints.push_back(points[0]);
        selected[0] = true;

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
    std::vector<Point> points;

    double closestDistance(const Point& point, const std::vector<Point>& sampledPoints) {
        double minDist = std::numeric_limits<double>::max();
        for (const Point& p : sampledPoints) {
            minDist = std::min(minDist, point.distanceTo(p));
        }
        return minDist;
    }
};

// Nearest Neighbor Search Algorithm
class NNSAlgorithm {
public:
    NNSAlgorithm(const std::vector<Point>& points) : points(points) {}

    Point findNearestNeighbor(const Point& queryPoint) {
        double minDist = std::numeric_limits<double>::max();
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
    std::vector<Point> points;
};

// Generate random points
std::vector<Point> generateRandomPoints(size_t numPoints) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    std::vector<Point> points;
    for (size_t i = 0; i < numPoints; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        points.push_back(Point(x, y, 0)); // 2D points for simplicity
    }
    return points;
}

// KD-Tree Node
struct KDNode {
    Point point;
    KDNode *left, *right;

    KDNode(Point p) : point(p), left(nullptr), right(nullptr) {}
};

// KD-Tree
class KDTree {
public:
    KDTree(const std::vector<Point>& points) {
        root = build(points, 0);
    }

    Point nearestNeighbor(const Point& target) {
        return nearest(root, target, 0).point;
    }

private:
    KDNode* root;

    KDNode* build(std::vector<Point> points, int depth) {
        if (points.empty()) return nullptr;

        int axis = depth % 3;
        std::sort(points.begin(), points.end(), [axis](const Point& a, const Point& b) {
            return (axis == 0) ? a.x < b.x : (axis == 1) ? a.y < b.y : a.z < b.z;
        });

        int mid = points.size() / 2;
        KDNode* node = new KDNode(points[mid]);
        std::vector<Point> leftPoints(points.begin(), points.begin() + mid);
        std::vector<Point> rightPoints(points.begin() + mid + 1, points.end());

        node->left = build(leftPoints, depth + 1);
        node->right = build(rightPoints, depth + 1);

        return node;
    }

    struct KDResult {
        Point point;
        double distance;
    };

    KDResult nearest(KDNode* node, const Point& target, int depth) {
        if (!node) return {Point(), std::numeric_limits<double>::max()};

        int axis = depth % 3;
        KDNode* nextBranch = nullptr;
        KDNode* otherBranch = nullptr;

        if ((axis == 0 && target.x < node->point.x) ||
            (axis == 1 && target.y < node->point.y) ||
            (axis == 2 && target.z < node->point.z)) {
            nextBranch = node->left;
            otherBranch = node->right;
        } else {
            nextBranch = node->right;
            otherBranch = node->left;
        }

        KDResult best = nearest(nextBranch, target, depth + 1);
        double dist = target.distanceTo(node->point);
        if (dist < best.distance) {
            best = {node->point, dist};
        }

        double splitDist = (axis == 0) ? fabs(target.x - node->point.x) :
                           (axis == 1) ? fabs(target.y - node->point.y) :
                                         fabs(target.z - node->point.z);

        if (splitDist < best.distance) {
            KDResult otherBest = nearest(otherBranch, target, depth + 1);
            if (otherBest.distance < best.distance) {
                best = otherBest;
            }
        }

        return best;
    }
};

// R-Tree-like implementation (brute-force for simplicity)
class RTree {
public:
    RTree(const std::vector<Point>& points) : points(points) {}

    Point nearestNeighbor(const Point& target) {
        double minDist = std::numeric_limits<double>::max();
        Point nearest;

        for (const auto& point : points) {
            double dist = target.distanceTo(point);
            if (dist < minDist) {
                minDist = dist;
                nearest = point;
            }
        }

        return nearest;
    }

private:
    std::vector<Point> points;
};

#endif // CUSTOMMETHODS_H
