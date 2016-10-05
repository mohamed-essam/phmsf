// Copyright 2016 Mohamed Essam
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>

#define PI 3.141569

struct Edge{
    unsigned int from, to, weight;
};

struct Graph{
    Edge* edges;
    int edgeCount;
};

// Graph Creation
int square(int a) {
    return a*a;
}

int getEuclidianDistance(int r1, int g1, int b1, int r2, int g2, int b2) {
    return 4 * sqrt(static_cast<double>(square(r1 - r2) +
     square(g1 - g2) + square(b1 - b2)));
}

Graph* createGraph(unsigned char* rgb, unsigned int maxWeight,
    unsigned int height, unsigned int width) {
    unsigned int edgeCount = 2 * width * height - width - height;
    Edge* graph = new Edge[edgeCount];
    unsigned int nextEdge = 0;
    unsigned int histogram[1800];
    memset(histogram, 0, sizeof(histogram));
    edgeCount = 0;
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            unsigned int idx1, idx2, idx3;
            idx1 = (i*width + j) * 3;
            idx3 = ((i + 1)*width + j) * 3;
            if (i + 1 < height) {
                unsigned int w = getEuclidianDistance(rgb[idx1], rgb[idx1 + 1],
                    rgb[idx1 + 2], rgb[idx3], rgb[idx3 + 1], rgb[idx3 + 2]);
                if (w <= maxWeight) {
                    graph[nextEdge++] = { idx1 / 3, idx3 / 3, w };
                    histogram[w]++;
                    edgeCount++;
                }
            }
            idx2 = (i*width + j + 1) * 3;
            if (j + 1 < width) {
                unsigned int w = getEuclidianDistance(rgb[idx1], rgb[idx1 + 1],
                    rgb[idx1 + 2], rgb[idx2], rgb[idx2 + 1], rgb[idx2 + 2]);
                if (w <= maxWeight) {
                    graph[nextEdge++] = { idx1 / 3, idx2 / 3, w };
                    histogram[w]++;
                    edgeCount++;
                }
            }
        }
    }
    for (int i = 1; i < 1800; i++) {
        histogram[i] += histogram[i - 1];
    }
    Edge* sortedGraph = new Edge[edgeCount];
    int st_pos = 0;
    for (int i = 0; i < edgeCount; i++) {
        int w = graph[i].weight;
        if (w == 0) {
            sortedGraph[st_pos++] = graph[i];
        } else {
            sortedGraph[histogram[w - 1]++] = graph[i];
        }
    }
    delete[] graph;
    Graph* ret = new Graph();
    ret->edgeCount = edgeCount;
    ret->edges = sortedGraph;
    return ret;
}

// Initial Regions
struct UnionFind{
    unsigned int* parent;
    unsigned int* rank;

    explicit UnionFind(unsigned int size) {
        parent = new unsigned int[size];
        rank = new unsigned int[size];
        for (unsigned int i = 0; i < size; i++) {
            rank[i] = 1;
            parent[i] = i;
        }
    }

    ~UnionFind() {
        delete[] parent;
        delete[] rank;
    }

    unsigned int find(unsigned int a) {
        return parent[a] = (a == parent[a] ? a : find(parent[a]));
    }

    unsigned int unite(unsigned int a, unsigned int b) {
        a = find(a), b = find(b);
        if (a == b) return a;
        if (rank[b] > rank[a]) std::swap(a, b);
        parent[b] = a;
        rank[a] += rank[b];
        return a;
    }

    bool isRepresetative(unsigned int x) {
        return x == find(x);
    }

    bool isUnited(unsigned int a, unsigned int b) {
        return find(a) == find(b);
    }
};

struct Region{
    unsigned int representative, size, credit;
};

struct RegionalData {
    UnionFind* uf;
    Region* regions;

    RegionalData(UnionFind* _uf, Region* _regions) {
        uf = _uf;
        regions = _regions;
    }
};

int computeCredit(int regionSize) {
    return sqrt(4 * PI * regionSize);
}

RegionalData* getInitialRegions(Edge* graph, int minWeight, int minRegionSize,
    int edgeCount, int width, int height) {
    UnionFind* uf = new UnionFind(width*height);
    for (int i = 0; i < edgeCount && graph[i].weight <= minWeight; i++) {
        if (graph[i].weight > minWeight)
            continue;
        uf->unite(graph[i].from, graph[i].to);
    }
    Region* regions = new Region[width*height];
    for (int i = 0; i < width*height; i++) {
        regions[i].size = regions[i].representative = regions[i].credit = 0;
    }
    for (int i = 0; i < width*height; i++) {
        regions[uf->find(i)].size++;
        regions[uf->find(i)].representative = uf->find(i);
    }
    for (int i = 0; i < width*height; i++) {
        if (regions[i].size >= minRegionSize)
            regions[i].credit = computeCredit(regions[i].size);
        else
            regions[i].credit = 1e9;
    }
    RegionalData* ret = new RegionalData(uf, regions);
    return ret;
}

// Region Expansion

struct SegmentationData{
    int* segment;
    int segmentCount;
};

SegmentationData expandRegions(RegionalData* rd, int minWeight,
    int width, int height, Edge* graph, int edgeCount) {
    for (int i = 0; i < edgeCount; i++) {
        if (graph[i].weight <= minWeight) continue;
        if (!rd->uf->isUnited(graph[i].from, graph[i].to)) {
            int a = rd->uf->find(graph[i].from);
            int b = rd->uf->find(graph[i].to);
            int credit = std::min(rd->regions[a].credit, rd->regions[b].credit);
            if (credit > graph[i].weight) {
                int s = rd->uf->unite(a, b);
                rd->regions[s].credit = credit - graph[i].weight;
            }
        }
    }
    int* repTable = new int[width*height];
    int idx = 0;
    for (int i = 0; i < width*height; i++) {
        if (rd->uf->isRepresetative(i)) {
            repTable[i] = idx++;
        }
    }
    int* regionals = new int[width*height];
    for (int i = 0; i < width*height; i++) {
        regionals[i] = repTable[rd->uf->find(i)];
    }
    SegmentationData ret = { regionals, idx };
    return ret;
}

// Main entry point

SegmentationData segmentImage(unsigned char* rgb, int height,
    int width, int minWeight, int maxWeight, int minRegionSize) {
    Graph* graph = createGraph(rgb, maxWeight, height, width);
    RegionalData* regional = getInitialRegions(graph->edges, minWeight,
        minRegionSize, graph->edgeCount, width, height);
    return expandRegions(regional, minWeight, width, height,
        graph->edges, graph->edgeCount);
}
