// Copyright 2016 Mohamed Essam
#include <math.h>
#include <algorithm>

#define PI 3.141569

struct Edge{
  unsigned int from, to, weight;
};

// Graph Creation
int square(int a) {
    return a*a;
}

int getEuclidianDistance(int r1, int g1, int b1, int r2, int g2, int b2) {
    return sqrt(square(r1-r2)+square(g1-g2)+square(b1-b2));
}

Edge* createGraph(unsigned char* rgb, unsigned int maxWeight,
        unsigned int height, unsigned int width) {
    unsigned int edgeCount = 2 * width * height - w - h;
    Edge* graph = new Edge[edgeCount];
    unsigned int nextEdge = 0;
    unsigned int histogram[450];
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            unsigned int idx1, idx2, idx3;
            idx1 = (i*width+j)*3;
            idx2 = (i*width+j+1)*3;
            if (j+1 < width) {
                unsigned int w = getEuclidianDistance(rgb[idx1], rgb[idx1+1],
                 rgb[idx1+2], rgb[idx2], rgb[idx2+1], rgb[idx2+2]);
                if (w > maxWeight) {
                    edgeCount--;
                } else {
                    graph[nextEdge++] = {idx1/3, idx2/3, w}
                    histogram[w]++;
                }
            }
            idx3 = ((i+1)*width+j)*3;
            if (i+1 < width) {
                unsigned int w = getEuclidianDistance(rgb[idx1], rgb[idx1+1],
                 rgb[idx1+2], rgb[idx3], rgb[idx3+1], rgb[idx3+2]);
                 if (w > maxWeight) {
                     edgeCount--;
                 } else {
                    graph[nextEdge++] = {idx1/3, idx3/3, w}
                    histogram[w]++;
                }
            }
        }
    }
    for (int i = 1; i < 450; i++) {
        histogram[i] += histogram[i-1];
    }
    Edge* sortedGraph = new Edge[edgeCount];
    int st_pos = 0;
    for (int i = 0; i < edgeCount; i++) {
        int w = graph[i].weight;
        if (w == 0) {
            sortedGraph[st_pos++] = graph[i];
        } else {
            sortedGraph[histogram[w-1]++] = graph[i];
        }
    }
    delete[] graph;
    return sortedGraph;
}

// Initial Regions
struct UnionFind{
    unsigned int* parent;
    unsigned int* rank;

    explicit UnionFind(unsigned int size) {
        parent = new unsigned int[size];
        rank = new unsigned int[size];
        for (unsigned int i = 0; i < size; i++) {
            rank[i] = 0;
            parent[i] = i;
        }
    }

    ~UnionFind() {
        delete[] parent;
    }

    unsigned int find(unsigned int a) {
        return a == parent[a] ? a : find(parent[a]);
    }

    unsigned int unite(unsigned int a, unsigned int b) {
        a = find(a), b = find(b);
        if (a == b) return a;
        if (rank[b] > rank[a]) swap(a, b);
        else if (rank[a] == rank[b]) rank[a]++;
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
};

int computeCredit(int regionSize) {
    return sqrt(4 * PI * regionSize);
}

RegionalData* getInitialRegions(Edge* graph, int minWeight, int minRegionSize,
        int edgeCount, int width, int height) {
    UnionFind* uf = new UnionFind(width*height);
    for (int i = 0; i < edgeCount && graph[i].weight <= minWeight; i++) {
        if (!uf.isUnited(graph[i].from, graph[i].to)) {
            uf.unite(graph[i].from, graph[i].to);
        }
    }
    Region* regions = new Region[width*height];
    for (int i = 0; i < width*height; i++) {
        regions[uf.find(i)].size++;
        regions[uf.find(i)].representative = uf.find(i);
    }
    for (int i = 0; i < width*height; i++) {
        if (region[i].size >= minRegionSize)
            regions[i].credit = computeCredit(regions[i].size);
        else
            regions[i].credit = 1e9;
    }
    RegionalData* ret = new RegionalData(uf, regions);
    return ret;
}
