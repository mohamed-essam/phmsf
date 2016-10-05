// Copyright 2016 Mohamed Essam
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <queue>

#define PI 3.141569

struct Edge{
    unsigned int from, to, weight;
};

struct Tile{
    std::vector<Edge> edges;
    std::queue<Edge> tileQ;
    unsigned int edgeCount;
    int histogram[1800];
    int start, end;
};

int square(int a) {
    return a*a;
}

int getEuclidianDistance(int r1, int g1, int b1, int r2, int g2, int b2) {
    return 4 * sqrt(static_cast<double>(square(r1 - r2) +
     square(g1 - g2) + square(b1 - b2)));
}

struct UnionFind{
    unsigned int* parent;
    unsigned int* rank;
    bool* marked;

    explicit UnionFind(unsigned int size) {
        parent = new unsigned int[size];
        rank = new unsigned int[size];
        marked = new bool[size];
        for (unsigned int i = 0; i < size; i++) {
            rank[i] = 1;
            parent[i] = i;
            marked[i] = false;
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

bool isInside(int a, int b, int c) {
    return a >= b && a <= c;
}

double computeCredit(int size) {
    return sqrt(4 * PI * size);
}

int* constructGraph(unsigned char* rgb, unsigned int height, unsigned int width,
    unsigned int minWeight, unsigned int maxWeight, int& segmentCount) {
    int threadCount = omp_get_max_threads();
    Tile* tiles = new Tile[threadCount];
    for (int i = 0; i < threadCount; i++) {
        tiles[i].edgeCount = tiles[i].end = 0;
        tiles[i].start = width*height;
        memset(tiles[i].histogram, 0, sizeof(tiles[i].histogram));
    }
#pragma omp parallel for schedule(static, width*height / omp_get_max_threads())
    for (int i = 0; i < width*height; i++) {
        int th = omp_get_thread_num();
        tiles[th].start = std::min(tiles[th].start, i);
        tiles[th].end = std::max(tiles[th].end, i);
        int j = i + 1;
        unsigned int w1 = getEuclidianDistance(rgb[i * 3], rgb[i * 3 + 1],
         rgb[i * 3 + 2], rgb[j * 3], rgb[j * 3 + 1], rgb[j * 3 + 2]);
        if (w1 <= maxWeight && j % width > 0) {
            tiles[th].edgeCount++;
            tiles[th].edges.push_back({ i, j, w1 });
            tiles[th].histogram[w1]++;
        }
        int k = i + width;
        unsigned int w2 = getEuclidianDistance(rgb[i * 3], rgb[i * 3 + 1],
         rgb[i * 3 + 2], rgb[k * 3], rgb[k * 3 + 1], rgb[k * 3 + 2]);
        if (w2 <= maxWeight && k / width < height) {
            tiles[th].edgeCount++;
            tiles[th].edges.push_back({ i, k, w2 });
            tiles[th].histogram[w2]++;
        }
    }
#pragma omp parallel
    {
        int th = omp_get_thread_num();
        for (int i = 1; i < 1800; i++) {
            tiles[th].histogram[i] += tiles[th].histogram[i - 1];
        }
        int st_num = 0;
        std::vector<Edge> temp(tiles[th].edgeCount);
        for (int i = 0; i < tiles[th].edgeCount; i++) {
            int w = tiles[th].edges[i].weight;
            if (w == 0) {
                temp[st_num++] = tiles[th].edges[i];
            } else {
                temp[tiles[th].histogram[w - 1]++] = tiles[th].edges[i];
            }
        }
        tiles[th].edges.clear();
        tiles[th].edges.shrink_to_fit();
        tiles[th].edges = temp;
    }
    UnionFind* uf = new UnionFind(width*height);
#pragma omp parallel
    {
        int th = omp_get_thread_num();
        for (int i = 0; i < tiles[th].edgeCount &&
         tiles[th].edges[i].weight <= minWeight; i++) {
            if (isInside(tiles[th].edges[i].from,
             tiles[th].start, tiles[th].end)
                && isInside(tiles[th].edges[i].to,
                 tiles[th].start, tiles[th].end)) {
                uf->unite(tiles[th].edges[i].from, tiles[th].edges[i].to);
            } else {
                int survivor = uf->unite(tiles[th].edges[i].from,
                 tiles[th].edges[i].to);
                uf->marked[survivor] = true;
            }
        }
    }
    unsigned int* minEdgeWeight = new unsigned int[width*height];
#pragma omp parallel
    {
        int th = omp_get_thread_num();
        for (int i = 0; i < tiles[th].edgeCount; i++) {
            if (!uf->isUnited(tiles[th].edges[i].from, tiles[th].edges[i].to)) {
                minEdgeWeight[uf->find(tiles[th].edges[i].from)] =
                std::min(minEdgeWeight[uf->find(tiles[th].edges[i].from)],
                 tiles[th].edges[i].weight);
                minEdgeWeight[uf->find(tiles[th].edges[i].to)] =
                 std::min(minEdgeWeight[uf->find(tiles[th].edges[i].to)],
                  tiles[th].edges[i].weight);
            }
        }
    }
    int* credit = new int[width*height];
#pragma omp parallel for
    for (int i = 0; i < width*height; i++) {
        if (uf->isRepresetative(i)) {
            credit[i] = computeCredit(uf->rank[i]) * minEdgeWeight[i];
        }
    }
#pragma omp parallel
    {
        int th = omp_get_thread_num();
        for (int i = 0; i < tiles[th].edgeCount; i++) {
            if (!(isInside(tiles[th].edges[i].from,
             tiles[th].start, tiles[th].end) &&
              isInside(tiles[th].edges[i].to,
               tiles[th].start, tiles[th].end))) {
                uf->marked[uf->find(tiles[th].edges[i].from)] =
                uf->marked[uf->find(tiles[th].edges[i].from)] = true;
            } else if (uf->marked[uf->find(tiles[th].edges[i].from)] &&
             uf->marked[uf->find(tiles[th].edges[i].to)]) {
                tiles[th].tileQ.push(tiles[th].edges[i]);
            } else {
                if (!uf->isUnited(tiles[th].edges[i].from,
                 tiles[th].edges[i].to)) {
                    int cr = std::min(credit[uf->find(tiles[th].edges[i].from)],
                     credit[uf->find(tiles[th].edges[i].to)]);
                    if (cr > tiles[th].edges[i].weight) {
                        int sur = uf->unite(uf->find(tiles[th].edges[i].from),
                         uf->find(tiles[th].edges[i].to));
                        credit[sur] = cr - tiles[th].edges[i].weight;
                        uf->marked[sur] =
                        uf->marked[uf->find(tiles[th].edges[i].from)] |
                         uf->marked[uf->find(tiles[th].edges[i].to)];
                    }
                }
            }
        }
    }
    for (int i = 0; i < threadCount; i++) {
        while (tiles[i].tileQ.size()) {
            if (!uf->isUnited(tiles[i].tileQ.front().from,
             tiles[i].tileQ.front().to)) {
                int cr = std::min(
                credit[uf->find(tiles[i].tileQ.front().from)],
                 credit[uf->find(tiles[i].tileQ.front().to)]);
                if (cr > tiles[i].tileQ.front().weight) {
                    int sur = uf->unite(
                    uf->find(tiles[i].tileQ.front().from),
                     uf->find(tiles[i].tileQ.front().to));
                    credit[sur] = cr - tiles[i].tileQ.front().weight;
                    uf->marked[sur] =
                    uf->marked[uf->find(tiles[i].tileQ.front().from)] |
                     uf->marked[uf->find(tiles[i].tileQ.front().to)];
                }
            }
            tiles[i].tileQ.pop();
        }
    }
    segmentCount = 0;
    int* repTable = new int[width*height];
    for (int i = 0; i < width*height; i++) {
        if (uf->isRepresetative(i)) repTable[i] = segmentCount++;
    }
    int* segment = new int[width*height];
#pragma omp parallel for
    for (int i = 0; i < width*height; i++) {
        segment[i] = repTable[uf->find(i)];
    }
    return segment;
}
