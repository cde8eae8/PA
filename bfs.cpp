#include <vector>
#include <iomanip>
#include <atomic>
#include <iostream>
#include <optional>
#include <cstdint>
#include <queue>
#include <assert.h>
#include "tbb/tbb.h"

struct Node {
    std::vector<uint64_t> edges;
    std::atomic<int> distance = -1;
    int parent = -1;
};

struct Graph {
    Graph(size_t size) : size{size}, nodes(size * size * size) { }

    bool valid_id(int i, int j, int k) {
        return !(i >= size || j >= size || k >= size || i < 0 || j < 0 || k < 0);
    }

    int get_id(size_t i, size_t j, size_t k) {
        return i * size * size + j * size + k;
    }

    size_t size;
    std::vector<Node> nodes;
};

void seq_bfs(Graph& graph) {
    std::queue<uint64_t> q;
    int id = graph.get_id(0, 0, 0);
    q.push(id);
    graph.nodes[id].distance = 0;
    while (!q.empty()) {
        int id = q.front();
        q.pop();
        int dist = graph.nodes[id].distance;
        for (auto e : graph.nodes[id].edges) {
            if (graph.nodes[e].distance == -1) {
                graph.nodes[e].distance = dist + 1;
                q.push(e);
            }
        }
    }
}

void par_bfs(Graph& graph) {
    struct Front {
        uint64_t id;
    };

    std::unique_ptr<Front[]> frontPtr = std::make_unique<Front[]>(graph.nodes.size());
    auto firstChildPosLocal = std::make_unique<uint64_t[]>(graph.nodes.size());
    auto firstChildPos = std::make_unique<uint64_t[]>(graph.nodes.size());
    std::unique_ptr<Front[]> front2Ptr = std::make_unique<Front[]>(graph.nodes.size());
    Front* front = frontPtr.get();
    Front* front2 = front2Ptr.get();

    size_t front_size = 1;
    front[0].id = 0;
    graph.nodes[0].distance = 0;

    while (front_size != 0) {
        size_t dist = graph.nodes[front[0].id].distance + 1;
        tbb::parallel_for(
             tbb::blocked_range<int>(0, front_size),
             [&](tbb::blocked_range<int> r) {
                for (int i = r.begin(); i != r.end(); ++i) {
                    uint64_t id = front[i].id;
                    size_t idx = 0;
                    for (auto e : graph.nodes[id].edges) {
                        Node& node = graph.nodes[e];
                        int v = -1;
                        if (node.distance.compare_exchange_strong(v, dist)) {
                            node.parent = id;
                            idx++;
                        }
                    }
                    firstChildPosLocal.get()[i] = idx;
                }
            });
        size_t new_front_size = tbb::parallel_scan(
                    tbb::blocked_range<int>(0, front_size),
                    0,
                    [&](const tbb::blocked_range<int>& r, int sum, bool is_final_scan) {
                        for (int i = r.begin(); i != r.end(); ++i) {
                            if (is_final_scan) {
                                firstChildPos.get()[i] = sum;
                            }
                            sum += firstChildPosLocal.get()[i];
                        }
                        return sum;
                    },
                    [](int left, int right) {
                        return left + right;
                    });
        tbb::parallel_for(
             tbb::blocked_range<int>(0, front_size),
             [&](tbb::blocked_range<int> r) {
                for (int i = r.begin(); i != r.end(); ++i) {
                    uint64_t id = front[i].id;
                    size_t start = firstChildPos.get()[i];
                    int idx = 0;
                    for (auto e : graph.nodes[id].edges) {
                        if (graph.nodes[e].parent == id) {
                            front2[start + idx++].id = e;
                        }
                    }
                }
            });

        std::swap(front, front2);
        front_size = new_front_size;
    }
}

bool check(Graph& graph) {
    bool ok = true;
    for (int i = 0; i < graph.size; ++i) {
        for (int j = 0; j < graph.size; ++j) {
            for (int k = 0; k < graph.size; ++k) {
                 if (graph.nodes[graph.get_id(i, j, k)].distance != i + j + k) {
                     std::cout << "failed for " << i << " " << j << " " << k << " " <<
                                  graph.nodes[graph.get_id(i, j, k)].distance << std::endl;
                     ok = false;
                 }
            }
        }
    }
    return ok;
}

struct BenchInfo {
    uint64_t lastDelta{};
    uint64_t total{};
};

struct ns {
    ns(uint64_t v) : value{v} { }

    uint64_t value;
};

std::ostream& operator<<(std::ostream& out, ns ns) {
    return out << ns.value / 1000000000 << "." << std::setw(9) << std::setfill('0') << ns.value % 1000000000;
}

template <typename T>
void bench(BenchInfo& b, T const& f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto finish = std::chrono::high_resolution_clock::now();
    b.lastDelta = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();
}

void clearGraph(Graph &g) {
    for (auto& node : g.nodes) {
        node.distance = -1;
        node.parent = -1;
    }
}

int main() {
    size_t size = 400;
    size_t iters = 5;
    std::cout << "arguments: " << size << std::endl;

    std::unordered_map<std::string, BenchInfo> deltas;
    deltas.insert({"par", BenchInfo{}});
    deltas.insert({"seq", BenchInfo{}});

    Graph graph(size);

    auto add_edge = [&graph](int from, int i, int j, int k) {
        if (!graph.valid_id(i, j, k)) {
            return;
        }
        size_t edge = graph.get_id(i, j, k);
        graph.nodes[from].edges.push_back(edge);
    };

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            for (int k = 0; k < size; ++k) {
                int id = graph.get_id(i, j, k);
                Node& node = graph.nodes[id];

                for (int di : {-1, 1}) {
                    add_edge(id, i + di, j, k);
                    add_edge(id, i, j + di, k);
                    add_edge(id, i, j, k + di);
                }
            }
        }
    }

    for (int iter = 0; iter < iters; ++iter) {
        clearGraph(graph);

        bench(deltas.at("par"),
                [&graph]() {
                    par_bfs(graph);
                });
        if (!check(graph)) {
            std::abort();
        }

        clearGraph(graph);

        bench(deltas.at("seq"),
                [&graph]() {
                    seq_bfs(graph);
                });
        if (!check(graph)) {
            std::abort();
        }

        for (auto&& p : deltas) {
            p.second.total += p.second.lastDelta;
        }
        std::cout << "iter " << iter << ":\n";
        for (auto&& p : deltas) {
            std::cout << "\t" << p.first << ": " << ns{p.second.lastDelta} << std::endl;
        }
    }
    std::cout << "Total:\n";
    for (auto p : deltas) {
        std::cout << p.first << ": " << ns{p.second.total / iters} << std::endl;
    }
}
