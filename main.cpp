#include <iostream>
#include <fstream>
#include <future>
#include <vector>
#include <map>
#include <iomanip>
#include <random>
#include <ostream>
#include <execution>
#include <assert.h>
#include "tbb/tbb.h"

using namespace std;

std::ostream& operator<<(std::ostream& out, std::pair<int*, int*> v) {
	out << "{";
	for (int* p = v.first; p < v.second; ++p) {
		out << *p;
		if (p + 1 != v.second) {
			out << ", ";
		}
	}
	return out << "}";
}

template <typename Exec>
std::pair<int*, int*> partition(Exec exec, int* begin, int* end) {
    size_t size = end - begin;
    int* pivotPtr = begin + size / 2;
    int pivot = *pivotPtr;
    std::swap(*(pivotPtr), *(begin));
    auto it = std::partition(exec, begin + 1, end, [pivot](int value) { return value < pivot; });
    --it;
    std::swap(*begin, *it);
    return {it, it+1};
}

void qsort_serial(int* begin, int* end) {
    size_t size = end - begin;
//    std::cout << size << std::endl;
    if (size <= 1) {
        return;
    }
    auto [it, drop_eq_it] = partition(std::execution::seq, begin, end);
    qsort_serial(begin, it);
    qsort_serial(drop_eq_it, end);
}

//void qsort_serial(int* begin, int* end) {
//    size_t size = end - begin;
//	if (size <= 1) {
//		return;
//	}
//    auto i = partition(begin, 0, size - 1);
////	int* pivotPtr = begin + size / 2;
////	int pivot = *pivotPtr;
////	std::swap(*(pivotPtr), *(begin));
////    auto it = std::partition(begin + 1, end, [pivot](int value) {
////		return value <= pivot;
////	});
////	--it;
////	std::swap(*begin, *it);
////	qsort_serial(begin, it);
////	qsort_serial(it + 1, end);

//    qsort_serial(begin, begin + i);
//    qsort_serial(begin + i + 1, end);
//}

std::ofstream out("log");

void qsort(int* begin, int* end, size_t t, size_t level = 0) {
    size_t size = end - begin;
    if (size <= std::max<size_t>(1, t)) {
        qsort_serial(begin, end);
        return;
    }
    auto [it, drop_eq_it] = partition(std::execution::seq, begin, end);
//    assert(begin <= it);
//    assert(it < drop_eq_it);
//    assert(drop_eq_it <= end);

    tbb::parallel_invoke(
        [=]() { qsort(begin, it, t, level + 1); },
        [=]() { qsort(drop_eq_it, end, t, level + 1); }
    );

//    assert(a1.load());
//    assert(a2.load());
//    assert(std::is_sorted(begin, end));
}

struct BenchInfo {
    explicit BenchInfo(size_t size)
        : values(size), lastDelta{}, total{}
    {
    }

    std::vector<int> values;
    uint64_t lastDelta;
    uint64_t total;
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
    std::vector<int>& data = b.values;
    f(data.data(), data.data() + data.size());
    auto finish = std::chrono::high_resolution_clock::now();
    b.lastDelta = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count();
}

int main() {
    size_t size = 100000000;
    std::default_random_engine g{};
    std::uniform_int_distribution<int> d{0, 1000000};

    size_t iters = 5;

    std::unordered_map<std::string, BenchInfo> deltas;
    deltas.insert({"par", BenchInfo{size}});
    deltas.insert({"ser", BenchInfo{size}});
    deltas.insert({"std", BenchInfo{size}});

    std::vector<std::reference_wrapper<std::vector<int>>> datas{
        std::ref(deltas.at("par").values),
        std::ref(deltas.at("ser").values),
        std::ref(deltas.at("std").values)
    };

    for (int iter = 0; iter < iters; ++iter) {
        for (size_t j = 0; j < size; ++j) {
            datas.front().get()[j] = d(g);
        }

        for (auto& data : datas) {
            std::copy(datas.front().get().begin(), datas.front().get().end(), data.get().begin());
        }

        bench(deltas.at("par"),
                [](int* begin, int* end) {
                    qsort(begin, end, 1000);
                });

        bench(deltas.at("ser"),
                [](int* begin, int* end) {
                    qsort_serial(begin, end);
                });

        bench(deltas.at("std"),
                [](int* begin, int* end) {
                    std::sort(begin, end);
                });

        if (deltas.at("std").values != deltas.at("par").values) {
            auto& v1 = deltas.at("std").values;
            auto& v2 = deltas.at("par").values;
            auto& v3 = deltas.at("ser").values;
            for (size_t i = 0; i < v1.size(); ++i) {
                if (v1[i] != v2[i]) {
                    std::cout << i << " " << v1[i] << " " << v2[i] << " " << v3[i] << std::endl;
                }
            }
            std::cout << "par failed" << std::endl;
        }
        if (deltas.at("std").values != deltas.at("ser").values) {
            std::cout << "seq failed" << std::endl;
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
