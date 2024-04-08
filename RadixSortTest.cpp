#include <iostream>
#include <vector>
#include <array>

template<size_t N>
int getMax(const std::vector<std::array<int, N>>& arr) {
    int max = arr[0][0];
    for (const auto& a : arr) {
        for (int i = 0; i < N; ++i) {
            if (a[i] > max)
                max = a[i];
        }
    }
    return max;
}

template<size_t N>
void countSort(std::vector<std::array<int, N>>& arr, int exp) {
    std::vector<std::array<int, N>> output(arr.size());
    int count[10] = {0};

    for (const auto& a : arr)
        ++count[(a[N - 1] / exp) % 10];

    for (int i = 1; i < 10; ++i)
        count[i] += count[i - 1];

    for (int i = arr.size() - 1; i >= 0; --i) {
        output[count[(arr[i][N - 1] / exp) % 10] - 1] = arr[i];
        --count[(arr[i][N - 1] / exp) % 10];
    }

    for (size_t i = 0; i < arr.size(); ++i)
        arr[i] = output[i];
}

template<size_t N>
void radixSort(std::vector<std::array<int, N>>& arr) {
    int max = getMax(arr);
    for (int exp = 1; max / exp > 0; exp *= 10)
        countSort(arr, exp);
}

int main() {
    // Example usage
    constexpr size_t N = 5; // Size of each array
    std::vector<std::array<int, N>> arr = {
        {170, 45, 75, 90, 802},
        {2, 10, 5, 8, 9},
        {100, 1000, 10000, 1, 2}
    };

    radixSort(arr);

    for (const auto& a : arr) {
        for (int i = 0; i < N; ++i)
            std::cout << a[i] << " ";
        std::cout << std::endl;
    }

    return 0;
}
