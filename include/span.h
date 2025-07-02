#ifndef SPAN_H
#define SPAN_H

#include <cstddef>
#include <iterator>
#include <stdexcept>

template<typename T>
class Span {
    public:
        using value_type = T;
        using iterator = T*;
        using const_iterator = const T*;

        Span(T* data, size_t size) : data_(data), size_(size) {}
        Span() : data_(NULL), size_(0){}

        T& operator[](size_t i) const { return data_[i]; }

        T& at(size_t i) const {
            if (i >= size_) throw std::out_of_range("Span index out of range");
            return data_[i];
        }

        T* data() const { return data_; }

        size_t size() const { return size_; }
        bool empty() const { return size_ == 0; }

        auto rbegin() const { return std::reverse_iterator<T*>(data_ + size_); }
        auto rend() const { return std::reverse_iterator<T*>(data_); }

        auto begin() const { return data_; }
        auto end() const { return data_ + size_; }

    private:
        T* data_;
        size_t size_;
};

template <typename T>
bool operator==(const Span<T>& a, const Span<T>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

#endif
