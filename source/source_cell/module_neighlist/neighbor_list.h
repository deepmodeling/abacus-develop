#pragma once

#include <vector>
#include<iostream>
#include <cmath>
#include <cassert>
#include <mutex>

class PageAllocator {
public:
    struct Page {
        std::vector<int> data;
        int capacity;
        int offset;
    };

    std::vector<Page> pages;
    int pgsize;

    static constexpr int DEFAULT_PGSIZE = 1024;

    // Default constructor
    PageAllocator() : pgsize(DEFAULT_PGSIZE) {
        if (pgsize > 0) new_page_unlocked();
    }

    PageAllocator(int pgsize_) : pgsize(pgsize_) {
        new_page_unlocked();
    }

    ~PageAllocator() {
        // no manual delete[]; vectors clean themselves
    }

    PageAllocator(const PageAllocator&) = delete;
    PageAllocator& operator=(const PageAllocator&) = delete;

    // Allow move; mutex state is not moved.
    PageAllocator(PageAllocator&& other) noexcept {
        std::lock_guard<std::mutex> lock(other.mutex_);
        pages = std::move(other.pages);
        pgsize = other.pgsize;
    }

    PageAllocator& operator=(PageAllocator&& other) noexcept {
        if (this != &other) {
            std::lock(mutex_, other.mutex_);
            std::lock_guard<std::mutex> lock_this(mutex_, std::adopt_lock);
            std::lock_guard<std::mutex> lock_other(other.mutex_, std::adopt_lock);
            pages = std::move(other.pages);
            pgsize = other.pgsize;
        }
        return *this;
    }

    void new_page() {
        std::lock_guard<std::mutex> lock(mutex_);
        new_page_unlocked();
    }

    int* allocate(int n) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (n <= 0) return nullptr;
        // reject requests larger than a single page
        if (n > pgsize) {
            std::cerr << "PageAllocator::allocate error: request " << n << " larger than page size " << pgsize << std::endl;
            return nullptr;
        }
        if (pages.empty()) new_page_unlocked();
        Page& p = pages.back();
        if (p.offset + n > p.capacity) {
            new_page_unlocked();
            return allocate_unlocked(n);
        }
        int* ptr = p.data.data() + p.offset;
        p.offset += n;
        return ptr;
    }

    void reset() {
        std::lock_guard<std::mutex> lock(mutex_);
        pages.clear();
        if (pgsize > 0) new_page_unlocked();
    }

private:
    std::mutex mutex_;

    void new_page_unlocked() {
        Page p;
        p.capacity = pgsize;
        p.offset = 0;
        p.data.resize(pgsize);
        pages.push_back(std::move(p));
    }

    int* allocate_unlocked(int n) {
        if (n <= 0) return nullptr;
        Page& p = pages.back();
        assert(p.offset + n <= p.capacity);
        int* ptr = p.data.data() + p.offset;
        p.offset += n;
        return ptr;
    }
};

//////////////////////////////////////////////////////////////
// Neighbor List
//////////////////////////////////////////////////////////////

class NeighborList {
public:
    NeighborList() = default;
    ~NeighborList() = default;

    int nlocal;

    std::vector<int> numneigh;
    std::vector<int*> firstneigh;

    PageAllocator allocator;

    void initialize(int n, int pgsize) {
        nlocal = n;
    allocator = PageAllocator(pgsize);
    // ensure neighbor containers are sized and initialized
    numneigh.assign(n, 0);
    firstneigh.assign(n, nullptr);
    }

    void reset() {
        allocator.reset();
    }
};
