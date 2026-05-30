#pragma once

#include <vector>
#include<iostream>
#include <cmath>
#include <cassert>

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
        if (pgsize > 0) new_page();
    }

    PageAllocator(int pgsize_) : pgsize(pgsize_) {
        new_page();
    }

    ~PageAllocator() {
        // no manual delete[]; vectors clean themselves
    }

    PageAllocator(const PageAllocator&) = delete;
    PageAllocator& operator=(const PageAllocator&) = delete;

    // Allow move
    PageAllocator(PageAllocator&&) = default;
    PageAllocator& operator=(PageAllocator&&) = default;

    void new_page() {
        Page p;
        p.capacity = pgsize;
        p.offset = 0;
        p.data.resize(pgsize);
        pages.push_back(std::move(p));
    }

    int* allocate(int n) {
        if (n <= 0) return nullptr;
        // reject requests larger than a single page
        if (n > pgsize) {
            std::cerr << "PageAllocator::allocate error: request " << n << " larger than page size " << pgsize << std::endl;
            return nullptr;
        }
        if (pages.empty()) new_page();
        Page& p = pages.back();
        if (p.offset + n > p.capacity) {
            new_page();
            return allocate(n);
        }
        int* ptr = p.data.data() + p.offset;
        p.offset += n;
        return ptr;
    }

    void reset() {
        pages.resize(1);
        pages[0].offset = 0;
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