#ifndef PAGE_ALLOCATOR_H
#define PAGE_ALLOCATOR_H

#include <vector>

class PageAllocator
{
public:
    enum { default_pgsize = 1024 };

    PageAllocator() = default;
    PageAllocator(const PageAllocator&) = delete;
    PageAllocator& operator=(const PageAllocator&) = delete;
    PageAllocator(PageAllocator&&) = default;
    PageAllocator& operator=(PageAllocator&&) = default;

    void initialize(int pgsize);
    int* allocate(int n);
    void reset();
    int get_pgsize() const;

private:
    struct Page
    {
        int capacity = 0;
        int offset = 0;
        std::vector<int> data;
    };

    std::vector<Page> pages_;
    int pgsize_ = default_pgsize;

    void new_page_();
};

#endif // PAGE_ALLOCATOR_H
