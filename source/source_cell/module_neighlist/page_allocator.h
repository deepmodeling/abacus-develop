#ifndef PAGE_ALLOCATOR_H
#define PAGE_ALLOCATOR_H

#include <vector>
#include "source_base/tool_quit.h"

class PageAllocator
{
public:
    enum { default_pgsize = 1024 };

    PageAllocator() : pgsize_(default_pgsize)
    {
        new_page_();
    }

    explicit PageAllocator(int pgsize) : pgsize_(pgsize)
    {
        new_page_();
    }

    ~PageAllocator() = default;

    PageAllocator(const PageAllocator&) = delete;
    PageAllocator& operator=(const PageAllocator&) = delete;
    PageAllocator(PageAllocator&&) = default;
    PageAllocator& operator=(PageAllocator&&) = default;

    int* allocate(int n)
    {
        if (n <= 0) return nullptr;

        if (n > pgsize_)
        {
            ModuleBase::WARNING_QUIT(
                "PageAllocator::allocate",
                "request " + std::to_string(n) + " larger than page size " + std::to_string(pgsize_)
            );
        }

        if (pages_.empty()) new_page_();

        Page& p = pages_.back();

        if (p.offset + n > p.capacity)
        {
            new_page_();
            return allocate(n);
        }

        int* ptr = p.data.data() + p.offset;
        p.offset += n;

        return ptr;
    }

    void reset()
    {
        pages_.resize(1);
        pages_[0].offset = 0;
    }

    int get_pgsize() const { return pgsize_; }

private:
    struct Page
    {
        int capacity = 0;
        int offset = 0;
        std::vector<int> data;
    };

    std::vector<Page> pages_;
    int pgsize_ = 0;

    void new_page_()
    {
        Page p;
        p.capacity = pgsize_;
        p.offset = 0;
        p.data.resize(pgsize_);
        pages_.push_back(std::move(p));
    }
};

#endif // PAGE_ALLOCATOR_H