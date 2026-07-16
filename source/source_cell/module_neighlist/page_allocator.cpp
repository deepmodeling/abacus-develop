#include "page_allocator.h"
#include "source_base/tool_quit.h"

#include <stdexcept>
#include <utility>

PageAllocator::PageAllocator() : pgsize_(default_pgsize)
{
    new_page_();
}

PageAllocator::PageAllocator(int pgsize) : pgsize_(pgsize)
{
    if (pgsize_ <= 0)
    {
        throw std::invalid_argument("PageAllocator page size must be positive.");
    }
    new_page_();
}

PageAllocator::~PageAllocator() = default;

int* PageAllocator::allocate(int n)
{
    if (n <= 0)
    {
        return nullptr;
    }

    if (n > pgsize_)
    {
        ModuleBase::WARNING_QUIT(
            "PageAllocator::allocate",
            "request " + std::to_string(n) + " larger than page size " + std::to_string(pgsize_)
        );
    }

    if (pages_.empty())
    {
        new_page_();
    }

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

void PageAllocator::reset()
{
    pages_.resize(1);
    pages_[0].offset = 0;
}

int PageAllocator::get_pgsize() const
{
    return pgsize_;
}

void PageAllocator::new_page_()
{
    if (pgsize_ <= 0)
    {
        throw std::invalid_argument("PageAllocator page size must be positive.");
    }
    Page p;
    p.capacity = pgsize_;
    p.offset = 0;
    p.data.resize(pgsize_);
    pages_.push_back(std::move(p));
}
