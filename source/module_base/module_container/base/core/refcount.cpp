#include <base/core/refcount.h>
#include "refcount.h"

namespace container{
namespace base {

counted_base::counted_base() : ref_(1) {}

void counted_base::ref() const {
    ref_.fetch_add(1, std::memory_order_relaxed);
}

bool counted_base::unref() const {
    if (ref_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
        delete this;
        return true;
    }
    return false;
}

int_fast32_t counted_base::ref_count() const {
    return ref_.load(std::memory_order_acquire);
}

bool counted_base::ref_count_is_one() const {
    return ref_count() == 1;
}

} // namespace base
} // namespace container