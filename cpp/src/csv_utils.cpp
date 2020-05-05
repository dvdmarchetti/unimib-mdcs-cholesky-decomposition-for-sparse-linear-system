#include "csv_utils.h"

std::string join(std::vector<std::string> const& elements, const char *const delimiter) {
    std::ostringstream os;
    auto b = begin(elements), e = end(elements);

    if (b != e) {
        std::copy(b, prev(e), std::ostream_iterator<std::string>(os, delimiter));
        b = prev(e);
    }
    if (b != e) {
        os << *b;
    }

    return os.str();
};
