#pragma once

#include <bhtsne/tsne.h>
#include <map>

namespace bhtsne
{
    void applyCommandlineOptions(TSNE & tsne, const std::map<std::string, std::string> & options);
}
