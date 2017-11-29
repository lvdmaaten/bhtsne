#pragma once

#include <bhtsne/tsne.h>

namespace bhtsne {
    void applyCommandlineOptions(TSNE & tsne, int argc, char * argv[]);
}
