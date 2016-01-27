# SYNPOSIS

```C++:
#include "geom.h"
...
fquaternion q(0, 1, 0, 0);
fquaternion p = exp(0.2 * log(q));
std::cout << p*p*p*p*p;
```
