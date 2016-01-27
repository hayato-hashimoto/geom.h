# SYNPOSIS

```C++
#include "geom.h"

geom::fquaternion q(0, 1, 0, 0), r(cos(0.5), sin(0.5), 0, 0), s(cos(0.2), 0, sin(0.2), 0);
geom::fquaternion p = exp(0.2f * log(q));
std::cout << p*p*p*p*p << std::endl;

geom::unit_slerp_interpolater<float> interp(r, s);
for (int i = 0; i <= 10; i++) {
  geom::fquaternion rs = interp(0.1 * i);
  /* conjugation of q by rs */
  geom::fquaternion t  = qrot(rs, q);
  std::cout << t << std::endl;
}

geom::ftransform f = geom::translate(0.2f, 0.2f, 0.2f) * geom::rotateY(0.5f);
std::cout << f * q << std::endl;
```
