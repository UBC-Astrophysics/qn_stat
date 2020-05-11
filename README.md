# qn_stat
C implementation of Rousseeuw and Croux' Qn robust measure of scale

Python Wrapper:  you may have to add the current directory in front of 'robust_c.so' in robust.py: replace 'robust_c.so' with './robust_c.so'
~~~
  In [11]: from robust import qn_calc, hlqest
           import numpy as np
           from robust import qn_calc, hlqest
           import numpy as np
           
    0.002 seconds
 
In [12]: for i in range(10): print(qn_calc(np.random.normal(size=10000)))
         0.127 seconds
         
   0.9988082584224582
   1.0005799701197495
   0.994399639915227
   0.9858328854703164
   1.003603277603626
   1.0119214941582413
   1.0005374885852722
   1.0184855066541487
   0.9978532195511246
   0.9998017180779251
~~~

Original article
========
[Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute deviation. Journal of the American Statistical association, 88(424), 1273-1283.](http://wis.kuleuven.be/stat/robust/papers/publications-1993/rousseeuwcroux-alternativestomedianad-jasa-1993.pdf)

[J.F. Monahan, (1984). Algorithm 616: fast computation of the Hodges-Lehmann location estimator, ACM TOMS 10, 265-270](https://dl.acm.org/doi/abs/10.1145/1271.319414) [Code](http://netlib.org/toms/616.gz)
