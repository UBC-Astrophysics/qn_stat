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

Or you can try the robusttest.py program
~~~
heyl@evi qn_stat-master% python robusttest.py 
Mean Qn value:     1.00005
Mean Std value:    0.991275
Mean MADvalue:     0.993385
Std of Qn:         0.0802589
Std of MAD:        0.115005
Std of std:        0.0709179 (0.0707107)
Efficiency of Qn:  0.883614
Efficiency of MAD: 0.616652
Mean median value: -0.00205408
Mean HL value:     -0.00103553
Mean mean value:   -0.00128172
Std of med:        0.124821 (0.0797885)
Std of hl:         0.103014
Std of mean:       0.100378 (0.1)
Efficiency of med: 0.80418 (0.797885)
Efficiency of hl:  0.974417
heyl@evi qn_stat-master% 
~~~
Original article
========
[Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute deviation. Journal of the American Statistical association, 88(424), 1273-1283.](http://wis.kuleuven.be/stat/robust/papers/publications-1993/rousseeuwcroux-alternativestomedianad-jasa-1993.pdf)

[J.F. Monahan, (1984). Algorithm 616: fast computation of the Hodges-Lehmann location estimator, ACM TOMS 10, 265-270](https://dl.acm.org/doi/abs/10.1145/1271.319414) [Code](http://netlib.org/toms/616.gz)
