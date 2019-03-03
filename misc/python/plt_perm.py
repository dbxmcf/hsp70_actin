#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

ppn=np.array([1,2,4,8])
t=np.array([249,261,285,415])

plt.plot(ppn,t,'r-^')
plt.show()
