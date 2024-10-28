import matplotlib.pyplot as plt
import numpy as np

x = np.arange(100,dtype=np.float)
print (x.dtype)
m = (x > 10) & (x < 30)
x = np.ma.MaskedArray(x, mask=m)
y = np.ma.MaskedArray(x*x, mask=m)
e = np.ma.MaskedArray(100.*np.ones_like(x), mask=m)
print (x.dtype,y.dtype,e.dtype)
plt.plot(x,y)
plt.errorbar(x,y,yerr=e)
plt.show()
