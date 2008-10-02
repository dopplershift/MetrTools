import sys
import numpy as N
import pylab as P
import matplotlib as M
import oban

lam_star = N.linspace(0,10,100)
delta_n = 2.30564
L = 2 * delta_n
lam = lam_star * L
kappa0 = oban.calc_barnes_param(delta_n)

D0 = N.exp(-kappa0 * (N.pi / lam)**2).reshape(-1,1)
gamma = N.array([1.0, 0.4, 0.2]).reshape(1,-1)
D1 = D0 * (1.0 + D0**(gamma - 1.0) - D0**gamma)
D2 = D0 * (1.0 + (1.0 - D0) + (1.0 - D0)**2)

P.plot(lam_star, D0, '--', label=r'$D_0$')
specs = ('c', 'y', 'r')
for ind, g in enumerate(gamma.flatten()):
  P.plot(lam_star, D1[:,ind], specs[ind], label=r'$D_1, \gamma =%.1f$' % g)
P.plot(lam_star, D2, '-.', label=r'$D_2$')

P.xlabel(r'$\lambda/\rm{L}$')
P.ylabel('Filter Response')
P.title('Comparison of Barnes Filter Response Functions')
P.legend(loc = 'lower right')

if len(sys.argv) > 1 and sys.argv[1].startswith('silent'):
  P.savefig('response_funcs.png', dpi=300)
else:
  P.show()
