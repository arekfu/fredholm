# Copyright 2021 Davide Mancusi <davide.mancusi@cea.fr>
#
# This file is part of fredholm.
#
# fredholm is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# fredholm is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# fredholm.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt

theory = np.loadtxt('theory.dat')
mc_pn = np.loadtxt('fredholm_pn.dat')
mc_el = np.loadtxt('fredholm_el.dat')

fig, ax = plt.subplots(2, 1, sharex=True)

ax[0].errorbar(mc_pn[:, 0], mc_pn[:, 1], mc_pn[:, 2],
               label=r'Monte Carlo, positive-negative sampling', marker='.',
               linestyle='', zorder=1)
ax[0].errorbar(mc_el[:, 0], mc_el[:, 1], mc_el[:, 2],
               label=r'Monte Carlo, linear-exponential sampling', marker='.',
               linestyle='', zorder=2)
ax[0].plot(theory[:, 0], theory[:, 1], label='theory', zorder=3)
ax[0].legend()
ax[0].set_ylabel('density')

diff_pn = (mc_pn[:, 1] - theory[:, 1]) / mc_pn[:, 2]
diff_el = (mc_el[:, 1] - theory[:, 1]) / mc_el[:, 2]
ax[1].plot(mc_pn[:, 0], diff_pn, '.')
ax[1].plot(mc_el[:, 0], diff_el, '.')
ax[1].set_xlabel('position')
ax[1].set_ylabel('Student $t$ variable')
ax[1].set_ylim(-4, 4)

fig.tight_layout()
fig.savefig('fredholm.pdf', bbox_inches='tight')
