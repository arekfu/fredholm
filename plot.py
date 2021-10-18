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
mc = np.loadtxt('fredholm.dat')

fig, ax = plt.subplots(2, 1, sharex=True)

ax[0].errorbar(mc[:, 0], mc[:, 1], mc[:, 2], label='Monte Carlo', marker='.',
             linestyle='')
ax[0].plot(theory[:, 0], theory[:, 1], label='theory')
ax[0].legend()
ax[0].set_ylabel('density')

diff = (mc[:, 1] - theory[:, 1]) / mc[:, 2]
ax[1].plot(mc[:, 0], diff, '.')
ax[1].set_xlabel('position')
ax[1].set_ylabel('Student $t$ variable')
ax[1].set_ylim(-3, 3)

fig.tight_layout()
fig.savefig('fredholm.pdf', bbox_inches='tight')
