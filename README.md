# `triversity`: Diversity Measures on Tripartite Graphs

`triversity` is an R package for the computation of diversity measures on tripartite graphs. First, it implements a parametrized family of such diversity measures which apply on probability distributions. Sometimes called "True Diversity", this family contains famous measures such as the Richness, the Shannon entropy, the Herfindahl-Hirschman index, and the Berger-Parker index. Second, the package allows to apply these measures on probability distributions resulting from random walks between the levels of tripartite graphs. By defining an initial distribution at a given level of the graph and a path to follow between the three levels, the probability of the walker's position within the final level is then computed, thus providing a particular instance of diversity to measure.

### Clone
```
git clone https://github.com/Lamarche-Perrin/triversity
```

### Authors
This package has been developed by researchers of the [Complex Networks](http://www.complexnetworks.fr/) team, within the [Computer Science Laboratory of Paris 6](https://www.lip6.fr/), for the [AlgoDiv](http://algodiv.huma-num.fr/) project, founded by the [French National Agency of Research](http://www.agence-nationale-recherche.fr/) under grant ANR-15-CE38-0001.

List of main collaborators:
- [Robin Lamarche-Perrin](https://www-complexnetworks.lip6.fr/~lamarche/)
- Lionel Tabourier
- Fabien Tarissan
- Raphaël Fournier S'niehotta
- Rémy Cazabet

### License
Copyright © 2017 Robin Lamarche-Perrin (<Robin.Lamarche-Perrin@lip6.fr>)  
`triversity` is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
