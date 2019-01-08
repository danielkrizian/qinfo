
# Functions from [information theory](https://en.wikipedia.org/wiki/Information_theory)

Among others, mutual information between continuous variables.
Discretization of continuous space is implemented with [Darbellay-Vajda adaptive partitioning](http://bsp.teithe.gr/members/downloads/mutin/references/00761290.pdf) algorithm.

Area of application: [Feature selection](https://en.wikipedia.org/wiki/Feature_selection)

```
q)x:copula[10000;.8]  / generate 10k observations of correlated variables
q).info.mut.f x       / mutual information. Better than correlation: detects any kind of (non-linear) relationship
0.4956356
```
