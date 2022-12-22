# QForm

QForm is a library designed for two purposes:
* To work with the intersection forms of 4-manifolds given a multisection.
* To generate random, valid multisections for further study.

The main interest to the public is likely the former point. This library gives an easy way to compute the intersection of simple closed curves representing classes in $H_1(S_g)$ for some closed oriented surface $S_g$ of genus $g$. It also allows this information to be converted into a matrix representation of the intersection pairing (a symplectic form on $H_1(S_g)$). If these simple closed curves represent, say, a multisection, then QForm also allows for the computation of the intersection form of the resulting 4-manifold.

The latter was helpful in the context of the UGA REU. We were interested in learning what types of bilinear forms can be recognized by means of certain multisections (namely (2,0)-multisections), and so this "random generation" was helpful to experimentally verify that certain forms arise in this way. As such, QForm provides a rudimentary way of generating random symplectic matrices preserving the intersection pairing in order to generate random multisections.
