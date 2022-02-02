# plaisir

**p**retty(,) **l**ame **a**nd **i**nefficient **s**olver for **i**nventory **r**outing

This repository contains the source code of plaisir,
a solver for the Inventory Routing Problem,
used for the
[IRP Track of the 12th DIMACS Implementation Challenge](http://dimacs.rutgers.edu/programs/challenge/vrp/irp/).

## Important license notice

plaisir itself uses the MIT License, see the LICENSE file.

However, plaisir includes an adapted version of the
[LKH TSP solver](http://webhotel4.ruc.dk/~keld/research/LKH/)
as third-party code (included in the directory `lkh-sys/lkh/src`).
The code is distributed for academic and non-commercial use only.
The author, Keld Helsgaun, reserves all rights to the code.

## Dependencies:

* Rust 1.58.1
* gcc (or any ANSI C Compiler)
* [Gurobi 9.5](https://www.gurobi.com/)
