#ifndef LIBLKH_H
#define LIBLKH_H

struct NodeCoords {
    double x;
    double y;
};

/*
 * Parameters:
 * dimension = number of nodes
 * coords = coordinates of all nodes
 *
 * Return value:
 * a tour (int array of dimension elements)
 */
int const *run(int dimension, struct NodeCoords const * coords);

#endif
