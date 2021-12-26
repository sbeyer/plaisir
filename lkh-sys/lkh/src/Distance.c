#include "LKH.h"

/*
 * Functions for computing distances (see TSPLIB).
 *
 * The appropriate function is referenced by the function pointer Distance.
 */

int Distance(Node * Na, Node * Nb)
{
    double xd = Na->X - Nb->X, yd = Na->Y - Nb->Y;
    return (int) (sqrt(xd * xd + yd * yd) + 0.5);
}
