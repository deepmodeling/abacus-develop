#include "mymath.h"

#include "timer.h"

namespace ModuleBase
{

void heapAjust(double *r, int *ind, int s, int m)
{
    int j = 0, ic = 0;
    double rc = 0.0;
    rc = r[s];
    ic = ind[s];

    for (j = 2 * s; j <= m; j *= 2)
    {
        if (j < m && (r[j] < r[j + 1]))
            j++;

        if (!(rc < r[j]))
            break;

        r[s] = r[j];

        ind[s] = ind[j];

        s = j;
    }

    r[s] = rc;

    ind[s] = ic;
    return;
}

void heapsort(const int n, double *r, int *ind)
{
    ModuleBase::timer::tick("mymath", "heapsort");
    int i = 0, ic = 0
    double rc = 0.0;

    if (ind[0] == 0)
    {
        for (i = 0; i < n; i++)
        {
            ind[i] = i;
        }
    }

    for (i = n / 2; i >= 0; i--)
    {
        heapAjust(r, ind, i, n - 1);
    }

    for (i = n - 1; i > 0; i--)
    {
        rc = r[0];
        r[0] = r[i];
        r[i] = rc;
        ic = ind[0];
        ind[0] = ind[i];
        ind[i] = ic;
        heapAjust(r, ind, 0, i - 1);
    }
    ModuleBase::timer::tick("mymath", "heapsort");
    return;
}

/*--------------------------------------------------------------------
      subroutine  hpsort(n,ra,ind)
c---------------------------------------------------------------------
c sort an array ra(1:n) into ascending order using heapsort algorithm.
c n is input, ra is replaced on output by its sorted rearrangement.
c create an index table (ind) by making an exchange in the index array
c whenever an exchange is made on the sorted data array (ra).
c in case of equal values in the data array (ra) the values in the
c index array (ind) are used to order the entries.
c if on input ind(1)  = 0 then indices are initialized in the routine,
c if on input ind(1) != 0 then indices are assumed to have been
c    initialized before entering the routine and these
c    indices are carried around during the sorting process
c
c no work space needed !
c free us from machine-dependent sorting-routines !
c
c adapted from Numerical Recipes pg. 329 (new edition)
*********************************************************************/

// from hpsort.f90
void hpsort(int n, double *ra, int *ind)
{
    int i = 0, ir = 0, j = 0, k = 0, iind = 0;
    double rra = 0.0;

    if (ind[0] == 0)
    {
        for (i = 1; i <= n; i++)
            ind[i - 1] = i;
    }

    if (n < 2)
        return; // nothing to order

    k = n / 2;

    ir = n - 1;

    while (true)
    {
        if (k > 0) // still in hiring phase
        {
            k = k - 1;
            rra = ra[k];
            iind = ind[k];
        }
        else // in retirement-promotion phase.
        {
            rra = ra[ir]; // clear a space at the end of the array
            iind = ind[ir]; //
            ra[ir] = ra[0]; // retire the top of the heap into it
            ind[ir] = ind[0]; //
            ir = ir - 1; // decrease the size of the corporation

            if (ir == 0) // done with the last promotion
            {
                ra[0] = rra; // the least competent worker at all //
                ind[0] = iind; //
                return;
            }
        }

        i = k; // wheter in hiring or promotion phase, we

        j = k + k + 1; // set up to place rra in its proper level

        while (j <= ir)
        {
            if (j < ir)
            {
                if (ra[j] < ra[j + 1]) // compare to better underling
                {
                    j = j + 1;
                }
                else if (ra[j] == ra[j + 1])
                {
                    if (ind[j] < ind[j + 1])
                        j = j + 1;
                }
            }

            if (rra < ra[j]) // demote rra
            {
                ra[i] = ra[j];
                ind[i] = ind[j];
                i = j;
                j = j + j + 1;
            }
            else if (rra == ra[j])
            {
                if (iind < ind[j]) // demote rra
                {
                    ra[i] = ra[j];
                    ind[i] = ind[j];
                    i = j;
                    j = j + j + 1;
                }
                else
                    j = ir + 1; // set j to terminate do-while loop
            }
            else // this is the right place for rra
                j = ir + 1; // set j to terminate do-while loop
        }

        ra[i] = rra;

        ind[i] = iind;
    }
}

} // namespace ModuleBase
