#ifndef ParamFitHeadings
#define ParamFitHeadings

#include "ParamFitDS.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"

int str4cmp(char* A, char* B);

int strendswith(char* A, char* B);

int TypeCompare(char* T1a, char* T1b, char* T1c, char* T1d, char* T2a,
                char* T2b, char* T2c, char* T2d, int order, int strict);

dmat CompExclMatrix(prmtop *tp, imat *partnb);

imat BuildMatchMatrix(mmsys *myconf, nmroper *myop);

double EvalNMROperation(nmroper *myop, double r);

void FitParams(prmset *mp, trajcon *tj);

#endif
