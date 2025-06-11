#include <dsdp/dsdp5.h>

int main(void) {
    DSDP solver;
    SDPCone cone;
    int info;
    info = DSDPCreate(4, &solver);
    if(info) exit(1);
    info = DSDPSetDualObjective(solver, 1, 2.0);
    if(info) exit(1);
    info = DSDPSetDualObjective(solver, 2, 0.0);
    if(info) exit(1);
    info = DSDPSetDualObjective(solver, 3, 2.0);
    if(info) exit(1);
    info = DSDPSetDualObjective(solver, 4, -1.0);
    if(info) exit(1);
  
    // creating cone with all constraints and C matrix
    info = DSDPCreateSDPCone(solver, 1, &cone);
    if(info) exit(1);
  
    // setting the C matrix
    double val0[2] = {-2, 1}, val1[2] = {1, 1}, val2[2] = {1, -1};
    int ind1[2] = {0, 2}, ind2[2] = {0, 5}, ind3[2] = {0, 9}, ind4[2] = {0, 14}, ind5[2] = {0, 20};
    info = SDPConeSetASparseVecMat(cone, 0, 0, 6, 
      1.0, 0,
      ind1, val0, 2);
    if(info) exit(1);
    info = SDPConeViewDataMatrix(cone, 0, 0);
    if(info) exit(1);
    printf("\n");
    info = SDPConeSetASparseVecMat(cone, 0, 1, 6,
      1.0, 0,
      ind2, val1, 2);
    if(info) exit(1);
    info = SDPConeViewDataMatrix(cone, 0, 1);
    if(info) exit(1);
    printf("\n");
    info = SDPConeSetASparseVecMat(cone, 0, 2, 6,
      1.0, 0,
      ind3, val2, 2);
    if(info) exit(1);
    info = SDPConeViewDataMatrix(cone, 0, 2);
    if(info) exit(1);
    printf("\n");
    info = SDPConeSetASparseVecMat(cone, 0, 3, 6,
      1.0, 0,
      ind4, val1, 2);
    if(info) exit(1);
    info = SDPConeViewDataMatrix(cone, 0, 3);
    if(info) exit(1);
    printf("\n");
    info = SDPConeSetASparseVecMat(cone, 0, 4, 6,
      1.0, 0,
      ind5, val2, 2);
    if(info) exit(1);
    info = SDPConeViewDataMatrix(cone, 0, 4);
    if(info) exit(1);
    printf("\n");
    info = DSDPSetGapTolerance(solver, 1e-7);
    if(info) exit(1);
    DSDPSetStandardMonitor(solver, 5);
    info = DSDPSetup(solver);
    if(info) exit(1);
    info = DSDPSolve(solver);
    if(info) exit(1);
    info = DSDPComputeX(solver);
    if(info) exit(1);
    int len;
    double * result = (double*) malloc(30*sizeof(double));
    info = SDPConeGetXArray(cone, 0, &result, &len);
    if(info) exit(1);
    for(int i = 0; i < len; i++)
      printf("%lf ", result[i]);
    double x;
    info = DSDPGetDObjective(solver, &x);
    if(info) exit(1);
    // x = -1/(k-1)
    printf("\nobjective: %lf\n", x);
    DSDPTerminationReason reason;
    info = DSDPStopReason(solver,&reason); 
    DSDPSolutionType pdfeasible;
    DSDPGetSolutionType(solver, &pdfeasible);
    if(pdfeasible == DSDP_PDFEASIBLE)
      printf("feasible\n");
    if(reason == DSDP_CONVERGED)
      printf("converged\n");
    return 0;
}