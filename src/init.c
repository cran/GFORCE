#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    kmeans_dual_solution_primal_min_R
    kmeans_dual_solution_primal_min_nok_R
    FORCE_R
    FORCE_adapt_R
    hclust_R
    hclust_agglomerate_R
    kmeans_pp_R
    scio_column_R

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gamma_alternative_estimator_par_R(void *, void *, void *, void *, void *, void *);
extern void gamma_alternative_estimator_R(void *, void *, void *, void *, void *, void *);
extern void v_measure(void *, void *, void *, void *);
extern void v_measure_par(void *, void *, void *, void *);
extern void kmeans_pp_R(void *, void *, void *, void *,void *, void *, void *, void *); // 8
extern void FORCE_R(void *, void *, void *, void *,void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *,void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *,void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *,void *, void *, void *, void *, void *); //39
extern void FORCE_adapt_R(void *, void *, void *, void *,void *, void *, void *, void *, void *, void *,
                        void *, void *, void *, void *,void *, void *, void *, void *, void *, void *,
                        void *, void *, void *, void *,void *, void *, void *, void *, void *, void *,
                        void *, void *, void *, void *,void *, void *, void *, void *); //38
extern void hclust_agglomerate_R(void *, void *, void *, void *,void *); //5
extern void hclust_R(void *, void *, void *, void *,void *); //5

extern void kmeans_dual_solution_primal_min_R(void *, void *, void *, void *,void *,
                                              void *, void *, void *, void *,void *); //10
extern void kmeans_dual_solution_primal_min_nok_R(void *, void *, void *, void *,void *,void *,void *); //7
extern void scio_column_R(void *, void *, void *, void *,void *,void *,void *); //7

/* .C calls to test hooks */
extern void test_smoothed_gradient_S_base(void *, void *, void *, void *,void *); //5
extern void test_smoothed_gradient_X_base(void *, void *, void *, void *,void *,void *); //6
extern void test_smoothed_gradient(void *, void *,void *,void *, void *, void *,void *,void *); // 8
extern void test_project_C_perpendicular(void *, void *, void *, void *,void *); //5
extern void test_project_C_perpendicular_nok(void *, void *, void *, void *); //4
extern void test_project_E(void *, void *, void *, void *,void *); //5
extern void test_smoothed_objective(void *, void *, void *, void *,void *,
                                    void *, void *,void *); //8
extern void test_clust_to_opt_val(void *, void *, void *, void *,void *); //5

extern void test_daps(void *, void *, void *); //3
extern void test_dsmtd(void *, void *, void *, void *); //4
extern void test_dvexp(void *, void *); //2
extern void test_dsumv(void *, void *, void *); //3
extern void test_dtrace(void *, void *, void *); //3
extern void test_dcsum(void *, void *, void *); //3
extern void test_dxpyez(void *, void *, void *, void *); //4

static const R_CMethodDef CEntries[] = {
    {"gamma_alternative_estimator_par_R", (DL_FUNC) &gamma_alternative_estimator_par_R, 6},
    {"gamma_alternative_estimator_R",     (DL_FUNC) &gamma_alternative_estimator_R,     6},
    {"v_measure",                         (DL_FUNC) &v_measure,                         4},
    {"v_measure_par",                     (DL_FUNC) &v_measure_par,                     4},
    {"kmeans_pp_R",                       (DL_FUNC) &kmeans_pp_R,                       8},
    {"FORCE_R",                           (DL_FUNC) &FORCE_R,                          39},
    {"FORCE_adapt_R",                     (DL_FUNC) &FORCE_adapt_R,                    38},
    {"hclust_agglomerate_R",              (DL_FUNC) &hclust_agglomerate_R,              5},
    {"hclust_R",                          (DL_FUNC) &hclust_R,                          5},
    {"kmeans_dual_solution_primal_min_R", (DL_FUNC) &kmeans_dual_solution_primal_min_R, 10},
    {"kmeans_dual_solution_primal_min_nok_R", (DL_FUNC) &kmeans_dual_solution_primal_min_nok_R, 7},
    {"scio_column_R",                     (DL_FUNC) &scio_column_R,                     7},
    {"test_smoothed_gradient_S_base",     (DL_FUNC) &test_smoothed_gradient_S_base,     5},
    {"test_smoothed_gradient_X_base",     (DL_FUNC) &test_smoothed_gradient_X_base,     6},
    {"test_smoothed_gradient",            (DL_FUNC) &test_smoothed_gradient,            8},
    {"test_project_C_perpendicular",      (DL_FUNC) &test_project_C_perpendicular,      5},
    {"test_project_C_perpendicular_nok",  (DL_FUNC) &test_project_C_perpendicular_nok,  4},
    {"test_project_E",                    (DL_FUNC) &test_project_E,                    5},
    {"test_smoothed_objective",           (DL_FUNC) &test_smoothed_objective,           8},
    {"test_clust_to_opt_val",             (DL_FUNC) &test_clust_to_opt_val,             5},
    {"test_daps",                         (DL_FUNC) &test_daps,                         3},
    {"test_dsmtd",                        (DL_FUNC) &test_dsmtd,                        4},
    {"test_dvexp",                        (DL_FUNC) &test_dvexp,                        2},
    {"test_dsumv",                        (DL_FUNC) &test_dsumv,                        3},
    {"test_dtrace",                       (DL_FUNC) &test_dtrace,                       3},
    {"test_dcsum",                        (DL_FUNC) &test_dcsum,                        3},
    {"test_dxpyez",                       (DL_FUNC) &test_dxpyez,                       4},
    {NULL, NULL, 0}
};

void R_init_GFORCE(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
