/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_sfun.h"
#include "c1_test.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "test_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c1_debug_family_names[9] = { "A", "A_pseudo", "nargin",
  "nargout", "u", "k", "l", "b", "omega" };

/* Function Declarations */
static void initialize_c1_test(SFc1_testInstanceStruct *chartInstance);
static void initialize_params_c1_test(SFc1_testInstanceStruct *chartInstance);
static void enable_c1_test(SFc1_testInstanceStruct *chartInstance);
static void disable_c1_test(SFc1_testInstanceStruct *chartInstance);
static void c1_update_debugger_state_c1_test(SFc1_testInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_test(SFc1_testInstanceStruct
  *chartInstance);
static void set_sim_state_c1_test(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_st);
static void finalize_c1_test(SFc1_testInstanceStruct *chartInstance);
static void sf_c1_test(SFc1_testInstanceStruct *chartInstance);
static void c1_chartstep_c1_test(SFc1_testInstanceStruct *chartInstance);
static void initSimStructsc1_test(SFc1_testInstanceStruct *chartInstance);
static void registerMessagesc1_test(SFc1_testInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static void c1_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_omega, const char_T *c1_identifier, real_T c1_y[6]);
static void c1_b_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[6]);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static real_T c1_c_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_d_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[24]);
static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static void c1_e_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[24]);
static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static void c1_info_helper(c1_ResolvedFunctionInfo c1_info[122]);
static void c1_b_info_helper(c1_ResolvedFunctionInfo c1_info[122]);
static void c1_eml_scalar_eg(SFc1_testInstanceStruct *chartInstance);
static void c1_realmin(SFc1_testInstanceStruct *chartInstance);
static void c1_eps(SFc1_testInstanceStruct *chartInstance);
static void c1_eml_matlab_zgetrf(SFc1_testInstanceStruct *chartInstance, real_T
  c1_A[16], real_T c1_b_A[16], int32_T c1_ipiv[4], int32_T *c1_info);
static void c1_check_forloop_overflow_error(SFc1_testInstanceStruct
  *chartInstance, boolean_T c1_overflow);
static void c1_eml_xger(SFc1_testInstanceStruct *chartInstance, int32_T c1_m,
  int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0, int32_T c1_iy0, real_T c1_A[16],
  int32_T c1_ia0, real_T c1_b_A[16]);
static void c1_eml_warning(SFc1_testInstanceStruct *chartInstance);
static void c1_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A[16],
  real_T c1_B[24], real_T c1_b_B[24]);
static void c1_below_threshold(SFc1_testInstanceStruct *chartInstance);
static void c1_b_eml_scalar_eg(SFc1_testInstanceStruct *chartInstance);
static void c1_b_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A
  [16], real_T c1_B[24], real_T c1_b_B[24]);
static void c1_c_eml_scalar_eg(SFc1_testInstanceStruct *chartInstance);
static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_f_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_g_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_test, const char_T *c1_identifier);
static uint8_T c1_h_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_eml_matlab_zgetrf(SFc1_testInstanceStruct *chartInstance,
  real_T c1_A[16], int32_T c1_ipiv[4], int32_T *c1_info);
static void c1_b_eml_xger(SFc1_testInstanceStruct *chartInstance, int32_T c1_m,
  int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0, int32_T c1_iy0, real_T c1_A[16],
  int32_T c1_ia0);
static void c1_c_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A
  [16], real_T c1_B[24]);
static void c1_d_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A
  [16], real_T c1_B[24]);
static void init_dsm_address_info(SFc1_testInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c1_test(SFc1_testInstanceStruct *chartInstance)
{
  chartInstance->c1_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c1_is_active_c1_test = 0U;
}

static void initialize_params_c1_test(SFc1_testInstanceStruct *chartInstance)
{
}

static void enable_c1_test(SFc1_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c1_test(SFc1_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c1_update_debugger_state_c1_test(SFc1_testInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c1_test(SFc1_testInstanceStruct
  *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  int32_T c1_i0;
  real_T c1_u[6];
  const mxArray *c1_b_y = NULL;
  uint8_T c1_hoistedGlobal;
  uint8_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  real_T (*c1_omega)[6];
  c1_omega = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellarray(2), FALSE);
  for (c1_i0 = 0; c1_i0 < 6; c1_i0++) {
    c1_u[c1_i0] = (*c1_omega)[c1_i0];
  }

  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_hoistedGlobal = chartInstance->c1_is_active_c1_test;
  c1_b_u = c1_hoistedGlobal;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  sf_mex_assign(&c1_st, c1_y, FALSE);
  return c1_st;
}

static void set_sim_state_c1_test(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_st)
{
  const mxArray *c1_u;
  real_T c1_dv0[6];
  int32_T c1_i1;
  real_T (*c1_omega)[6];
  c1_omega = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c1_doneDoubleBufferReInit = TRUE;
  c1_u = sf_mex_dup(c1_st);
  c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u, 0)),
                      "omega", c1_dv0);
  for (c1_i1 = 0; c1_i1 < 6; c1_i1++) {
    (*c1_omega)[c1_i1] = c1_dv0[c1_i1];
  }

  chartInstance->c1_is_active_c1_test = c1_g_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_u, 1)), "is_active_c1_test");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_test(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_test(SFc1_testInstanceStruct *chartInstance)
{
}

static void sf_c1_test(SFc1_testInstanceStruct *chartInstance)
{
  int32_T c1_i2;
  int32_T c1_i3;
  real_T *c1_k;
  real_T *c1_l;
  real_T *c1_b;
  real_T (*c1_omega)[6];
  real_T (*c1_u)[4];
  c1_b = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c1_l = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c1_omega = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  for (c1_i2 = 0; c1_i2 < 4; c1_i2++) {
    _SFD_DATA_RANGE_CHECK((*c1_u)[c1_i2], 0U);
  }

  for (c1_i3 = 0; c1_i3 < 6; c1_i3++) {
    _SFD_DATA_RANGE_CHECK((*c1_omega)[c1_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c1_k, 2U);
  _SFD_DATA_RANGE_CHECK(*c1_l, 3U);
  _SFD_DATA_RANGE_CHECK(*c1_b, 4U);
  chartInstance->c1_sfEvent = CALL_EVENT;
  c1_chartstep_c1_test(chartInstance);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_testMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c1_chartstep_c1_test(SFc1_testInstanceStruct *chartInstance)
{
  real_T c1_hoistedGlobal;
  real_T c1_b_hoistedGlobal;
  real_T c1_c_hoistedGlobal;
  int32_T c1_i4;
  real_T c1_u[4];
  real_T c1_k;
  real_T c1_l;
  real_T c1_b;
  uint32_T c1_debug_family_var_map[9];
  real_T c1_A[24];
  real_T c1_A_pseudo[24];
  real_T c1_nargin = 4.0;
  real_T c1_nargout = 1.0;
  real_T c1_omega[6];
  real_T c1_b_b;
  real_T c1_y;
  real_T c1_a;
  real_T c1_c_b;
  real_T c1_b_y;
  real_T c1_d_b;
  real_T c1_c_y;
  real_T c1_b_a;
  real_T c1_e_b;
  real_T c1_d_y;
  real_T c1_c_a;
  real_T c1_f_b;
  real_T c1_e_y;
  real_T c1_g_b;
  real_T c1_f_y;
  real_T c1_d_a;
  real_T c1_h_b;
  real_T c1_g_y;
  real_T c1_i_b;
  real_T c1_h_y;
  real_T c1_e_a;
  real_T c1_j_b;
  real_T c1_i_y;
  real_T c1_f_a;
  real_T c1_k_b;
  real_T c1_j_y;
  real_T c1_l_b;
  real_T c1_k_y;
  real_T c1_g_a;
  real_T c1_m_b;
  real_T c1_l_y;
  real_T c1_n_b;
  real_T c1_m_y;
  real_T c1_h_a;
  real_T c1_o_b;
  real_T c1_n_y;
  real_T c1_p_b;
  real_T c1_o_y;
  real_T c1_i_a;
  real_T c1_q_b;
  real_T c1_p_y;
  real_T c1_r_b;
  real_T c1_q_y;
  real_T c1_j_a;
  real_T c1_s_b;
  real_T c1_r_y;
  int32_T c1_i5;
  real_T c1_k_a[24];
  int32_T c1_i6;
  int32_T c1_i7;
  int32_T c1_i8;
  int32_T c1_i9;
  real_T c1_t_b[24];
  int32_T c1_i10;
  int32_T c1_i11;
  int32_T c1_i12;
  int32_T c1_i13;
  real_T c1_s_y[16];
  int32_T c1_i14;
  int32_T c1_i15;
  int32_T c1_i16;
  int32_T c1_i17;
  int32_T c1_i18;
  int32_T c1_i19;
  int32_T c1_i20;
  int32_T c1_i21;
  int32_T c1_i22;
  int32_T c1_i23;
  real_T c1_b_A[16];
  int32_T c1_i24;
  int32_T c1_i25;
  int32_T c1_i26;
  int32_T c1_i27;
  int32_T c1_info;
  int32_T c1_ipiv[4];
  int32_T c1_b_info;
  int32_T c1_c_info;
  int32_T c1_d_info;
  int32_T c1_i;
  int32_T c1_b_i;
  int32_T c1_ip;
  int32_T c1_j;
  int32_T c1_b_j;
  real_T c1_temp;
  int32_T c1_i28;
  real_T c1_c_A[16];
  int32_T c1_i29;
  real_T c1_d_A[16];
  int32_T c1_i30;
  int32_T c1_i31;
  int32_T c1_i32;
  int32_T c1_i33;
  int32_T c1_i34;
  int32_T c1_i35;
  real_T c1_u_b[4];
  int32_T c1_i36;
  int32_T c1_i37;
  int32_T c1_i38;
  real_T c1_C[6];
  int32_T c1_i39;
  int32_T c1_i40;
  int32_T c1_i41;
  int32_T c1_i42;
  int32_T c1_i43;
  int32_T c1_i44;
  int32_T c1_i45;
  real_T (*c1_b_omega)[6];
  real_T *c1_v_b;
  real_T *c1_b_l;
  real_T *c1_b_k;
  real_T (*c1_b_u)[4];
  c1_v_b = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c1_b_l = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_b_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c1_b_omega = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c1_b_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
  c1_hoistedGlobal = *c1_b_k;
  c1_b_hoistedGlobal = *c1_b_l;
  c1_c_hoistedGlobal = *c1_v_b;
  for (c1_i4 = 0; c1_i4 < 4; c1_i4++) {
    c1_u[c1_i4] = (*c1_b_u)[c1_i4];
  }

  c1_k = c1_hoistedGlobal;
  c1_l = c1_b_hoistedGlobal;
  c1_b = c1_c_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 9U, 9U, c1_debug_family_names,
    c1_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_A, 0U, c1_e_sf_marshallOut,
    c1_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_A_pseudo, 1U, c1_d_sf_marshallOut,
    c1_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargin, 2U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c1_nargout, 3U, c1_b_sf_marshallOut,
    c1_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c1_u, 4U, c1_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_k, 5U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_l, 6U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c1_b, 7U, c1_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c1_omega, 8U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 4);
  c1_b_b = c1_k;
  c1_y = -0.5 * c1_b_b;
  c1_a = c1_y;
  c1_c_b = c1_l;
  c1_b_y = c1_a * c1_c_b;
  c1_d_b = c1_k;
  c1_c_y = 0.5 * c1_d_b;
  c1_b_a = c1_c_y;
  c1_e_b = c1_l;
  c1_d_y = c1_b_a * c1_e_b;
  c1_c_a = c1_k;
  c1_f_b = c1_l;
  c1_e_y = c1_c_a * c1_f_b;
  c1_g_b = c1_k;
  c1_f_y = 0.5 * c1_g_b;
  c1_d_a = c1_f_y;
  c1_h_b = c1_l;
  c1_g_y = c1_d_a * c1_h_b;
  c1_i_b = c1_k;
  c1_h_y = -0.5 * c1_i_b;
  c1_e_a = c1_h_y;
  c1_j_b = c1_l;
  c1_i_y = c1_e_a * c1_j_b;
  c1_f_a = -c1_k;
  c1_k_b = c1_l;
  c1_j_y = c1_f_a * c1_k_b;
  c1_l_b = c1_k;
  c1_k_y = 0.8660254037844386 * c1_l_b;
  c1_g_a = c1_k_y;
  c1_m_b = c1_l;
  c1_l_y = c1_g_a * c1_m_b;
  c1_n_b = c1_k;
  c1_m_y = 0.8660254037844386 * c1_n_b;
  c1_h_a = c1_m_y;
  c1_o_b = c1_l;
  c1_n_y = c1_h_a * c1_o_b;
  c1_p_b = c1_k;
  c1_o_y = -0.8660254037844386 * c1_p_b;
  c1_i_a = c1_o_y;
  c1_q_b = c1_l;
  c1_p_y = c1_i_a * c1_q_b;
  c1_r_b = c1_k;
  c1_q_y = -0.8660254037844386 * c1_r_b;
  c1_j_a = c1_q_y;
  c1_s_b = c1_l;
  c1_r_y = c1_j_a * c1_s_b;
  c1_A[0] = c1_k;
  c1_A[4] = c1_k;
  c1_A[8] = c1_k;
  c1_A[12] = c1_k;
  c1_A[16] = c1_k;
  c1_A[20] = c1_k;
  c1_A[1] = c1_b_y;
  c1_A[5] = c1_d_y;
  c1_A[9] = c1_e_y;
  c1_A[13] = c1_g_y;
  c1_A[17] = c1_i_y;
  c1_A[21] = c1_j_y;
  c1_A[2] = c1_l_y;
  c1_A[6] = c1_n_y;
  c1_A[10] = 0.0;
  c1_A[14] = c1_p_y;
  c1_A[18] = c1_r_y;
  c1_A[22] = 0.0;
  c1_A[3] = c1_b;
  c1_A[7] = -c1_b;
  c1_A[11] = c1_b;
  c1_A[15] = -c1_b;
  c1_A[19] = c1_b;
  c1_A[23] = -c1_b;
  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 9);
  for (c1_i5 = 0; c1_i5 < 24; c1_i5++) {
    c1_k_a[c1_i5] = c1_A[c1_i5];
  }

  c1_i6 = 0;
  for (c1_i7 = 0; c1_i7 < 4; c1_i7++) {
    c1_i8 = 0;
    for (c1_i9 = 0; c1_i9 < 6; c1_i9++) {
      c1_t_b[c1_i9 + c1_i6] = c1_A[c1_i8 + c1_i7];
      c1_i8 += 4;
    }

    c1_i6 += 6;
  }

  c1_eml_scalar_eg(chartInstance);
  c1_eml_scalar_eg(chartInstance);
  for (c1_i10 = 0; c1_i10 < 4; c1_i10++) {
    c1_i11 = 0;
    c1_i12 = 0;
    for (c1_i13 = 0; c1_i13 < 4; c1_i13++) {
      c1_s_y[c1_i11 + c1_i10] = 0.0;
      c1_i14 = 0;
      for (c1_i15 = 0; c1_i15 < 6; c1_i15++) {
        c1_s_y[c1_i11 + c1_i10] += c1_k_a[c1_i14 + c1_i10] * c1_t_b[c1_i15 +
          c1_i12];
        c1_i14 += 4;
      }

      c1_i11 += 4;
      c1_i12 += 6;
    }
  }

  c1_i16 = 0;
  for (c1_i17 = 0; c1_i17 < 4; c1_i17++) {
    c1_i18 = 0;
    for (c1_i19 = 0; c1_i19 < 6; c1_i19++) {
      c1_t_b[c1_i19 + c1_i16] = c1_A[c1_i18 + c1_i17];
      c1_i18 += 4;
    }

    c1_i16 += 6;
  }

  c1_i20 = 0;
  for (c1_i21 = 0; c1_i21 < 4; c1_i21++) {
    c1_i22 = 0;
    for (c1_i23 = 0; c1_i23 < 4; c1_i23++) {
      c1_b_A[c1_i23 + c1_i20] = c1_s_y[c1_i22 + c1_i21];
      c1_i22 += 4;
    }

    c1_i20 += 4;
  }

  c1_i24 = 0;
  for (c1_i25 = 0; c1_i25 < 6; c1_i25++) {
    c1_i26 = 0;
    for (c1_i27 = 0; c1_i27 < 4; c1_i27++) {
      c1_k_a[c1_i27 + c1_i24] = c1_t_b[c1_i26 + c1_i25];
      c1_i26 += 6;
    }

    c1_i24 += 4;
  }

  c1_b_eml_matlab_zgetrf(chartInstance, c1_b_A, c1_ipiv, &c1_info);
  c1_b_info = c1_info;
  c1_c_info = c1_b_info;
  c1_d_info = c1_c_info;
  if (c1_d_info > 0) {
    c1_eml_warning(chartInstance);
  }

  for (c1_i = 1; c1_i < 5; c1_i++) {
    c1_b_i = c1_i;
    if (c1_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_i), 1, 4, 1, 0) - 1] != c1_b_i) {
      c1_ip = c1_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 4, 1, 0) - 1];
      for (c1_j = 1; c1_j < 7; c1_j++) {
        c1_b_j = c1_j;
        c1_temp = c1_k_a[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 4, 1, 0) +
                          ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 6, 2, 0) - 1) << 2)) - 1];
        c1_k_a[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_b_i), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK(
                   "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 6, 2,
                   0) - 1) << 2)) - 1] = c1_k_a[(_SFD_EML_ARRAY_BOUNDS_CHECK("",
          (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_ip), 1, 4, 1, 0) +
          ((_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_j), 1, 6, 2, 0) - 1) << 2)) - 1];
        c1_k_a[(_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
                  (real_T)c1_ip), 1, 4, 1, 0) + ((_SFD_EML_ARRAY_BOUNDS_CHECK("",
                   (int32_T)_SFD_INTEGER_CHECK("", (real_T)c1_b_j), 1, 6, 2, 0)
                  - 1) << 2)) - 1] = c1_temp;
      }
    }
  }

  for (c1_i28 = 0; c1_i28 < 16; c1_i28++) {
    c1_c_A[c1_i28] = c1_b_A[c1_i28];
  }

  c1_c_eml_xtrsm(chartInstance, c1_c_A, c1_k_a);
  for (c1_i29 = 0; c1_i29 < 16; c1_i29++) {
    c1_d_A[c1_i29] = c1_b_A[c1_i29];
  }

  c1_d_eml_xtrsm(chartInstance, c1_d_A, c1_k_a);
  c1_i30 = 0;
  for (c1_i31 = 0; c1_i31 < 4; c1_i31++) {
    c1_i32 = 0;
    for (c1_i33 = 0; c1_i33 < 6; c1_i33++) {
      c1_A_pseudo[c1_i33 + c1_i30] = c1_k_a[c1_i32 + c1_i31];
      c1_i32 += 4;
    }

    c1_i30 += 6;
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, 11);
  for (c1_i34 = 0; c1_i34 < 24; c1_i34++) {
    c1_t_b[c1_i34] = c1_A_pseudo[c1_i34];
  }

  for (c1_i35 = 0; c1_i35 < 4; c1_i35++) {
    c1_u_b[c1_i35] = c1_u[c1_i35];
  }

  c1_c_eml_scalar_eg(chartInstance);
  c1_c_eml_scalar_eg(chartInstance);
  for (c1_i36 = 0; c1_i36 < 6; c1_i36++) {
    c1_omega[c1_i36] = 0.0;
  }

  for (c1_i37 = 0; c1_i37 < 6; c1_i37++) {
    c1_omega[c1_i37] = 0.0;
  }

  for (c1_i38 = 0; c1_i38 < 6; c1_i38++) {
    c1_C[c1_i38] = c1_omega[c1_i38];
  }

  for (c1_i39 = 0; c1_i39 < 6; c1_i39++) {
    c1_omega[c1_i39] = c1_C[c1_i39];
  }

  for (c1_i40 = 0; c1_i40 < 6; c1_i40++) {
    c1_C[c1_i40] = c1_omega[c1_i40];
  }

  for (c1_i41 = 0; c1_i41 < 6; c1_i41++) {
    c1_omega[c1_i41] = c1_C[c1_i41];
  }

  for (c1_i42 = 0; c1_i42 < 6; c1_i42++) {
    c1_omega[c1_i42] = 0.0;
    c1_i43 = 0;
    for (c1_i44 = 0; c1_i44 < 4; c1_i44++) {
      c1_omega[c1_i42] += c1_t_b[c1_i43 + c1_i42] * c1_u_b[c1_i44];
      c1_i43 += 6;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c1_sfEvent, -11);
  _SFD_SYMBOL_SCOPE_POP();
  for (c1_i45 = 0; c1_i45 < 6; c1_i45++) {
    (*c1_b_omega)[c1_i45] = c1_omega[c1_i45];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c1_sfEvent);
}

static void initSimStructsc1_test(SFc1_testInstanceStruct *chartInstance)
{
}

static void registerMessagesc1_test(SFc1_testInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber)
{
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i46;
  real_T c1_b_inData[6];
  int32_T c1_i47;
  real_T c1_u[6];
  const mxArray *c1_y = NULL;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i46 = 0; c1_i46 < 6; c1_i46++) {
    c1_b_inData[c1_i46] = (*(real_T (*)[6])c1_inData)[c1_i46];
  }

  for (c1_i47 = 0; c1_i47 < 6; c1_i47++) {
    c1_u[c1_i47] = c1_b_inData[c1_i47];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static void c1_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_omega, const char_T *c1_identifier, real_T c1_y[6])
{
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_omega), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_omega);
}

static void c1_b_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[6])
{
  real_T c1_dv1[6];
  int32_T c1_i48;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv1, 1, 0, 0U, 1, 0U, 1, 6);
  for (c1_i48 = 0; c1_i48 < 6; c1_i48++) {
    c1_y[c1_i48] = c1_dv1[c1_i48];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_omega;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[6];
  int32_T c1_i49;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_omega = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_omega), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_omega);
  for (c1_i49 = 0; c1_i49 < 6; c1_i49++) {
    (*(real_T (*)[6])c1_outData)[c1_i49] = c1_y[c1_i49];
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static const mxArray *c1_c_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i50;
  real_T c1_b_inData[4];
  int32_T c1_i51;
  real_T c1_u[4];
  const mxArray *c1_y = NULL;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  for (c1_i50 = 0; c1_i50 < 4; c1_i50++) {
    c1_b_inData[c1_i50] = (*(real_T (*)[4])c1_inData)[c1_i50];
  }

  for (c1_i51 = 0; c1_i51 < 4; c1_i51++) {
    c1_u[c1_i51] = c1_b_inData[c1_i51];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static real_T c1_c_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d0;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d0, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_nargout;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_nargout = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_nargout), &c1_thisId);
  sf_mex_destroy(&c1_nargout);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_d_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i52;
  int32_T c1_i53;
  int32_T c1_i54;
  real_T c1_b_inData[24];
  int32_T c1_i55;
  int32_T c1_i56;
  int32_T c1_i57;
  real_T c1_u[24];
  const mxArray *c1_y = NULL;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i52 = 0;
  for (c1_i53 = 0; c1_i53 < 4; c1_i53++) {
    for (c1_i54 = 0; c1_i54 < 6; c1_i54++) {
      c1_b_inData[c1_i54 + c1_i52] = (*(real_T (*)[24])c1_inData)[c1_i54 +
        c1_i52];
    }

    c1_i52 += 6;
  }

  c1_i55 = 0;
  for (c1_i56 = 0; c1_i56 < 4; c1_i56++) {
    for (c1_i57 = 0; c1_i57 < 6; c1_i57++) {
      c1_u[c1_i57 + c1_i55] = c1_b_inData[c1_i57 + c1_i55];
    }

    c1_i55 += 6;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 6, 4), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static void c1_d_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[24])
{
  real_T c1_dv2[24];
  int32_T c1_i58;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv2, 1, 0, 0U, 1, 0U, 2, 6, 4);
  for (c1_i58 = 0; c1_i58 < 24; c1_i58++) {
    c1_y[c1_i58] = c1_dv2[c1_i58];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_A_pseudo;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[24];
  int32_T c1_i59;
  int32_T c1_i60;
  int32_T c1_i61;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_A_pseudo = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_A_pseudo), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_A_pseudo);
  c1_i59 = 0;
  for (c1_i60 = 0; c1_i60 < 4; c1_i60++) {
    for (c1_i61 = 0; c1_i61 < 6; c1_i61++) {
      (*(real_T (*)[24])c1_outData)[c1_i61 + c1_i59] = c1_y[c1_i61 + c1_i59];
    }

    c1_i59 += 6;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

static const mxArray *c1_e_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_i62;
  int32_T c1_i63;
  int32_T c1_i64;
  real_T c1_b_inData[24];
  int32_T c1_i65;
  int32_T c1_i66;
  int32_T c1_i67;
  real_T c1_u[24];
  const mxArray *c1_y = NULL;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_i62 = 0;
  for (c1_i63 = 0; c1_i63 < 6; c1_i63++) {
    for (c1_i64 = 0; c1_i64 < 4; c1_i64++) {
      c1_b_inData[c1_i64 + c1_i62] = (*(real_T (*)[24])c1_inData)[c1_i64 +
        c1_i62];
    }

    c1_i62 += 4;
  }

  c1_i65 = 0;
  for (c1_i66 = 0; c1_i66 < 6; c1_i66++) {
    for (c1_i67 = 0; c1_i67 < 4; c1_i67++) {
      c1_u[c1_i67 + c1_i65] = c1_b_inData[c1_i67 + c1_i65];
    }

    c1_i65 += 4;
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 0, 0U, 1U, 0U, 2, 4, 6), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static void c1_e_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance, const
  mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId, real_T c1_y[24])
{
  real_T c1_dv3[24];
  int32_T c1_i68;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), c1_dv3, 1, 0, 0U, 1, 0U, 2, 4, 6);
  for (c1_i68 = 0; c1_i68 < 24; c1_i68++) {
    c1_y[c1_i68] = c1_dv3[c1_i68];
  }

  sf_mex_destroy(&c1_u);
}

static void c1_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_A;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y[24];
  int32_T c1_i69;
  int32_T c1_i70;
  int32_T c1_i71;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_A = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_A), &c1_thisId, c1_y);
  sf_mex_destroy(&c1_A);
  c1_i69 = 0;
  for (c1_i70 = 0; c1_i70 < 6; c1_i70++) {
    for (c1_i71 = 0; c1_i71 < 4; c1_i71++) {
      (*(real_T (*)[24])c1_outData)[c1_i71 + c1_i69] = c1_y[c1_i71 + c1_i69];
    }

    c1_i69 += 4;
  }

  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_test_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo;
  c1_ResolvedFunctionInfo c1_info[122];
  const mxArray *c1_m0 = NULL;
  int32_T c1_i72;
  c1_ResolvedFunctionInfo *c1_r0;
  c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  c1_info_helper(c1_info);
  c1_b_info_helper(c1_info);
  sf_mex_assign(&c1_m0, sf_mex_createstruct("nameCaptureInfo", 1, 122), FALSE);
  for (c1_i72 = 0; c1_i72 < 122; c1_i72++) {
    c1_r0 = &c1_info[c1_i72];
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c1_r0->context)), "context", "nameCaptureInfo",
                    c1_i72);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c1_r0->name)), "name", "nameCaptureInfo", c1_i72);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c1_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c1_i72);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c1_r0->resolved)), "resolved", "nameCaptureInfo",
                    c1_i72);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c1_i72);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c1_i72);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c1_i72);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c1_i72);
  }

  sf_mex_assign(&c1_nameCaptureInfo, c1_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static void c1_info_helper(c1_ResolvedFunctionInfo c1_info[122])
{
  c1_info[0].context = "";
  c1_info[0].name = "mtimes";
  c1_info[0].dominantType = "double";
  c1_info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[0].fileTimeLo = 1289519692U;
  c1_info[0].fileTimeHi = 0U;
  c1_info[0].mFileTimeLo = 0U;
  c1_info[0].mFileTimeHi = 0U;
  c1_info[1].context = "";
  c1_info[1].name = "sqrt";
  c1_info[1].dominantType = "double";
  c1_info[1].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c1_info[1].fileTimeLo = 1343830386U;
  c1_info[1].fileTimeHi = 0U;
  c1_info[1].mFileTimeLo = 0U;
  c1_info[1].mFileTimeHi = 0U;
  c1_info[2].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c1_info[2].name = "eml_error";
  c1_info[2].dominantType = "char";
  c1_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  c1_info[2].fileTimeLo = 1343830358U;
  c1_info[2].fileTimeHi = 0U;
  c1_info[2].mFileTimeLo = 0U;
  c1_info[2].mFileTimeHi = 0U;
  c1_info[3].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c1_info[3].name = "eml_scalar_sqrt";
  c1_info[3].dominantType = "double";
  c1_info[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  c1_info[3].fileTimeLo = 1286818738U;
  c1_info[3].fileTimeHi = 0U;
  c1_info[3].mFileTimeLo = 0U;
  c1_info[3].mFileTimeHi = 0U;
  c1_info[4].context = "";
  c1_info[4].name = "mrdivide";
  c1_info[4].dominantType = "double";
  c1_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c1_info[4].fileTimeLo = 1357951548U;
  c1_info[4].fileTimeHi = 0U;
  c1_info[4].mFileTimeLo = 1319729966U;
  c1_info[4].mFileTimeHi = 0U;
  c1_info[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c1_info[5].name = "rdivide";
  c1_info[5].dominantType = "double";
  c1_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c1_info[5].fileTimeLo = 1346510388U;
  c1_info[5].fileTimeHi = 0U;
  c1_info[5].mFileTimeLo = 0U;
  c1_info[5].mFileTimeHi = 0U;
  c1_info[6].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c1_info[6].name = "eml_scalexp_compatible";
  c1_info[6].dominantType = "double";
  c1_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  c1_info[6].fileTimeLo = 1286818796U;
  c1_info[6].fileTimeHi = 0U;
  c1_info[6].mFileTimeLo = 0U;
  c1_info[6].mFileTimeHi = 0U;
  c1_info[7].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c1_info[7].name = "eml_div";
  c1_info[7].dominantType = "double";
  c1_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c1_info[7].fileTimeLo = 1313347810U;
  c1_info[7].fileTimeHi = 0U;
  c1_info[7].mFileTimeLo = 0U;
  c1_info[7].mFileTimeHi = 0U;
  c1_info[8].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[8].name = "eml_index_class";
  c1_info[8].dominantType = "";
  c1_info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[8].fileTimeLo = 1323170578U;
  c1_info[8].fileTimeHi = 0U;
  c1_info[8].mFileTimeLo = 0U;
  c1_info[8].mFileTimeHi = 0U;
  c1_info[9].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[9].name = "eml_scalar_eg";
  c1_info[9].dominantType = "double";
  c1_info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[9].fileTimeLo = 1286818796U;
  c1_info[9].fileTimeHi = 0U;
  c1_info[9].mFileTimeLo = 0U;
  c1_info[9].mFileTimeHi = 0U;
  c1_info[10].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[10].name = "eml_xgemm";
  c1_info[10].dominantType = "char";
  c1_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c1_info[10].fileTimeLo = 1299076772U;
  c1_info[10].fileTimeHi = 0U;
  c1_info[10].mFileTimeLo = 0U;
  c1_info[10].mFileTimeHi = 0U;
  c1_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c1_info[11].name = "eml_blas_inline";
  c1_info[11].dominantType = "";
  c1_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c1_info[11].fileTimeLo = 1299076768U;
  c1_info[11].fileTimeHi = 0U;
  c1_info[11].mFileTimeLo = 0U;
  c1_info[11].mFileTimeHi = 0U;
  c1_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold";
  c1_info[12].name = "mtimes";
  c1_info[12].dominantType = "double";
  c1_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[12].fileTimeLo = 1289519692U;
  c1_info[12].fileTimeHi = 0U;
  c1_info[12].mFileTimeLo = 0U;
  c1_info[12].mFileTimeHi = 0U;
  c1_info[13].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c1_info[13].name = "eml_index_class";
  c1_info[13].dominantType = "";
  c1_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[13].fileTimeLo = 1323170578U;
  c1_info[13].fileTimeHi = 0U;
  c1_info[13].mFileTimeLo = 0U;
  c1_info[13].mFileTimeHi = 0U;
  c1_info[14].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c1_info[14].name = "eml_scalar_eg";
  c1_info[14].dominantType = "double";
  c1_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[14].fileTimeLo = 1286818796U;
  c1_info[14].fileTimeHi = 0U;
  c1_info[14].mFileTimeLo = 0U;
  c1_info[14].mFileTimeHi = 0U;
  c1_info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c1_info[15].name = "eml_refblas_xgemm";
  c1_info[15].dominantType = "char";
  c1_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c1_info[15].fileTimeLo = 1299076774U;
  c1_info[15].fileTimeHi = 0U;
  c1_info[15].mFileTimeLo = 0U;
  c1_info[15].mFileTimeHi = 0U;
  c1_info[16].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c1_info[16].name = "mldivide";
  c1_info[16].dominantType = "double";
  c1_info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c1_info[16].fileTimeLo = 1357951548U;
  c1_info[16].fileTimeHi = 0U;
  c1_info[16].mFileTimeLo = 1319729966U;
  c1_info[16].mFileTimeHi = 0U;
  c1_info[17].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c1_info[17].name = "eml_lusolve";
  c1_info[17].dominantType = "double";
  c1_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c1_info[17].fileTimeLo = 1309451196U;
  c1_info[17].fileTimeHi = 0U;
  c1_info[17].mFileTimeLo = 0U;
  c1_info[17].mFileTimeHi = 0U;
  c1_info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c1_info[18].name = "eml_index_class";
  c1_info[18].dominantType = "";
  c1_info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[18].fileTimeLo = 1323170578U;
  c1_info[18].fileTimeHi = 0U;
  c1_info[18].mFileTimeLo = 0U;
  c1_info[18].mFileTimeHi = 0U;
  c1_info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c1_info[19].name = "eml_index_class";
  c1_info[19].dominantType = "";
  c1_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[19].fileTimeLo = 1323170578U;
  c1_info[19].fileTimeHi = 0U;
  c1_info[19].mFileTimeLo = 0U;
  c1_info[19].mFileTimeHi = 0U;
  c1_info[20].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c1_info[20].name = "eml_xgetrf";
  c1_info[20].dominantType = "double";
  c1_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c1_info[20].fileTimeLo = 1286818806U;
  c1_info[20].fileTimeHi = 0U;
  c1_info[20].mFileTimeLo = 0U;
  c1_info[20].mFileTimeHi = 0U;
  c1_info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c1_info[21].name = "eml_lapack_xgetrf";
  c1_info[21].dominantType = "double";
  c1_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c1_info[21].fileTimeLo = 1286818810U;
  c1_info[21].fileTimeHi = 0U;
  c1_info[21].mFileTimeLo = 0U;
  c1_info[21].mFileTimeHi = 0U;
  c1_info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c1_info[22].name = "eml_matlab_zgetrf";
  c1_info[22].dominantType = "double";
  c1_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[22].fileTimeLo = 1302688994U;
  c1_info[22].fileTimeHi = 0U;
  c1_info[22].mFileTimeLo = 0U;
  c1_info[22].mFileTimeHi = 0U;
  c1_info[23].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[23].name = "realmin";
  c1_info[23].dominantType = "char";
  c1_info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m";
  c1_info[23].fileTimeLo = 1307651242U;
  c1_info[23].fileTimeHi = 0U;
  c1_info[23].mFileTimeLo = 0U;
  c1_info[23].mFileTimeHi = 0U;
  c1_info[24].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m";
  c1_info[24].name = "eml_realmin";
  c1_info[24].dominantType = "char";
  c1_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m";
  c1_info[24].fileTimeLo = 1307651244U;
  c1_info[24].fileTimeHi = 0U;
  c1_info[24].mFileTimeLo = 0U;
  c1_info[24].mFileTimeHi = 0U;
  c1_info[25].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m";
  c1_info[25].name = "eml_float_model";
  c1_info[25].dominantType = "char";
  c1_info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m";
  c1_info[25].fileTimeLo = 1326727996U;
  c1_info[25].fileTimeHi = 0U;
  c1_info[25].mFileTimeLo = 0U;
  c1_info[25].mFileTimeHi = 0U;
  c1_info[26].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[26].name = "eps";
  c1_info[26].dominantType = "char";
  c1_info[26].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c1_info[26].fileTimeLo = 1326727996U;
  c1_info[26].fileTimeHi = 0U;
  c1_info[26].mFileTimeLo = 0U;
  c1_info[26].mFileTimeHi = 0U;
  c1_info[27].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c1_info[27].name = "eml_is_float_class";
  c1_info[27].dominantType = "char";
  c1_info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  c1_info[27].fileTimeLo = 1286818782U;
  c1_info[27].fileTimeHi = 0U;
  c1_info[27].mFileTimeLo = 0U;
  c1_info[27].mFileTimeHi = 0U;
  c1_info[28].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c1_info[28].name = "eml_eps";
  c1_info[28].dominantType = "char";
  c1_info[28].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m";
  c1_info[28].fileTimeLo = 1326727996U;
  c1_info[28].fileTimeHi = 0U;
  c1_info[28].mFileTimeLo = 0U;
  c1_info[28].mFileTimeHi = 0U;
  c1_info[29].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m";
  c1_info[29].name = "eml_float_model";
  c1_info[29].dominantType = "char";
  c1_info[29].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m";
  c1_info[29].fileTimeLo = 1326727996U;
  c1_info[29].fileTimeHi = 0U;
  c1_info[29].mFileTimeLo = 0U;
  c1_info[29].mFileTimeHi = 0U;
  c1_info[30].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[30].name = "min";
  c1_info[30].dominantType = "coder.internal.indexInt";
  c1_info[30].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c1_info[30].fileTimeLo = 1311255318U;
  c1_info[30].fileTimeHi = 0U;
  c1_info[30].mFileTimeLo = 0U;
  c1_info[30].mFileTimeHi = 0U;
  c1_info[31].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c1_info[31].name = "eml_min_or_max";
  c1_info[31].dominantType = "char";
  c1_info[31].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  c1_info[31].fileTimeLo = 1334071490U;
  c1_info[31].fileTimeHi = 0U;
  c1_info[31].mFileTimeLo = 0U;
  c1_info[31].mFileTimeHi = 0U;
  c1_info[32].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c1_info[32].name = "eml_scalar_eg";
  c1_info[32].dominantType = "coder.internal.indexInt";
  c1_info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[32].fileTimeLo = 1286818796U;
  c1_info[32].fileTimeHi = 0U;
  c1_info[32].mFileTimeLo = 0U;
  c1_info[32].mFileTimeHi = 0U;
  c1_info[33].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c1_info[33].name = "eml_scalexp_alloc";
  c1_info[33].dominantType = "coder.internal.indexInt";
  c1_info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c1_info[33].fileTimeLo = 1352424860U;
  c1_info[33].fileTimeHi = 0U;
  c1_info[33].mFileTimeLo = 0U;
  c1_info[33].mFileTimeHi = 0U;
  c1_info[34].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c1_info[34].name = "eml_index_class";
  c1_info[34].dominantType = "";
  c1_info[34].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[34].fileTimeLo = 1323170578U;
  c1_info[34].fileTimeHi = 0U;
  c1_info[34].mFileTimeLo = 0U;
  c1_info[34].mFileTimeHi = 0U;
  c1_info[35].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum";
  c1_info[35].name = "eml_scalar_eg";
  c1_info[35].dominantType = "coder.internal.indexInt";
  c1_info[35].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[35].fileTimeLo = 1286818796U;
  c1_info[35].fileTimeHi = 0U;
  c1_info[35].mFileTimeLo = 0U;
  c1_info[35].mFileTimeHi = 0U;
  c1_info[36].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[36].name = "colon";
  c1_info[36].dominantType = "double";
  c1_info[36].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c1_info[36].fileTimeLo = 1348191928U;
  c1_info[36].fileTimeHi = 0U;
  c1_info[36].mFileTimeLo = 0U;
  c1_info[36].mFileTimeHi = 0U;
  c1_info[37].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c1_info[37].name = "colon";
  c1_info[37].dominantType = "double";
  c1_info[37].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c1_info[37].fileTimeLo = 1348191928U;
  c1_info[37].fileTimeHi = 0U;
  c1_info[37].mFileTimeLo = 0U;
  c1_info[37].mFileTimeHi = 0U;
  c1_info[38].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c1_info[38].name = "floor";
  c1_info[38].dominantType = "double";
  c1_info[38].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c1_info[38].fileTimeLo = 1343830380U;
  c1_info[38].fileTimeHi = 0U;
  c1_info[38].mFileTimeLo = 0U;
  c1_info[38].mFileTimeHi = 0U;
  c1_info[39].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c1_info[39].name = "eml_scalar_floor";
  c1_info[39].dominantType = "double";
  c1_info[39].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  c1_info[39].fileTimeLo = 1286818726U;
  c1_info[39].fileTimeHi = 0U;
  c1_info[39].mFileTimeLo = 0U;
  c1_info[39].mFileTimeHi = 0U;
  c1_info[40].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  c1_info[40].name = "intmin";
  c1_info[40].dominantType = "char";
  c1_info[40].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c1_info[40].fileTimeLo = 1311255318U;
  c1_info[40].fileTimeHi = 0U;
  c1_info[40].mFileTimeLo = 0U;
  c1_info[40].mFileTimeHi = 0U;
  c1_info[41].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  c1_info[41].name = "intmax";
  c1_info[41].dominantType = "char";
  c1_info[41].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c1_info[41].fileTimeLo = 1311255316U;
  c1_info[41].fileTimeHi = 0U;
  c1_info[41].mFileTimeLo = 0U;
  c1_info[41].mFileTimeHi = 0U;
  c1_info[42].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  c1_info[42].name = "intmin";
  c1_info[42].dominantType = "char";
  c1_info[42].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c1_info[42].fileTimeLo = 1311255318U;
  c1_info[42].fileTimeHi = 0U;
  c1_info[42].mFileTimeLo = 0U;
  c1_info[42].mFileTimeHi = 0U;
  c1_info[43].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  c1_info[43].name = "intmax";
  c1_info[43].dominantType = "char";
  c1_info[43].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c1_info[43].fileTimeLo = 1311255316U;
  c1_info[43].fileTimeHi = 0U;
  c1_info[43].mFileTimeLo = 0U;
  c1_info[43].mFileTimeHi = 0U;
  c1_info[44].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  c1_info[44].name = "eml_isa_uint";
  c1_info[44].dominantType = "coder.internal.indexInt";
  c1_info[44].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m";
  c1_info[44].fileTimeLo = 1286818784U;
  c1_info[44].fileTimeHi = 0U;
  c1_info[44].mFileTimeLo = 0U;
  c1_info[44].mFileTimeHi = 0U;
  c1_info[45].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c1_info[45].name = "eml_unsigned_class";
  c1_info[45].dominantType = "char";
  c1_info[45].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  c1_info[45].fileTimeLo = 1323170580U;
  c1_info[45].fileTimeHi = 0U;
  c1_info[45].mFileTimeLo = 0U;
  c1_info[45].mFileTimeHi = 0U;
  c1_info[46].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  c1_info[46].name = "eml_index_class";
  c1_info[46].dominantType = "";
  c1_info[46].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[46].fileTimeLo = 1323170578U;
  c1_info[46].fileTimeHi = 0U;
  c1_info[46].mFileTimeLo = 0U;
  c1_info[46].mFileTimeHi = 0U;
  c1_info[47].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c1_info[47].name = "eml_index_class";
  c1_info[47].dominantType = "";
  c1_info[47].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[47].fileTimeLo = 1323170578U;
  c1_info[47].fileTimeHi = 0U;
  c1_info[47].mFileTimeLo = 0U;
  c1_info[47].mFileTimeHi = 0U;
  c1_info[48].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c1_info[48].name = "intmax";
  c1_info[48].dominantType = "char";
  c1_info[48].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c1_info[48].fileTimeLo = 1311255316U;
  c1_info[48].fileTimeHi = 0U;
  c1_info[48].mFileTimeLo = 0U;
  c1_info[48].mFileTimeHi = 0U;
  c1_info[49].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c1_info[49].name = "eml_isa_uint";
  c1_info[49].dominantType = "coder.internal.indexInt";
  c1_info[49].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m";
  c1_info[49].fileTimeLo = 1286818784U;
  c1_info[49].fileTimeHi = 0U;
  c1_info[49].mFileTimeLo = 0U;
  c1_info[49].mFileTimeHi = 0U;
  c1_info[50].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c1_info[50].name = "eml_index_plus";
  c1_info[50].dominantType = "double";
  c1_info[50].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[50].fileTimeLo = 1286818778U;
  c1_info[50].fileTimeHi = 0U;
  c1_info[50].mFileTimeLo = 0U;
  c1_info[50].mFileTimeHi = 0U;
  c1_info[51].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[51].name = "eml_index_class";
  c1_info[51].dominantType = "";
  c1_info[51].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[51].fileTimeLo = 1323170578U;
  c1_info[51].fileTimeHi = 0U;
  c1_info[51].mFileTimeLo = 0U;
  c1_info[51].mFileTimeHi = 0U;
  c1_info[52].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon";
  c1_info[52].name = "eml_int_forloop_overflow_check";
  c1_info[52].dominantType = "";
  c1_info[52].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c1_info[52].fileTimeLo = 1346510340U;
  c1_info[52].fileTimeHi = 0U;
  c1_info[52].mFileTimeLo = 0U;
  c1_info[52].mFileTimeHi = 0U;
  c1_info[53].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  c1_info[53].name = "intmax";
  c1_info[53].dominantType = "char";
  c1_info[53].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c1_info[53].fileTimeLo = 1311255316U;
  c1_info[53].fileTimeHi = 0U;
  c1_info[53].mFileTimeLo = 0U;
  c1_info[53].mFileTimeHi = 0U;
  c1_info[54].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[54].name = "eml_index_class";
  c1_info[54].dominantType = "";
  c1_info[54].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[54].fileTimeLo = 1323170578U;
  c1_info[54].fileTimeHi = 0U;
  c1_info[54].mFileTimeLo = 0U;
  c1_info[54].mFileTimeHi = 0U;
  c1_info[55].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[55].name = "eml_index_plus";
  c1_info[55].dominantType = "double";
  c1_info[55].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[55].fileTimeLo = 1286818778U;
  c1_info[55].fileTimeHi = 0U;
  c1_info[55].mFileTimeLo = 0U;
  c1_info[55].mFileTimeHi = 0U;
  c1_info[56].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[56].name = "eml_int_forloop_overflow_check";
  c1_info[56].dominantType = "";
  c1_info[56].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c1_info[56].fileTimeLo = 1346510340U;
  c1_info[56].fileTimeHi = 0U;
  c1_info[56].mFileTimeLo = 0U;
  c1_info[56].mFileTimeHi = 0U;
  c1_info[57].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[57].name = "eml_index_minus";
  c1_info[57].dominantType = "double";
  c1_info[57].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c1_info[57].fileTimeLo = 1286818778U;
  c1_info[57].fileTimeHi = 0U;
  c1_info[57].mFileTimeLo = 0U;
  c1_info[57].mFileTimeHi = 0U;
  c1_info[58].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c1_info[58].name = "eml_index_class";
  c1_info[58].dominantType = "";
  c1_info[58].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[58].fileTimeLo = 1323170578U;
  c1_info[58].fileTimeHi = 0U;
  c1_info[58].mFileTimeLo = 0U;
  c1_info[58].mFileTimeHi = 0U;
  c1_info[59].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[59].name = "eml_index_minus";
  c1_info[59].dominantType = "coder.internal.indexInt";
  c1_info[59].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c1_info[59].fileTimeLo = 1286818778U;
  c1_info[59].fileTimeHi = 0U;
  c1_info[59].mFileTimeLo = 0U;
  c1_info[59].mFileTimeHi = 0U;
  c1_info[60].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[60].name = "eml_index_times";
  c1_info[60].dominantType = "coder.internal.indexInt";
  c1_info[60].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c1_info[60].fileTimeLo = 1286818780U;
  c1_info[60].fileTimeHi = 0U;
  c1_info[60].mFileTimeLo = 0U;
  c1_info[60].mFileTimeHi = 0U;
  c1_info[61].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c1_info[61].name = "eml_index_class";
  c1_info[61].dominantType = "";
  c1_info[61].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[61].fileTimeLo = 1323170578U;
  c1_info[61].fileTimeHi = 0U;
  c1_info[61].mFileTimeLo = 0U;
  c1_info[61].mFileTimeHi = 0U;
  c1_info[62].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[62].name = "eml_index_plus";
  c1_info[62].dominantType = "coder.internal.indexInt";
  c1_info[62].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[62].fileTimeLo = 1286818778U;
  c1_info[62].fileTimeHi = 0U;
  c1_info[62].mFileTimeLo = 0U;
  c1_info[62].mFileTimeHi = 0U;
  c1_info[63].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[63].name = "eml_ixamax";
  c1_info[63].dominantType = "double";
  c1_info[63].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m";
  c1_info[63].fileTimeLo = 1299076770U;
  c1_info[63].fileTimeHi = 0U;
  c1_info[63].mFileTimeLo = 0U;
  c1_info[63].mFileTimeHi = 0U;
}

static void c1_b_info_helper(c1_ResolvedFunctionInfo c1_info[122])
{
  c1_info[64].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m";
  c1_info[64].name = "eml_blas_inline";
  c1_info[64].dominantType = "";
  c1_info[64].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c1_info[64].fileTimeLo = 1299076768U;
  c1_info[64].fileTimeHi = 0U;
  c1_info[64].mFileTimeLo = 0U;
  c1_info[64].mFileTimeHi = 0U;
  c1_info[65].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m!below_threshold";
  c1_info[65].name = "length";
  c1_info[65].dominantType = "double";
  c1_info[65].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  c1_info[65].fileTimeLo = 1303146206U;
  c1_info[65].fileTimeHi = 0U;
  c1_info[65].mFileTimeLo = 0U;
  c1_info[65].mFileTimeHi = 0U;
  c1_info[66].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength";
  c1_info[66].name = "eml_index_class";
  c1_info[66].dominantType = "";
  c1_info[66].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[66].fileTimeLo = 1323170578U;
  c1_info[66].fileTimeHi = 0U;
  c1_info[66].mFileTimeLo = 0U;
  c1_info[66].mFileTimeHi = 0U;
  c1_info[67].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m";
  c1_info[67].name = "eml_index_class";
  c1_info[67].dominantType = "";
  c1_info[67].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[67].fileTimeLo = 1323170578U;
  c1_info[67].fileTimeHi = 0U;
  c1_info[67].mFileTimeLo = 0U;
  c1_info[67].mFileTimeHi = 0U;
  c1_info[68].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m";
  c1_info[68].name = "eml_refblas_ixamax";
  c1_info[68].dominantType = "double";
  c1_info[68].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c1_info[68].fileTimeLo = 1299076770U;
  c1_info[68].fileTimeHi = 0U;
  c1_info[68].mFileTimeLo = 0U;
  c1_info[68].mFileTimeHi = 0U;
  c1_info[69].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c1_info[69].name = "eml_index_class";
  c1_info[69].dominantType = "";
  c1_info[69].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[69].fileTimeLo = 1323170578U;
  c1_info[69].fileTimeHi = 0U;
  c1_info[69].mFileTimeLo = 0U;
  c1_info[69].mFileTimeHi = 0U;
  c1_info[70].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c1_info[70].name = "eml_xcabs1";
  c1_info[70].dominantType = "double";
  c1_info[70].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c1_info[70].fileTimeLo = 1286818706U;
  c1_info[70].fileTimeHi = 0U;
  c1_info[70].mFileTimeLo = 0U;
  c1_info[70].mFileTimeHi = 0U;
  c1_info[71].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c1_info[71].name = "abs";
  c1_info[71].dominantType = "double";
  c1_info[71].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c1_info[71].fileTimeLo = 1343830366U;
  c1_info[71].fileTimeHi = 0U;
  c1_info[71].mFileTimeLo = 0U;
  c1_info[71].mFileTimeHi = 0U;
  c1_info[72].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c1_info[72].name = "eml_scalar_abs";
  c1_info[72].dominantType = "double";
  c1_info[72].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c1_info[72].fileTimeLo = 1286818712U;
  c1_info[72].fileTimeHi = 0U;
  c1_info[72].mFileTimeLo = 0U;
  c1_info[72].mFileTimeHi = 0U;
  c1_info[73].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c1_info[73].name = "eml_int_forloop_overflow_check";
  c1_info[73].dominantType = "";
  c1_info[73].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c1_info[73].fileTimeLo = 1346510340U;
  c1_info[73].fileTimeHi = 0U;
  c1_info[73].mFileTimeLo = 0U;
  c1_info[73].mFileTimeHi = 0U;
  c1_info[74].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c1_info[74].name = "eml_index_plus";
  c1_info[74].dominantType = "coder.internal.indexInt";
  c1_info[74].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[74].fileTimeLo = 1286818778U;
  c1_info[74].fileTimeHi = 0U;
  c1_info[74].mFileTimeLo = 0U;
  c1_info[74].mFileTimeHi = 0U;
  c1_info[75].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[75].name = "eml_xswap";
  c1_info[75].dominantType = "double";
  c1_info[75].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m";
  c1_info[75].fileTimeLo = 1299076778U;
  c1_info[75].fileTimeHi = 0U;
  c1_info[75].mFileTimeLo = 0U;
  c1_info[75].mFileTimeHi = 0U;
  c1_info[76].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m";
  c1_info[76].name = "eml_blas_inline";
  c1_info[76].dominantType = "";
  c1_info[76].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c1_info[76].fileTimeLo = 1299076768U;
  c1_info[76].fileTimeHi = 0U;
  c1_info[76].mFileTimeLo = 0U;
  c1_info[76].mFileTimeHi = 0U;
  c1_info[77].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m";
  c1_info[77].name = "eml_index_class";
  c1_info[77].dominantType = "";
  c1_info[77].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[77].fileTimeLo = 1323170578U;
  c1_info[77].fileTimeHi = 0U;
  c1_info[77].mFileTimeLo = 0U;
  c1_info[77].mFileTimeHi = 0U;
  c1_info[78].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m";
  c1_info[78].name = "eml_refblas_xswap";
  c1_info[78].dominantType = "double";
  c1_info[78].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c1_info[78].fileTimeLo = 1299076786U;
  c1_info[78].fileTimeHi = 0U;
  c1_info[78].mFileTimeLo = 0U;
  c1_info[78].mFileTimeHi = 0U;
  c1_info[79].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c1_info[79].name = "eml_index_class";
  c1_info[79].dominantType = "";
  c1_info[79].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[79].fileTimeLo = 1323170578U;
  c1_info[79].fileTimeHi = 0U;
  c1_info[79].mFileTimeLo = 0U;
  c1_info[79].mFileTimeHi = 0U;
  c1_info[80].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c1_info[80].name = "abs";
  c1_info[80].dominantType = "coder.internal.indexInt";
  c1_info[80].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c1_info[80].fileTimeLo = 1343830366U;
  c1_info[80].fileTimeHi = 0U;
  c1_info[80].mFileTimeLo = 0U;
  c1_info[80].mFileTimeHi = 0U;
  c1_info[81].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c1_info[81].name = "eml_scalar_abs";
  c1_info[81].dominantType = "coder.internal.indexInt";
  c1_info[81].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c1_info[81].fileTimeLo = 1286818712U;
  c1_info[81].fileTimeHi = 0U;
  c1_info[81].mFileTimeLo = 0U;
  c1_info[81].mFileTimeHi = 0U;
  c1_info[82].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c1_info[82].name = "eml_int_forloop_overflow_check";
  c1_info[82].dominantType = "";
  c1_info[82].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c1_info[82].fileTimeLo = 1346510340U;
  c1_info[82].fileTimeHi = 0U;
  c1_info[82].mFileTimeLo = 0U;
  c1_info[82].mFileTimeHi = 0U;
  c1_info[83].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c1_info[83].name = "eml_index_plus";
  c1_info[83].dominantType = "coder.internal.indexInt";
  c1_info[83].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[83].fileTimeLo = 1286818778U;
  c1_info[83].fileTimeHi = 0U;
  c1_info[83].mFileTimeLo = 0U;
  c1_info[83].mFileTimeHi = 0U;
  c1_info[84].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[84].name = "eml_div";
  c1_info[84].dominantType = "double";
  c1_info[84].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c1_info[84].fileTimeLo = 1313347810U;
  c1_info[84].fileTimeHi = 0U;
  c1_info[84].mFileTimeLo = 0U;
  c1_info[84].mFileTimeHi = 0U;
  c1_info[85].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c1_info[85].name = "eml_xgeru";
  c1_info[85].dominantType = "double";
  c1_info[85].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c1_info[85].fileTimeLo = 1299076774U;
  c1_info[85].fileTimeHi = 0U;
  c1_info[85].mFileTimeLo = 0U;
  c1_info[85].mFileTimeHi = 0U;
  c1_info[86].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c1_info[86].name = "eml_blas_inline";
  c1_info[86].dominantType = "";
  c1_info[86].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c1_info[86].fileTimeLo = 1299076768U;
  c1_info[86].fileTimeHi = 0U;
  c1_info[86].mFileTimeLo = 0U;
  c1_info[86].mFileTimeHi = 0U;
  c1_info[87].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c1_info[87].name = "eml_xger";
  c1_info[87].dominantType = "double";
  c1_info[87].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m";
  c1_info[87].fileTimeLo = 1299076774U;
  c1_info[87].fileTimeHi = 0U;
  c1_info[87].mFileTimeLo = 0U;
  c1_info[87].mFileTimeHi = 0U;
  c1_info[88].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m";
  c1_info[88].name = "eml_blas_inline";
  c1_info[88].dominantType = "";
  c1_info[88].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c1_info[88].fileTimeLo = 1299076768U;
  c1_info[88].fileTimeHi = 0U;
  c1_info[88].mFileTimeLo = 0U;
  c1_info[88].mFileTimeHi = 0U;
  c1_info[89].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold";
  c1_info[89].name = "intmax";
  c1_info[89].dominantType = "char";
  c1_info[89].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c1_info[89].fileTimeLo = 1311255316U;
  c1_info[89].fileTimeHi = 0U;
  c1_info[89].mFileTimeLo = 0U;
  c1_info[89].mFileTimeHi = 0U;
  c1_info[90].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold";
  c1_info[90].name = "min";
  c1_info[90].dominantType = "double";
  c1_info[90].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c1_info[90].fileTimeLo = 1311255318U;
  c1_info[90].fileTimeHi = 0U;
  c1_info[90].mFileTimeLo = 0U;
  c1_info[90].mFileTimeHi = 0U;
  c1_info[91].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c1_info[91].name = "eml_scalar_eg";
  c1_info[91].dominantType = "double";
  c1_info[91].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[91].fileTimeLo = 1286818796U;
  c1_info[91].fileTimeHi = 0U;
  c1_info[91].mFileTimeLo = 0U;
  c1_info[91].mFileTimeHi = 0U;
  c1_info[92].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c1_info[92].name = "eml_scalexp_alloc";
  c1_info[92].dominantType = "double";
  c1_info[92].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c1_info[92].fileTimeLo = 1352424860U;
  c1_info[92].fileTimeHi = 0U;
  c1_info[92].mFileTimeLo = 0U;
  c1_info[92].mFileTimeHi = 0U;
  c1_info[93].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum";
  c1_info[93].name = "eml_scalar_eg";
  c1_info[93].dominantType = "double";
  c1_info[93].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[93].fileTimeLo = 1286818796U;
  c1_info[93].fileTimeHi = 0U;
  c1_info[93].mFileTimeLo = 0U;
  c1_info[93].mFileTimeHi = 0U;
  c1_info[94].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold";
  c1_info[94].name = "mtimes";
  c1_info[94].dominantType = "double";
  c1_info[94].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[94].fileTimeLo = 1289519692U;
  c1_info[94].fileTimeHi = 0U;
  c1_info[94].mFileTimeLo = 0U;
  c1_info[94].mFileTimeHi = 0U;
  c1_info[95].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m";
  c1_info[95].name = "eml_index_class";
  c1_info[95].dominantType = "";
  c1_info[95].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[95].fileTimeLo = 1323170578U;
  c1_info[95].fileTimeHi = 0U;
  c1_info[95].mFileTimeLo = 0U;
  c1_info[95].mFileTimeHi = 0U;
  c1_info[96].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m";
  c1_info[96].name = "eml_refblas_xger";
  c1_info[96].dominantType = "double";
  c1_info[96].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m";
  c1_info[96].fileTimeLo = 1299076776U;
  c1_info[96].fileTimeHi = 0U;
  c1_info[96].mFileTimeLo = 0U;
  c1_info[96].mFileTimeHi = 0U;
  c1_info[97].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m";
  c1_info[97].name = "eml_refblas_xgerx";
  c1_info[97].dominantType = "char";
  c1_info[97].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[97].fileTimeLo = 1299076778U;
  c1_info[97].fileTimeHi = 0U;
  c1_info[97].mFileTimeLo = 0U;
  c1_info[97].mFileTimeHi = 0U;
  c1_info[98].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[98].name = "eml_index_class";
  c1_info[98].dominantType = "";
  c1_info[98].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[98].fileTimeLo = 1323170578U;
  c1_info[98].fileTimeHi = 0U;
  c1_info[98].mFileTimeLo = 0U;
  c1_info[98].mFileTimeHi = 0U;
  c1_info[99].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[99].name = "abs";
  c1_info[99].dominantType = "coder.internal.indexInt";
  c1_info[99].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c1_info[99].fileTimeLo = 1343830366U;
  c1_info[99].fileTimeHi = 0U;
  c1_info[99].mFileTimeLo = 0U;
  c1_info[99].mFileTimeHi = 0U;
  c1_info[100].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[100].name = "eml_index_minus";
  c1_info[100].dominantType = "double";
  c1_info[100].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c1_info[100].fileTimeLo = 1286818778U;
  c1_info[100].fileTimeHi = 0U;
  c1_info[100].mFileTimeLo = 0U;
  c1_info[100].mFileTimeHi = 0U;
  c1_info[101].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[101].name = "eml_int_forloop_overflow_check";
  c1_info[101].dominantType = "";
  c1_info[101].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c1_info[101].fileTimeLo = 1346510340U;
  c1_info[101].fileTimeHi = 0U;
  c1_info[101].mFileTimeLo = 0U;
  c1_info[101].mFileTimeHi = 0U;
  c1_info[102].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[102].name = "eml_index_plus";
  c1_info[102].dominantType = "double";
  c1_info[102].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[102].fileTimeLo = 1286818778U;
  c1_info[102].fileTimeHi = 0U;
  c1_info[102].mFileTimeLo = 0U;
  c1_info[102].mFileTimeHi = 0U;
  c1_info[103].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c1_info[103].name = "eml_index_plus";
  c1_info[103].dominantType = "coder.internal.indexInt";
  c1_info[103].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[103].fileTimeLo = 1286818778U;
  c1_info[103].fileTimeHi = 0U;
  c1_info[103].mFileTimeLo = 0U;
  c1_info[103].mFileTimeHi = 0U;
  c1_info[104].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular";
  c1_info[104].name = "eml_warning";
  c1_info[104].dominantType = "char";
  c1_info[104].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m";
  c1_info[104].fileTimeLo = 1286818802U;
  c1_info[104].fileTimeHi = 0U;
  c1_info[104].mFileTimeLo = 0U;
  c1_info[104].mFileTimeHi = 0U;
  c1_info[105].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c1_info[105].name = "eml_scalar_eg";
  c1_info[105].dominantType = "double";
  c1_info[105].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[105].fileTimeLo = 1286818796U;
  c1_info[105].fileTimeHi = 0U;
  c1_info[105].mFileTimeLo = 0U;
  c1_info[105].mFileTimeHi = 0U;
  c1_info[106].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c1_info[106].name = "eml_int_forloop_overflow_check";
  c1_info[106].dominantType = "";
  c1_info[106].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c1_info[106].fileTimeLo = 1346510340U;
  c1_info[106].fileTimeHi = 0U;
  c1_info[106].mFileTimeLo = 0U;
  c1_info[106].mFileTimeHi = 0U;
  c1_info[107].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c1_info[107].name = "eml_xtrsm";
  c1_info[107].dominantType = "char";
  c1_info[107].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m";
  c1_info[107].fileTimeLo = 1299076778U;
  c1_info[107].fileTimeHi = 0U;
  c1_info[107].mFileTimeLo = 0U;
  c1_info[107].mFileTimeHi = 0U;
  c1_info[108].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m";
  c1_info[108].name = "eml_blas_inline";
  c1_info[108].dominantType = "";
  c1_info[108].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c1_info[108].fileTimeLo = 1299076768U;
  c1_info[108].fileTimeHi = 0U;
  c1_info[108].mFileTimeLo = 0U;
  c1_info[108].mFileTimeHi = 0U;
  c1_info[109].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m!below_threshold";
  c1_info[109].name = "mtimes";
  c1_info[109].dominantType = "double";
  c1_info[109].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c1_info[109].fileTimeLo = 1289519692U;
  c1_info[109].fileTimeHi = 0U;
  c1_info[109].mFileTimeLo = 0U;
  c1_info[109].mFileTimeHi = 0U;
  c1_info[110].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c1_info[110].name = "eml_index_class";
  c1_info[110].dominantType = "";
  c1_info[110].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[110].fileTimeLo = 1323170578U;
  c1_info[110].fileTimeHi = 0U;
  c1_info[110].mFileTimeLo = 0U;
  c1_info[110].mFileTimeHi = 0U;
  c1_info[111].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c1_info[111].name = "eml_scalar_eg";
  c1_info[111].dominantType = "double";
  c1_info[111].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[111].fileTimeLo = 1286818796U;
  c1_info[111].fileTimeHi = 0U;
  c1_info[111].mFileTimeLo = 0U;
  c1_info[111].mFileTimeHi = 0U;
  c1_info[112].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c1_info[112].name = "eml_refblas_xtrsm";
  c1_info[112].dominantType = "char";
  c1_info[112].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[112].fileTimeLo = 1299076786U;
  c1_info[112].fileTimeHi = 0U;
  c1_info[112].mFileTimeLo = 0U;
  c1_info[112].mFileTimeHi = 0U;
  c1_info[113].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[113].name = "eml_scalar_eg";
  c1_info[113].dominantType = "double";
  c1_info[113].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c1_info[113].fileTimeLo = 1286818796U;
  c1_info[113].fileTimeHi = 0U;
  c1_info[113].mFileTimeLo = 0U;
  c1_info[113].mFileTimeHi = 0U;
  c1_info[114].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[114].name = "eml_index_minus";
  c1_info[114].dominantType = "double";
  c1_info[114].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c1_info[114].fileTimeLo = 1286818778U;
  c1_info[114].fileTimeHi = 0U;
  c1_info[114].mFileTimeLo = 0U;
  c1_info[114].mFileTimeHi = 0U;
  c1_info[115].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[115].name = "eml_index_class";
  c1_info[115].dominantType = "";
  c1_info[115].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c1_info[115].fileTimeLo = 1323170578U;
  c1_info[115].fileTimeHi = 0U;
  c1_info[115].mFileTimeLo = 0U;
  c1_info[115].mFileTimeHi = 0U;
  c1_info[116].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[116].name = "eml_int_forloop_overflow_check";
  c1_info[116].dominantType = "";
  c1_info[116].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c1_info[116].fileTimeLo = 1346510340U;
  c1_info[116].fileTimeHi = 0U;
  c1_info[116].mFileTimeLo = 0U;
  c1_info[116].mFileTimeHi = 0U;
  c1_info[117].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[117].name = "eml_index_times";
  c1_info[117].dominantType = "coder.internal.indexInt";
  c1_info[117].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c1_info[117].fileTimeLo = 1286818780U;
  c1_info[117].fileTimeHi = 0U;
  c1_info[117].mFileTimeLo = 0U;
  c1_info[117].mFileTimeHi = 0U;
  c1_info[118].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[118].name = "eml_index_plus";
  c1_info[118].dominantType = "coder.internal.indexInt";
  c1_info[118].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[118].fileTimeLo = 1286818778U;
  c1_info[118].fileTimeHi = 0U;
  c1_info[118].mFileTimeLo = 0U;
  c1_info[118].mFileTimeHi = 0U;
  c1_info[119].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[119].name = "eml_index_plus";
  c1_info[119].dominantType = "double";
  c1_info[119].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c1_info[119].fileTimeLo = 1286818778U;
  c1_info[119].fileTimeHi = 0U;
  c1_info[119].mFileTimeLo = 0U;
  c1_info[119].mFileTimeHi = 0U;
  c1_info[120].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  c1_info[120].name = "intmin";
  c1_info[120].dominantType = "char";
  c1_info[120].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c1_info[120].fileTimeLo = 1311255318U;
  c1_info[120].fileTimeHi = 0U;
  c1_info[120].mFileTimeLo = 0U;
  c1_info[120].mFileTimeHi = 0U;
  c1_info[121].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c1_info[121].name = "eml_div";
  c1_info[121].dominantType = "double";
  c1_info[121].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c1_info[121].fileTimeLo = 1313347810U;
  c1_info[121].fileTimeHi = 0U;
  c1_info[121].mFileTimeLo = 0U;
  c1_info[121].mFileTimeHi = 0U;
}

static void c1_eml_scalar_eg(SFc1_testInstanceStruct *chartInstance)
{
}

static void c1_realmin(SFc1_testInstanceStruct *chartInstance)
{
}

static void c1_eps(SFc1_testInstanceStruct *chartInstance)
{
}

static void c1_eml_matlab_zgetrf(SFc1_testInstanceStruct *chartInstance, real_T
  c1_A[16], real_T c1_b_A[16], int32_T c1_ipiv[4], int32_T *c1_info)
{
  int32_T c1_i73;
  for (c1_i73 = 0; c1_i73 < 16; c1_i73++) {
    c1_b_A[c1_i73] = c1_A[c1_i73];
  }

  c1_b_eml_matlab_zgetrf(chartInstance, c1_b_A, c1_ipiv, c1_info);
}

static void c1_check_forloop_overflow_error(SFc1_testInstanceStruct
  *chartInstance, boolean_T c1_overflow)
{
  int32_T c1_i74;
  static char_T c1_cv0[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c1_u[34];
  const mxArray *c1_y = NULL;
  int32_T c1_i75;
  static char_T c1_cv1[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c1_b_u[23];
  const mxArray *c1_b_y = NULL;
  if (!c1_overflow) {
  } else {
    for (c1_i74 = 0; c1_i74 < 34; c1_i74++) {
      c1_u[c1_i74] = c1_cv0[c1_i74];
    }

    c1_y = NULL;
    sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  FALSE);
    for (c1_i75 = 0; c1_i75 < 23; c1_i75++) {
      c1_b_u[c1_i75] = c1_cv1[c1_i75];
    }

    c1_b_y = NULL;
    sf_mex_assign(&c1_b_y, sf_mex_create("y", c1_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U,
      14, c1_y, 14, c1_b_y));
  }
}

static void c1_eml_xger(SFc1_testInstanceStruct *chartInstance, int32_T c1_m,
  int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0, int32_T c1_iy0, real_T c1_A[16],
  int32_T c1_ia0, real_T c1_b_A[16])
{
  int32_T c1_i76;
  for (c1_i76 = 0; c1_i76 < 16; c1_i76++) {
    c1_b_A[c1_i76] = c1_A[c1_i76];
  }

  c1_b_eml_xger(chartInstance, c1_m, c1_n, c1_alpha1, c1_ix0, c1_iy0, c1_b_A,
                c1_ia0);
}

static void c1_eml_warning(SFc1_testInstanceStruct *chartInstance)
{
  int32_T c1_i77;
  static char_T c1_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c1_u[27];
  const mxArray *c1_y = NULL;
  for (c1_i77 = 0; c1_i77 < 27; c1_i77++) {
    c1_u[c1_i77] = c1_varargin_1[c1_i77];
  }

  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", c1_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
    14, c1_y));
}

static void c1_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A[16],
  real_T c1_B[24], real_T c1_b_B[24])
{
  int32_T c1_i78;
  int32_T c1_i79;
  real_T c1_b_A[16];
  for (c1_i78 = 0; c1_i78 < 24; c1_i78++) {
    c1_b_B[c1_i78] = c1_B[c1_i78];
  }

  for (c1_i79 = 0; c1_i79 < 16; c1_i79++) {
    c1_b_A[c1_i79] = c1_A[c1_i79];
  }

  c1_c_eml_xtrsm(chartInstance, c1_b_A, c1_b_B);
}

static void c1_below_threshold(SFc1_testInstanceStruct *chartInstance)
{
}

static void c1_b_eml_scalar_eg(SFc1_testInstanceStruct *chartInstance)
{
}

static void c1_b_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A
  [16], real_T c1_B[24], real_T c1_b_B[24])
{
  int32_T c1_i80;
  int32_T c1_i81;
  real_T c1_b_A[16];
  for (c1_i80 = 0; c1_i80 < 24; c1_i80++) {
    c1_b_B[c1_i80] = c1_B[c1_i80];
  }

  for (c1_i81 = 0; c1_i81 < 16; c1_i81++) {
    c1_b_A[c1_i81] = c1_A[c1_i81];
  }

  c1_d_eml_xtrsm(chartInstance, c1_b_A, c1_b_B);
}

static void c1_c_eml_scalar_eg(SFc1_testInstanceStruct *chartInstance)
{
}

static const mxArray *c1_f_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static int32_T c1_f_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i82;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i82, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i82;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_b_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)chartInstanceVoid;
  c1_b_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_sfEvent),
    &c1_thisId);
  sf_mex_destroy(&c1_b_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_g_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_b_is_active_c1_test, const char_T *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_b_is_active_c1_test),
    &c1_thisId);
  sf_mex_destroy(&c1_b_is_active_c1_test);
  return c1_y;
}

static uint8_T c1_h_emlrt_marshallIn(SFc1_testInstanceStruct *chartInstance,
  const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_eml_matlab_zgetrf(SFc1_testInstanceStruct *chartInstance,
  real_T c1_A[16], int32_T c1_ipiv[4], int32_T *c1_info)
{
  int32_T c1_i83;
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_a;
  int32_T c1_jm1;
  int32_T c1_b;
  int32_T c1_mmj;
  int32_T c1_b_a;
  int32_T c1_c;
  int32_T c1_b_b;
  int32_T c1_jj;
  int32_T c1_c_a;
  int32_T c1_jp1j;
  int32_T c1_d_a;
  int32_T c1_b_c;
  int32_T c1_n;
  int32_T c1_ix0;
  int32_T c1_b_n;
  int32_T c1_b_ix0;
  int32_T c1_c_n;
  int32_T c1_c_ix0;
  int32_T c1_idxmax;
  int32_T c1_ix;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_c_x;
  real_T c1_y;
  real_T c1_d_x;
  real_T c1_e_x;
  real_T c1_b_y;
  real_T c1_smax;
  int32_T c1_d_n;
  int32_T c1_c_b;
  int32_T c1_d_b;
  boolean_T c1_overflow;
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_e_a;
  real_T c1_f_x;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_c_y;
  real_T c1_i_x;
  real_T c1_j_x;
  real_T c1_d_y;
  real_T c1_s;
  int32_T c1_f_a;
  int32_T c1_jpiv_offset;
  int32_T c1_g_a;
  int32_T c1_e_b;
  int32_T c1_jpiv;
  int32_T c1_h_a;
  int32_T c1_f_b;
  int32_T c1_c_c;
  int32_T c1_g_b;
  int32_T c1_jrow;
  int32_T c1_i_a;
  int32_T c1_h_b;
  int32_T c1_jprow;
  int32_T c1_d_ix0;
  int32_T c1_iy0;
  int32_T c1_e_ix0;
  int32_T c1_b_iy0;
  int32_T c1_f_ix0;
  int32_T c1_c_iy0;
  int32_T c1_b_ix;
  int32_T c1_iy;
  int32_T c1_c_k;
  real_T c1_temp;
  int32_T c1_j_a;
  int32_T c1_k_a;
  int32_T c1_b_jp1j;
  int32_T c1_l_a;
  int32_T c1_d_c;
  int32_T c1_m_a;
  int32_T c1_i_b;
  int32_T c1_i84;
  int32_T c1_n_a;
  int32_T c1_j_b;
  int32_T c1_o_a;
  int32_T c1_k_b;
  boolean_T c1_b_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  real_T c1_k_x;
  real_T c1_e_y;
  real_T c1_z;
  int32_T c1_l_b;
  int32_T c1_e_c;
  int32_T c1_p_a;
  int32_T c1_f_c;
  int32_T c1_q_a;
  int32_T c1_g_c;
  int32_T c1_m;
  int32_T c1_e_n;
  int32_T c1_g_ix0;
  int32_T c1_d_iy0;
  int32_T c1_ia0;
  real_T c1_d1;
  c1_realmin(chartInstance);
  c1_eps(chartInstance);
  for (c1_i83 = 0; c1_i83 < 4; c1_i83++) {
    c1_ipiv[c1_i83] = 1 + c1_i83;
  }

  *c1_info = 0;
  for (c1_j = 1; c1_j < 4; c1_j++) {
    c1_b_j = c1_j;
    c1_a = c1_b_j - 1;
    c1_jm1 = c1_a;
    c1_b = c1_b_j;
    c1_mmj = 4 - c1_b;
    c1_b_a = c1_jm1;
    c1_c = c1_b_a * 5;
    c1_b_b = c1_c + 1;
    c1_jj = c1_b_b;
    c1_c_a = c1_jj + 1;
    c1_jp1j = c1_c_a;
    c1_d_a = c1_mmj;
    c1_b_c = c1_d_a;
    c1_n = c1_b_c + 1;
    c1_ix0 = c1_jj;
    c1_b_n = c1_n;
    c1_b_ix0 = c1_ix0;
    c1_c_n = c1_b_n;
    c1_c_ix0 = c1_b_ix0;
    if (c1_c_n < 1) {
      c1_idxmax = 0;
    } else {
      c1_idxmax = 1;
      if (c1_c_n > 1) {
        c1_ix = c1_c_ix0;
        c1_x = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c1_ix), 1, 16, 1, 0) - 1];
        c1_b_x = c1_x;
        c1_c_x = c1_b_x;
        c1_y = muDoubleScalarAbs(c1_c_x);
        c1_d_x = 0.0;
        c1_e_x = c1_d_x;
        c1_b_y = muDoubleScalarAbs(c1_e_x);
        c1_smax = c1_y + c1_b_y;
        c1_d_n = c1_c_n;
        c1_c_b = c1_d_n;
        c1_d_b = c1_c_b;
        if (2 > c1_d_b) {
          c1_overflow = FALSE;
        } else {
          c1_overflow = (c1_d_b > 2147483646);
        }

        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }

        for (c1_k = 2; c1_k <= c1_d_n; c1_k++) {
          c1_b_k = c1_k;
          c1_e_a = c1_ix + 1;
          c1_ix = c1_e_a;
          c1_f_x = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_ix), 1, 16, 1, 0) - 1];
          c1_g_x = c1_f_x;
          c1_h_x = c1_g_x;
          c1_c_y = muDoubleScalarAbs(c1_h_x);
          c1_i_x = 0.0;
          c1_j_x = c1_i_x;
          c1_d_y = muDoubleScalarAbs(c1_j_x);
          c1_s = c1_c_y + c1_d_y;
          if (c1_s > c1_smax) {
            c1_idxmax = c1_b_k;
            c1_smax = c1_s;
          }
        }
      }
    }

    c1_f_a = c1_idxmax - 1;
    c1_jpiv_offset = c1_f_a;
    c1_g_a = c1_jj;
    c1_e_b = c1_jpiv_offset;
    c1_jpiv = c1_g_a + c1_e_b;
    if (c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_jpiv), 1, 16, 1, 0) - 1] != 0.0) {
      if (c1_jpiv_offset != 0) {
        c1_h_a = c1_b_j;
        c1_f_b = c1_jpiv_offset;
        c1_c_c = c1_h_a + c1_f_b;
        c1_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_j), 1, 4, 1, 0) - 1] = c1_c_c;
        c1_g_b = c1_jm1 + 1;
        c1_jrow = c1_g_b;
        c1_i_a = c1_jrow;
        c1_h_b = c1_jpiv_offset;
        c1_jprow = c1_i_a + c1_h_b;
        c1_d_ix0 = c1_jrow;
        c1_iy0 = c1_jprow;
        c1_e_ix0 = c1_d_ix0;
        c1_b_iy0 = c1_iy0;
        c1_f_ix0 = c1_e_ix0;
        c1_c_iy0 = c1_b_iy0;
        c1_b_ix = c1_f_ix0;
        c1_iy = c1_c_iy0;
        for (c1_c_k = 1; c1_c_k < 5; c1_c_k++) {
          c1_temp = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c1_b_ix), 1, 16, 1, 0) - 1];
          c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ix), 1, 16, 1, 0) - 1] =
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_iy), 1, 16, 1, 0) - 1];
          c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_iy), 1, 16, 1, 0) - 1] = c1_temp;
          c1_j_a = c1_b_ix + 4;
          c1_b_ix = c1_j_a;
          c1_k_a = c1_iy + 4;
          c1_iy = c1_k_a;
        }
      }

      c1_b_jp1j = c1_jp1j;
      c1_l_a = c1_mmj;
      c1_d_c = c1_l_a;
      c1_m_a = c1_jp1j;
      c1_i_b = c1_d_c - 1;
      c1_i84 = c1_m_a + c1_i_b;
      c1_n_a = c1_b_jp1j;
      c1_j_b = c1_i84;
      c1_o_a = c1_n_a;
      c1_k_b = c1_j_b;
      if (c1_o_a > c1_k_b) {
        c1_b_overflow = FALSE;
      } else {
        c1_b_overflow = (c1_k_b > 2147483646);
      }

      if (c1_b_overflow) {
        c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
      }

      for (c1_i = c1_b_jp1j; c1_i <= c1_i84; c1_i++) {
        c1_b_i = c1_i;
        c1_k_x = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_b_i), 1, 16, 1, 0) - 1];
        c1_e_y = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c1_jj), 1, 16, 1, 0) - 1];
        c1_z = c1_k_x / c1_e_y;
        c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_b_i), 1, 16, 1, 0) - 1] = c1_z;
      }
    } else {
      *c1_info = c1_b_j;
    }

    c1_l_b = c1_b_j;
    c1_e_c = 4 - c1_l_b;
    c1_p_a = c1_jj;
    c1_f_c = c1_p_a;
    c1_q_a = c1_jj;
    c1_g_c = c1_q_a;
    c1_m = c1_mmj;
    c1_e_n = c1_e_c;
    c1_g_ix0 = c1_jp1j;
    c1_d_iy0 = c1_f_c + 4;
    c1_ia0 = c1_g_c + 5;
    c1_d1 = -1.0;
    c1_b_eml_xger(chartInstance, c1_m, c1_e_n, c1_d1, c1_g_ix0, c1_d_iy0, c1_A,
                  c1_ia0);
  }

  if (*c1_info == 0) {
    if (!(c1_A[15] != 0.0)) {
      *c1_info = 4;
    }
  }
}

static void c1_b_eml_xger(SFc1_testInstanceStruct *chartInstance, int32_T c1_m,
  int32_T c1_n, real_T c1_alpha1, int32_T c1_ix0, int32_T c1_iy0, real_T c1_A[16],
  int32_T c1_ia0)
{
  int32_T c1_b_m;
  int32_T c1_b_n;
  real_T c1_b_alpha1;
  int32_T c1_b_ix0;
  int32_T c1_b_iy0;
  int32_T c1_b_ia0;
  int32_T c1_c_m;
  int32_T c1_c_n;
  real_T c1_c_alpha1;
  int32_T c1_c_ix0;
  int32_T c1_c_iy0;
  int32_T c1_c_ia0;
  int32_T c1_d_m;
  int32_T c1_d_n;
  real_T c1_d_alpha1;
  int32_T c1_d_ix0;
  int32_T c1_d_iy0;
  int32_T c1_d_ia0;
  int32_T c1_ixstart;
  int32_T c1_a;
  int32_T c1_jA;
  int32_T c1_jy;
  int32_T c1_e_n;
  int32_T c1_b;
  int32_T c1_b_b;
  boolean_T c1_overflow;
  int32_T c1_j;
  real_T c1_yjy;
  real_T c1_temp;
  int32_T c1_ix;
  int32_T c1_c_b;
  int32_T c1_i85;
  int32_T c1_b_a;
  int32_T c1_d_b;
  int32_T c1_i86;
  int32_T c1_c_a;
  int32_T c1_e_b;
  int32_T c1_d_a;
  int32_T c1_f_b;
  boolean_T c1_b_overflow;
  int32_T c1_ijA;
  int32_T c1_b_ijA;
  int32_T c1_e_a;
  int32_T c1_f_a;
  int32_T c1_g_a;
  c1_b_m = c1_m;
  c1_b_n = c1_n;
  c1_b_alpha1 = c1_alpha1;
  c1_b_ix0 = c1_ix0;
  c1_b_iy0 = c1_iy0;
  c1_b_ia0 = c1_ia0;
  c1_c_m = c1_b_m;
  c1_c_n = c1_b_n;
  c1_c_alpha1 = c1_b_alpha1;
  c1_c_ix0 = c1_b_ix0;
  c1_c_iy0 = c1_b_iy0;
  c1_c_ia0 = c1_b_ia0;
  c1_d_m = c1_c_m;
  c1_d_n = c1_c_n;
  c1_d_alpha1 = c1_c_alpha1;
  c1_d_ix0 = c1_c_ix0;
  c1_d_iy0 = c1_c_iy0;
  c1_d_ia0 = c1_c_ia0;
  if (c1_d_alpha1 == 0.0) {
  } else {
    c1_ixstart = c1_d_ix0;
    c1_a = c1_d_ia0 - 1;
    c1_jA = c1_a;
    c1_jy = c1_d_iy0;
    c1_e_n = c1_d_n;
    c1_b = c1_e_n;
    c1_b_b = c1_b;
    if (1 > c1_b_b) {
      c1_overflow = FALSE;
    } else {
      c1_overflow = (c1_b_b > 2147483646);
    }

    if (c1_overflow) {
      c1_check_forloop_overflow_error(chartInstance, c1_overflow);
    }

    for (c1_j = 1; c1_j <= c1_e_n; c1_j++) {
      c1_yjy = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c1_jy), 1, 16, 1, 0) - 1];
      if (c1_yjy != 0.0) {
        c1_temp = c1_yjy * c1_d_alpha1;
        c1_ix = c1_ixstart;
        c1_c_b = c1_jA + 1;
        c1_i85 = c1_c_b;
        c1_b_a = c1_d_m;
        c1_d_b = c1_jA;
        c1_i86 = c1_b_a + c1_d_b;
        c1_c_a = c1_i85;
        c1_e_b = c1_i86;
        c1_d_a = c1_c_a;
        c1_f_b = c1_e_b;
        if (c1_d_a > c1_f_b) {
          c1_b_overflow = FALSE;
        } else {
          c1_b_overflow = (c1_f_b > 2147483646);
        }

        if (c1_b_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_b_overflow);
        }

        for (c1_ijA = c1_i85; c1_ijA <= c1_i86; c1_ijA++) {
          c1_b_ijA = c1_ijA;
          c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ijA), 1, 16, 1, 0) - 1] =
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_b_ijA), 1, 16, 1, 0) - 1] +
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_ix), 1, 16, 1, 0) - 1] * c1_temp;
          c1_e_a = c1_ix + 1;
          c1_ix = c1_e_a;
        }
      }

      c1_f_a = c1_jy + 4;
      c1_jy = c1_f_a;
      c1_g_a = c1_jA + 4;
      c1_jA = c1_g_a;
    }
  }
}

static void c1_c_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A
  [16], real_T c1_B[24])
{
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_a;
  int32_T c1_c;
  int32_T c1_b;
  int32_T c1_b_c;
  int32_T c1_b_b;
  int32_T c1_jBcol;
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_b_a;
  int32_T c1_c_c;
  int32_T c1_c_b;
  int32_T c1_d_c;
  int32_T c1_d_b;
  int32_T c1_kAcol;
  int32_T c1_c_a;
  int32_T c1_e_b;
  int32_T c1_e_c;
  int32_T c1_d_a;
  int32_T c1_i87;
  boolean_T c1_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  int32_T c1_e_a;
  int32_T c1_f_b;
  int32_T c1_f_c;
  int32_T c1_f_a;
  int32_T c1_g_b;
  int32_T c1_g_c;
  int32_T c1_g_a;
  int32_T c1_h_b;
  int32_T c1_h_c;
  int32_T c1_h_a;
  int32_T c1_i_b;
  int32_T c1_i_c;
  c1_below_threshold(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_j = 1; c1_j < 7; c1_j++) {
    c1_b_j = c1_j;
    c1_a = c1_b_j;
    c1_c = c1_a;
    c1_b = c1_c - 1;
    c1_b_c = c1_b << 2;
    c1_b_b = c1_b_c;
    c1_jBcol = c1_b_b;
    for (c1_k = 1; c1_k < 5; c1_k++) {
      c1_b_k = c1_k;
      c1_b_a = c1_b_k;
      c1_c_c = c1_b_a;
      c1_c_b = c1_c_c - 1;
      c1_d_c = c1_c_b << 2;
      c1_d_b = c1_d_c;
      c1_kAcol = c1_d_b;
      c1_c_a = c1_b_k;
      c1_e_b = c1_jBcol;
      c1_e_c = c1_c_a + c1_e_b;
      if (c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_e_c), 1, 24, 1, 0) - 1] != 0.0) {
        c1_d_a = c1_b_k;
        c1_i87 = c1_d_a;
        c1_overflow = FALSE;
        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }

        for (c1_i = c1_i87 + 1; c1_i < 5; c1_i++) {
          c1_b_i = c1_i;
          c1_e_a = c1_b_i;
          c1_f_b = c1_jBcol;
          c1_f_c = c1_e_a + c1_f_b;
          c1_f_a = c1_b_i;
          c1_g_b = c1_jBcol;
          c1_g_c = c1_f_a + c1_g_b;
          c1_g_a = c1_b_k;
          c1_h_b = c1_jBcol;
          c1_h_c = c1_g_a + c1_h_b;
          c1_h_a = c1_b_i;
          c1_i_b = c1_kAcol;
          c1_i_c = c1_h_a + c1_i_b;
          c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_f_c), 1, 24, 1, 0) - 1] =
            c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_g_c), 1, 24, 1, 0) - 1] -
            c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_h_c), 1, 24, 1, 0) - 1] *
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_i_c), 1, 16, 1, 0) - 1];
        }
      }
    }
  }
}

static void c1_d_eml_xtrsm(SFc1_testInstanceStruct *chartInstance, real_T c1_A
  [16], real_T c1_B[24])
{
  int32_T c1_j;
  int32_T c1_b_j;
  int32_T c1_a;
  int32_T c1_c;
  int32_T c1_b;
  int32_T c1_b_c;
  int32_T c1_b_b;
  int32_T c1_jBcol;
  int32_T c1_k;
  int32_T c1_b_k;
  int32_T c1_b_a;
  int32_T c1_c_c;
  int32_T c1_c_b;
  int32_T c1_d_c;
  int32_T c1_d_b;
  int32_T c1_kAcol;
  int32_T c1_c_a;
  int32_T c1_e_b;
  int32_T c1_e_c;
  int32_T c1_d_a;
  int32_T c1_f_b;
  int32_T c1_f_c;
  int32_T c1_e_a;
  int32_T c1_g_b;
  int32_T c1_g_c;
  int32_T c1_f_a;
  int32_T c1_h_b;
  int32_T c1_h_c;
  real_T c1_x;
  real_T c1_y;
  real_T c1_z;
  int32_T c1_g_a;
  int32_T c1_i88;
  int32_T c1_i_b;
  int32_T c1_j_b;
  boolean_T c1_overflow;
  int32_T c1_i;
  int32_T c1_b_i;
  int32_T c1_h_a;
  int32_T c1_k_b;
  int32_T c1_i_c;
  int32_T c1_i_a;
  int32_T c1_l_b;
  int32_T c1_j_c;
  int32_T c1_j_a;
  int32_T c1_m_b;
  int32_T c1_k_c;
  int32_T c1_k_a;
  int32_T c1_n_b;
  int32_T c1_l_c;
  c1_below_threshold(chartInstance);
  c1_b_eml_scalar_eg(chartInstance);
  for (c1_j = 1; c1_j < 7; c1_j++) {
    c1_b_j = c1_j;
    c1_a = c1_b_j;
    c1_c = c1_a;
    c1_b = c1_c - 1;
    c1_b_c = c1_b << 2;
    c1_b_b = c1_b_c;
    c1_jBcol = c1_b_b;
    for (c1_k = 4; c1_k > 0; c1_k--) {
      c1_b_k = c1_k;
      c1_b_a = c1_b_k;
      c1_c_c = c1_b_a;
      c1_c_b = c1_c_c - 1;
      c1_d_c = c1_c_b << 2;
      c1_d_b = c1_d_c;
      c1_kAcol = c1_d_b;
      c1_c_a = c1_b_k;
      c1_e_b = c1_jBcol;
      c1_e_c = c1_c_a + c1_e_b;
      if (c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_e_c), 1, 24, 1, 0) - 1] != 0.0) {
        c1_d_a = c1_b_k;
        c1_f_b = c1_jBcol;
        c1_f_c = c1_d_a + c1_f_b;
        c1_e_a = c1_b_k;
        c1_g_b = c1_jBcol;
        c1_g_c = c1_e_a + c1_g_b;
        c1_f_a = c1_b_k;
        c1_h_b = c1_kAcol;
        c1_h_c = c1_f_a + c1_h_b;
        c1_x = c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c1_g_c), 1, 24, 1, 0) - 1];
        c1_y = c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c1_h_c), 1, 16, 1, 0) - 1];
        c1_z = c1_x / c1_y;
        c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c1_f_c), 1, 24, 1, 0) - 1] = c1_z;
        c1_g_a = c1_b_k - 1;
        c1_i88 = c1_g_a;
        c1_i_b = c1_i88;
        c1_j_b = c1_i_b;
        if (1 > c1_j_b) {
          c1_overflow = FALSE;
        } else {
          c1_overflow = (c1_j_b > 2147483646);
        }

        if (c1_overflow) {
          c1_check_forloop_overflow_error(chartInstance, c1_overflow);
        }

        for (c1_i = 1; c1_i <= c1_i88; c1_i++) {
          c1_b_i = c1_i;
          c1_h_a = c1_b_i;
          c1_k_b = c1_jBcol;
          c1_i_c = c1_h_a + c1_k_b;
          c1_i_a = c1_b_i;
          c1_l_b = c1_jBcol;
          c1_j_c = c1_i_a + c1_l_b;
          c1_j_a = c1_b_k;
          c1_m_b = c1_jBcol;
          c1_k_c = c1_j_a + c1_m_b;
          c1_k_a = c1_b_i;
          c1_n_b = c1_kAcol;
          c1_l_c = c1_k_a + c1_n_b;
          c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_i_c), 1, 24, 1, 0) - 1] =
            c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_j_c), 1, 24, 1, 0) - 1] -
            c1_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_k_c), 1, 24, 1, 0) - 1] *
            c1_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c1_l_c), 1, 16, 1, 0) - 1];
        }
      }
    }
  }
}

static void init_dsm_address_info(SFc1_testInstanceStruct *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c1_test_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(106092327U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(102684612U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1708578213U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1679719394U);
}

mxArray *sf_c1_test_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("trNaZf6lUT0VkqbCbrSFSB");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c1_test_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c1_test(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"omega\",},{M[8],M[0],T\"is_active_c1_test\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_test_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_testInstanceStruct *chartInstance;
    chartInstance = (SFc1_testInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _testMachineNumber_,
           1,
           1,
           1,
           5,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_testMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_testMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _testMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"u");
          _SFD_SET_DATA_PROPS(1,2,0,1,"omega");
          _SFD_SET_DATA_PROPS(2,1,1,0,"k");
          _SFD_SET_DATA_PROPS(3,1,1,0,"l");
          _SFD_SET_DATA_PROPS(4,1,1,0,"b");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,266);
        _SFD_TRANS_COV_WTS(0,0,0,1,0);
        if (chartAlreadyPresent==0) {
          _SFD_TRANS_COV_MAPS(0,
                              0,NULL,NULL,
                              0,NULL,NULL,
                              1,NULL,NULL,
                              0,NULL,NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)
            c1_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c1_k;
          real_T *c1_l;
          real_T *c1_b;
          real_T (*c1_u)[4];
          real_T (*c1_omega)[6];
          c1_b = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c1_l = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c1_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c1_omega = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
          c1_u = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c1_u);
          _SFD_SET_DATA_VALUE_PTR(1U, *c1_omega);
          _SFD_SET_DATA_VALUE_PTR(2U, c1_k);
          _SFD_SET_DATA_VALUE_PTR(3U, c1_l);
          _SFD_SET_DATA_VALUE_PTR(4U, c1_b);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _testMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "mWYR3yrLHULIRIwlwvdjd";
}

static void sf_opaque_initialize_c1_test(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc1_testInstanceStruct*) chartInstanceVar)->S,0);
  initialize_params_c1_test((SFc1_testInstanceStruct*) chartInstanceVar);
  initialize_c1_test((SFc1_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c1_test(void *chartInstanceVar)
{
  enable_c1_test((SFc1_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_test(void *chartInstanceVar)
{
  disable_c1_test((SFc1_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c1_test(void *chartInstanceVar)
{
  sf_c1_test((SFc1_testInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c1_test(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c1_test((SFc1_testInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_test();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c1_test(SimStruct* S, const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_test();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c1_test((SFc1_testInstanceStruct*)chartInfo->chartInstance,
                        mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c1_test(SimStruct* S)
{
  return sf_internal_get_sim_state_c1_test(S);
}

static void sf_opaque_set_sim_state_c1_test(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c1_test(S, st);
}

static void sf_opaque_terminate_c1_test(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_testInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_optimization_info();
    }

    finalize_c1_test((SFc1_testInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_test((SFc1_testInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_test(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c1_test((SFc1_testInstanceStruct*)(((ChartInfoStruct *)
      ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c1_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,1,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,1);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 4; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(404780754U));
  ssSetChecksum1(S,(1082443012U));
  ssSetChecksum2(S,(4270730837U));
  ssSetChecksum3(S,(365230326U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_test(SimStruct *S)
{
  SFc1_testInstanceStruct *chartInstance;
  chartInstance = (SFc1_testInstanceStruct *)utMalloc(sizeof
    (SFc1_testInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_testInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c1_test;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c1_test;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c1_test;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_test;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_test;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c1_test;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c1_test;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c1_test;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_test;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_test;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_test;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c1_test_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_test(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_test(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_test(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_test_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
