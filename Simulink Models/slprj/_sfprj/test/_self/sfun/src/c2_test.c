/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_sfun.h"
#include "c2_test.h"
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
static const char * c2_debug_family_names[20] = { "omega_squared", "M_rb",
  "C_rb", "RotWind", "wind_b", "tau_b", "nargin", "nargout", "nu", "eta",
  "omega", "g", "m", "I", "b", "k", "l", "wind", "windDir", "nu_dot" };

static const char * c2_b_debug_family_names[4] = { "nargin", "nargout", "r", "S"
};

static const char * c2_c_debug_family_names[14] = { "nu1", "nu2", "Msym", "M11",
  "M12", "M21", "M22", "dt_dnu1", "dt_dnu2", "nargin", "nargout", "M", "nu", "C"
};

/* Function Declarations */
static void initialize_c2_test(SFc2_testInstanceStruct *chartInstance);
static void initialize_params_c2_test(SFc2_testInstanceStruct *chartInstance);
static void enable_c2_test(SFc2_testInstanceStruct *chartInstance);
static void disable_c2_test(SFc2_testInstanceStruct *chartInstance);
static void c2_update_debugger_state_c2_test(SFc2_testInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c2_test(SFc2_testInstanceStruct
  *chartInstance);
static void set_sim_state_c2_test(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_st);
static void finalize_c2_test(SFc2_testInstanceStruct *chartInstance);
static void sf_c2_test(SFc2_testInstanceStruct *chartInstance);
static void c2_chartstep_c2_test(SFc2_testInstanceStruct *chartInstance);
static void initSimStructsc2_test(SFc2_testInstanceStruct *chartInstance);
static void registerMessagesc2_test(SFc2_testInstanceStruct *chartInstance);
static void c2_m2c(SFc2_testInstanceStruct *chartInstance, real_T c2_M[36],
                   real_T c2_nu[6], real_T c2_C[36]);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static void c2_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_nu_dot, const char_T *c2_identifier, real_T c2_y[6]);
static void c2_b_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[6]);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static real_T c2_c_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_d_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[2]);
static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_e_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[4]);
static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_f_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[36]);
static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_g_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9]);
static void c2_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static void c2_h_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[3]);
static void c2_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static void c2_info_helper(c2_ResolvedFunctionInfo c2_info[159]);
static void c2_b_info_helper(c2_ResolvedFunctionInfo c2_info[159]);
static void c2_c_info_helper(c2_ResolvedFunctionInfo c2_info[159]);
static void c2_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance);
static void c2_b_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance);
static real_T c2_sum(SFc2_testInstanceStruct *chartInstance, real_T c2_x[6]);
static void c2_c_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance);
static void c2_d_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance);
static void c2_realmin(SFc2_testInstanceStruct *chartInstance);
static void c2_eps(SFc2_testInstanceStruct *chartInstance);
static void c2_eml_matlab_zgetrf(SFc2_testInstanceStruct *chartInstance, real_T
  c2_A[36], real_T c2_b_A[36], int32_T c2_ipiv[6], int32_T *c2_info);
static void c2_check_forloop_overflow_error(SFc2_testInstanceStruct
  *chartInstance, boolean_T c2_overflow);
static void c2_eml_xger(SFc2_testInstanceStruct *chartInstance, int32_T c2_m,
  int32_T c2_n, real_T c2_alpha1, int32_T c2_ix0, int32_T c2_iy0, real_T c2_A[36],
  int32_T c2_ia0, real_T c2_b_A[36]);
static void c2_eml_warning(SFc2_testInstanceStruct *chartInstance);
static void c2_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A[36],
  real_T c2_B[6], real_T c2_b_B[6]);
static void c2_below_threshold(SFc2_testInstanceStruct *chartInstance);
static void c2_e_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance);
static void c2_b_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A
  [36], real_T c2_B[6], real_T c2_b_B[6]);
static const mxArray *c2_h_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_i_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_j_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_b_is_active_c2_test, const char_T *c2_identifier);
static uint8_T c2_k_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_eml_matlab_zgetrf(SFc2_testInstanceStruct *chartInstance,
  real_T c2_A[36], int32_T c2_ipiv[6], int32_T *c2_info);
static void c2_b_eml_xger(SFc2_testInstanceStruct *chartInstance, int32_T c2_m,
  int32_T c2_n, real_T c2_alpha1, int32_T c2_ix0, int32_T c2_iy0, real_T c2_A[36],
  int32_T c2_ia0);
static void c2_c_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A
  [36], real_T c2_B[6]);
static void c2_d_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A
  [36], real_T c2_B[6]);
static void init_dsm_address_info(SFc2_testInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c2_test(SFc2_testInstanceStruct *chartInstance)
{
  chartInstance->c2_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c2_is_active_c2_test = 0U;
}

static void initialize_params_c2_test(SFc2_testInstanceStruct *chartInstance)
{
}

static void enable_c2_test(SFc2_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c2_test(SFc2_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c2_update_debugger_state_c2_test(SFc2_testInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c2_test(SFc2_testInstanceStruct
  *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  int32_T c2_i0;
  real_T c2_u[6];
  const mxArray *c2_b_y = NULL;
  uint8_T c2_hoistedGlobal;
  uint8_T c2_b_u;
  const mxArray *c2_c_y = NULL;
  real_T (*c2_nu_dot)[6];
  c2_nu_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c2_st = NULL;
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellarray(2), FALSE);
  for (c2_i0 = 0; c2_i0 < 6; c2_i0++) {
    c2_u[c2_i0] = (*c2_nu_dot)[c2_i0];
  }

  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_hoistedGlobal = chartInstance->c2_is_active_c2_test;
  c2_b_u = c2_hoistedGlobal;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  sf_mex_assign(&c2_st, c2_y, FALSE);
  return c2_st;
}

static void set_sim_state_c2_test(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_st)
{
  const mxArray *c2_u;
  real_T c2_dv0[6];
  int32_T c2_i1;
  real_T (*c2_nu_dot)[6];
  c2_nu_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c2_doneDoubleBufferReInit = TRUE;
  c2_u = sf_mex_dup(c2_st);
  c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 0)),
                      "nu_dot", c2_dv0);
  for (c2_i1 = 0; c2_i1 < 6; c2_i1++) {
    (*c2_nu_dot)[c2_i1] = c2_dv0[c2_i1];
  }

  chartInstance->c2_is_active_c2_test = c2_j_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c2_u, 1)), "is_active_c2_test");
  sf_mex_destroy(&c2_u);
  c2_update_debugger_state_c2_test(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_test(SFc2_testInstanceStruct *chartInstance)
{
}

static void sf_c2_test(SFc2_testInstanceStruct *chartInstance)
{
  int32_T c2_i2;
  int32_T c2_i3;
  int32_T c2_i4;
  int32_T c2_i5;
  int32_T c2_i6;
  int32_T c2_i7;
  real_T *c2_g;
  real_T *c2_m;
  real_T *c2_b;
  real_T *c2_k;
  real_T *c2_l;
  real_T *c2_windDir;
  real_T (*c2_wind)[2];
  real_T (*c2_I)[9];
  real_T (*c2_omega)[6];
  real_T (*c2_eta)[6];
  real_T (*c2_nu_dot)[6];
  real_T (*c2_nu)[6];
  c2_windDir = (real_T *)ssGetInputPortSignal(chartInstance->S, 10);
  c2_wind = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 9);
  c2_l = (real_T *)ssGetInputPortSignal(chartInstance->S, 8);
  c2_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 7);
  c2_b = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c2_I = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 5);
  c2_m = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c2_g = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c2_omega = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 2);
  c2_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c2_nu_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c2_nu = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  for (c2_i2 = 0; c2_i2 < 6; c2_i2++) {
    _SFD_DATA_RANGE_CHECK((*c2_nu)[c2_i2], 0U);
  }

  for (c2_i3 = 0; c2_i3 < 6; c2_i3++) {
    _SFD_DATA_RANGE_CHECK((*c2_nu_dot)[c2_i3], 1U);
  }

  for (c2_i4 = 0; c2_i4 < 6; c2_i4++) {
    _SFD_DATA_RANGE_CHECK((*c2_eta)[c2_i4], 2U);
  }

  for (c2_i5 = 0; c2_i5 < 6; c2_i5++) {
    _SFD_DATA_RANGE_CHECK((*c2_omega)[c2_i5], 3U);
  }

  _SFD_DATA_RANGE_CHECK(*c2_g, 4U);
  _SFD_DATA_RANGE_CHECK(*c2_m, 5U);
  for (c2_i6 = 0; c2_i6 < 9; c2_i6++) {
    _SFD_DATA_RANGE_CHECK((*c2_I)[c2_i6], 6U);
  }

  _SFD_DATA_RANGE_CHECK(*c2_b, 7U);
  _SFD_DATA_RANGE_CHECK(*c2_k, 8U);
  _SFD_DATA_RANGE_CHECK(*c2_l, 9U);
  for (c2_i7 = 0; c2_i7 < 2; c2_i7++) {
    _SFD_DATA_RANGE_CHECK((*c2_wind)[c2_i7], 10U);
  }

  _SFD_DATA_RANGE_CHECK(*c2_windDir, 11U);
  chartInstance->c2_sfEvent = CALL_EVENT;
  c2_chartstep_c2_test(chartInstance);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_testMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c2_chartstep_c2_test(SFc2_testInstanceStruct *chartInstance)
{
  real_T c2_hoistedGlobal;
  real_T c2_b_hoistedGlobal;
  real_T c2_c_hoistedGlobal;
  real_T c2_d_hoistedGlobal;
  real_T c2_e_hoistedGlobal;
  real_T c2_f_hoistedGlobal;
  int32_T c2_i8;
  real_T c2_nu[6];
  int32_T c2_i9;
  real_T c2_eta[6];
  int32_T c2_i10;
  real_T c2_omega[6];
  real_T c2_g;
  real_T c2_m;
  int32_T c2_i11;
  real_T c2_I[9];
  real_T c2_b;
  real_T c2_k;
  real_T c2_l;
  int32_T c2_i12;
  real_T c2_wind[2];
  real_T c2_windDir;
  uint32_T c2_debug_family_var_map[20];
  real_T c2_omega_squared[6];
  real_T c2_M_rb[36];
  real_T c2_C_rb[36];
  real_T c2_RotWind[4];
  real_T c2_wind_b[2];
  real_T c2_tau_b[6];
  real_T c2_nargin = 11.0;
  real_T c2_nargout = 1.0;
  real_T c2_nu_dot[6];
  int32_T c2_i13;
  real_T c2_a;
  int32_T c2_i14;
  static real_T c2_b_b[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

  real_T c2_y[9];
  int32_T c2_i15;
  int32_T c2_i16;
  int32_T c2_i17;
  int32_T c2_i18;
  int32_T c2_i19;
  int32_T c2_i20;
  int32_T c2_i21;
  int32_T c2_i22;
  int32_T c2_i23;
  int32_T c2_i24;
  int32_T c2_i25;
  int32_T c2_i26;
  int32_T c2_i27;
  int32_T c2_i28;
  int32_T c2_i29;
  real_T c2_b_M_rb[36];
  int32_T c2_i30;
  real_T c2_b_nu[6];
  real_T c2_dv1[36];
  int32_T c2_i31;
  real_T c2_x;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_d_x;
  real_T c2_e_x;
  real_T c2_f_x;
  real_T c2_g_x;
  real_T c2_h_x;
  int32_T c2_i32;
  real_T c2_b_a[4];
  int32_T c2_i33;
  real_T c2_c_b[2];
  int32_T c2_i34;
  int32_T c2_i35;
  int32_T c2_i36;
  real_T c2_C[2];
  int32_T c2_i37;
  int32_T c2_i38;
  int32_T c2_i39;
  int32_T c2_i40;
  int32_T c2_i41;
  int32_T c2_i42;
  real_T c2_c_a;
  real_T c2_d_b;
  real_T c2_b_y;
  real_T c2_i_x;
  real_T c2_j_x;
  real_T c2_d_a;
  real_T c2_e_b;
  real_T c2_c_y;
  real_T c2_e_a;
  real_T c2_f_b;
  real_T c2_d_y;
  real_T c2_k_x;
  real_T c2_l_x;
  real_T c2_f_a;
  real_T c2_g_b;
  real_T c2_e_y;
  real_T c2_m_x;
  real_T c2_n_x;
  real_T c2_g_a;
  real_T c2_h_b;
  real_T c2_f_y;
  real_T c2_h_a;
  real_T c2_i_b;
  real_T c2_g_y;
  real_T c2_o_x;
  real_T c2_p_x;
  real_T c2_i_a;
  real_T c2_j_b;
  real_T c2_h_y;
  real_T c2_q_x;
  real_T c2_r_x;
  real_T c2_j_a;
  real_T c2_k_b;
  real_T c2_i_y;
  real_T c2_k_a;
  int32_T c2_i43;
  real_T c2_b_omega_squared[6];
  real_T c2_l_b;
  real_T c2_j_y;
  real_T c2_l_a;
  real_T c2_m_b;
  real_T c2_k_y;
  real_T c2_m_a;
  int32_T c2_i44;
  static real_T c2_n_b[6] = { -0.5, 0.5, 1.0, 0.5, -0.5, -1.0 };

  real_T c2_l_y[6];
  int32_T c2_i45;
  real_T c2_o_b[6];
  real_T c2_m_y;
  int32_T c2_b_k;
  int32_T c2_c_k;
  real_T c2_p_b;
  real_T c2_n_y;
  real_T c2_n_a;
  real_T c2_q_b;
  real_T c2_o_y;
  real_T c2_o_a;
  int32_T c2_i46;
  static real_T c2_r_b[6] = { 1.0, 1.0, 0.0, -1.0, -1.0, 0.0 };

  int32_T c2_i47;
  real_T c2_p_y;
  int32_T c2_d_k;
  int32_T c2_e_k;
  real_T c2_p_a;
  int32_T c2_i48;
  static real_T c2_s_b[6] = { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };

  int32_T c2_i49;
  real_T c2_q_y;
  int32_T c2_f_k;
  int32_T c2_g_k;
  int32_T c2_i50;
  real_T c2_q_a[36];
  int32_T c2_i51;
  int32_T c2_i52;
  real_T c2_r_y[6];
  int32_T c2_i53;
  int32_T c2_i54;
  int32_T c2_i55;
  int32_T c2_i56;
  int32_T c2_info;
  int32_T c2_ipiv[6];
  int32_T c2_b_info;
  int32_T c2_c_info;
  int32_T c2_d_info;
  int32_T c2_i57;
  int32_T c2_i;
  int32_T c2_b_i;
  int32_T c2_ip;
  real_T c2_temp;
  int32_T c2_i58;
  real_T c2_dv2[36];
  int32_T c2_i59;
  real_T c2_dv3[36];
  int32_T c2_i60;
  real_T c2_dv4[36];
  int32_T c2_i61;
  real_T c2_dv5[36];
  int32_T c2_i62;
  real_T (*c2_b_nu_dot)[6];
  real_T *c2_b_windDir;
  real_T *c2_b_l;
  real_T *c2_h_k;
  real_T *c2_t_b;
  real_T *c2_b_m;
  real_T *c2_b_g;
  real_T (*c2_b_wind)[2];
  real_T (*c2_b_I)[9];
  real_T (*c2_b_omega)[6];
  real_T (*c2_b_eta)[6];
  real_T (*c2_c_nu)[6];
  c2_b_windDir = (real_T *)ssGetInputPortSignal(chartInstance->S, 10);
  c2_b_wind = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 9);
  c2_b_l = (real_T *)ssGetInputPortSignal(chartInstance->S, 8);
  c2_h_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 7);
  c2_t_b = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c2_b_I = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 5);
  c2_b_m = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c2_b_g = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c2_b_omega = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 2);
  c2_b_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c2_b_nu_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c2_c_nu = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
  c2_hoistedGlobal = *c2_b_g;
  c2_b_hoistedGlobal = *c2_b_m;
  c2_c_hoistedGlobal = *c2_t_b;
  c2_d_hoistedGlobal = *c2_h_k;
  c2_e_hoistedGlobal = *c2_b_l;
  c2_f_hoistedGlobal = *c2_b_windDir;
  for (c2_i8 = 0; c2_i8 < 6; c2_i8++) {
    c2_nu[c2_i8] = (*c2_c_nu)[c2_i8];
  }

  for (c2_i9 = 0; c2_i9 < 6; c2_i9++) {
    c2_eta[c2_i9] = (*c2_b_eta)[c2_i9];
  }

  for (c2_i10 = 0; c2_i10 < 6; c2_i10++) {
    c2_omega[c2_i10] = (*c2_b_omega)[c2_i10];
  }

  c2_g = c2_hoistedGlobal;
  c2_m = c2_b_hoistedGlobal;
  for (c2_i11 = 0; c2_i11 < 9; c2_i11++) {
    c2_I[c2_i11] = (*c2_b_I)[c2_i11];
  }

  c2_b = c2_c_hoistedGlobal;
  c2_k = c2_d_hoistedGlobal;
  c2_l = c2_e_hoistedGlobal;
  for (c2_i12 = 0; c2_i12 < 2; c2_i12++) {
    c2_wind[c2_i12] = (*c2_b_wind)[c2_i12];
  }

  c2_windDir = c2_f_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 20U, 20U, c2_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_omega_squared, 0U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_M_rb, 1U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_C_rb, 2U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_RotWind, 3U, c2_e_sf_marshallOut,
    c2_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_wind_b, 4U, c2_c_sf_marshallOut,
    c2_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_tau_b, 5U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 6U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 7U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_nu, 8U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_eta, 9U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_omega, 10U, c2_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_g, 11U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_m, 12U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_I, 13U, c2_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_b, 14U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_k, 15U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_l, 16U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c2_wind, 17U, c2_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c2_windDir, 18U, c2_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_nu_dot, 19U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 4);
  for (c2_i13 = 0; c2_i13 < 6; c2_i13++) {
    c2_omega_squared[c2_i13] = c2_omega[c2_i13];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 6);
  c2_a = c2_m;
  for (c2_i14 = 0; c2_i14 < 9; c2_i14++) {
    c2_y[c2_i14] = c2_a * c2_b_b[c2_i14];
  }

  c2_i15 = 0;
  c2_i16 = 0;
  for (c2_i17 = 0; c2_i17 < 3; c2_i17++) {
    for (c2_i18 = 0; c2_i18 < 3; c2_i18++) {
      c2_M_rb[c2_i18 + c2_i15] = c2_y[c2_i18 + c2_i16];
    }

    c2_i15 += 6;
    c2_i16 += 3;
  }

  c2_i19 = 0;
  for (c2_i20 = 0; c2_i20 < 3; c2_i20++) {
    for (c2_i21 = 0; c2_i21 < 3; c2_i21++) {
      c2_M_rb[(c2_i21 + c2_i19) + 18] = 0.0;
    }

    c2_i19 += 6;
  }

  c2_i22 = 0;
  for (c2_i23 = 0; c2_i23 < 3; c2_i23++) {
    for (c2_i24 = 0; c2_i24 < 3; c2_i24++) {
      c2_M_rb[(c2_i24 + c2_i22) + 3] = 0.0;
    }

    c2_i22 += 6;
  }

  c2_i25 = 0;
  c2_i26 = 0;
  for (c2_i27 = 0; c2_i27 < 3; c2_i27++) {
    for (c2_i28 = 0; c2_i28 < 3; c2_i28++) {
      c2_M_rb[(c2_i28 + c2_i25) + 21] = c2_I[c2_i28 + c2_i26];
    }

    c2_i25 += 6;
    c2_i26 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 12);
  for (c2_i29 = 0; c2_i29 < 36; c2_i29++) {
    c2_b_M_rb[c2_i29] = c2_M_rb[c2_i29];
  }

  for (c2_i30 = 0; c2_i30 < 6; c2_i30++) {
    c2_b_nu[c2_i30] = c2_nu[c2_i30];
  }

  c2_m2c(chartInstance, c2_b_M_rb, c2_b_nu, c2_dv1);
  for (c2_i31 = 0; c2_i31 < 36; c2_i31++) {
    c2_C_rb[c2_i31] = c2_dv1[c2_i31];
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 14);
  c2_x = c2_windDir;
  c2_b_x = c2_x;
  c2_b_x = muDoubleScalarCos(c2_b_x);
  c2_c_x = c2_windDir;
  c2_d_x = c2_c_x;
  c2_d_x = muDoubleScalarSin(c2_d_x);
  c2_e_x = c2_windDir;
  c2_f_x = c2_e_x;
  c2_f_x = muDoubleScalarSin(c2_f_x);
  c2_g_x = c2_windDir;
  c2_h_x = c2_g_x;
  c2_h_x = muDoubleScalarCos(c2_h_x);
  c2_RotWind[0] = c2_b_x;
  c2_RotWind[2] = c2_d_x;
  c2_RotWind[1] = -c2_f_x;
  c2_RotWind[3] = c2_h_x;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 16);
  for (c2_i32 = 0; c2_i32 < 4; c2_i32++) {
    c2_b_a[c2_i32] = c2_RotWind[c2_i32];
  }

  for (c2_i33 = 0; c2_i33 < 2; c2_i33++) {
    c2_c_b[c2_i33] = c2_wind[c2_i33];
  }

  c2_b_eml_scalar_eg(chartInstance);
  c2_b_eml_scalar_eg(chartInstance);
  for (c2_i34 = 0; c2_i34 < 2; c2_i34++) {
    c2_wind_b[c2_i34] = 0.0;
  }

  for (c2_i35 = 0; c2_i35 < 2; c2_i35++) {
    c2_wind_b[c2_i35] = 0.0;
  }

  for (c2_i36 = 0; c2_i36 < 2; c2_i36++) {
    c2_C[c2_i36] = c2_wind_b[c2_i36];
  }

  for (c2_i37 = 0; c2_i37 < 2; c2_i37++) {
    c2_wind_b[c2_i37] = c2_C[c2_i37];
  }

  for (c2_i38 = 0; c2_i38 < 2; c2_i38++) {
    c2_C[c2_i38] = c2_wind_b[c2_i38];
  }

  for (c2_i39 = 0; c2_i39 < 2; c2_i39++) {
    c2_wind_b[c2_i39] = c2_C[c2_i39];
  }

  for (c2_i40 = 0; c2_i40 < 2; c2_i40++) {
    c2_wind_b[c2_i40] = 0.0;
    c2_i41 = 0;
    for (c2_i42 = 0; c2_i42 < 2; c2_i42++) {
      c2_wind_b[c2_i40] += c2_b_a[c2_i41 + c2_i40] * c2_c_b[c2_i42];
      c2_i41 += 2;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 18);
  c2_c_a = -c2_m;
  c2_d_b = c2_g;
  c2_b_y = c2_c_a * c2_d_b;
  c2_i_x = c2_eta[4];
  c2_j_x = c2_i_x;
  c2_j_x = muDoubleScalarSin(c2_j_x);
  c2_d_a = c2_b_y;
  c2_e_b = c2_j_x;
  c2_c_y = c2_d_a * c2_e_b;
  c2_e_a = c2_m;
  c2_f_b = c2_g;
  c2_d_y = c2_e_a * c2_f_b;
  c2_k_x = c2_eta[4];
  c2_l_x = c2_k_x;
  c2_l_x = muDoubleScalarCos(c2_l_x);
  c2_f_a = c2_d_y;
  c2_g_b = c2_l_x;
  c2_e_y = c2_f_a * c2_g_b;
  c2_m_x = c2_eta[3];
  c2_n_x = c2_m_x;
  c2_n_x = muDoubleScalarSin(c2_n_x);
  c2_g_a = c2_e_y;
  c2_h_b = c2_n_x;
  c2_f_y = c2_g_a * c2_h_b;
  c2_h_a = c2_m;
  c2_i_b = c2_g;
  c2_g_y = c2_h_a * c2_i_b;
  c2_o_x = c2_eta[4];
  c2_p_x = c2_o_x;
  c2_p_x = muDoubleScalarCos(c2_p_x);
  c2_i_a = c2_g_y;
  c2_j_b = c2_p_x;
  c2_h_y = c2_i_a * c2_j_b;
  c2_q_x = c2_eta[3];
  c2_r_x = c2_q_x;
  c2_r_x = muDoubleScalarCos(c2_r_x);
  c2_j_a = c2_h_y;
  c2_k_b = c2_r_x;
  c2_i_y = c2_j_a * c2_k_b;
  c2_k_a = c2_k;
  for (c2_i43 = 0; c2_i43 < 6; c2_i43++) {
    c2_b_omega_squared[c2_i43] = c2_omega_squared[c2_i43];
  }

  c2_l_b = c2_sum(chartInstance, c2_b_omega_squared);
  c2_j_y = c2_k_a * c2_l_b;
  c2_l_a = c2_k;
  c2_m_b = c2_l;
  c2_k_y = c2_l_a * c2_m_b;
  c2_m_a = c2_k_y;
  for (c2_i44 = 0; c2_i44 < 6; c2_i44++) {
    c2_l_y[c2_i44] = c2_m_a * c2_n_b[c2_i44];
  }

  for (c2_i45 = 0; c2_i45 < 6; c2_i45++) {
    c2_o_b[c2_i45] = c2_omega_squared[c2_i45];
  }

  c2_c_eml_scalar_eg(chartInstance);
  c2_c_eml_scalar_eg(chartInstance);
  c2_m_y = 0.0;
  for (c2_b_k = 1; c2_b_k < 7; c2_b_k++) {
    c2_c_k = c2_b_k;
    c2_m_y += c2_l_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c2_c_k), 1, 6, 1, 0) - 1] *
      c2_o_b[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c2_c_k), 1, 6, 1, 0) - 1];
  }

  c2_p_b = c2_k;
  c2_n_y = 0.8660254037844386 * c2_p_b;
  c2_n_a = c2_n_y;
  c2_q_b = c2_l;
  c2_o_y = c2_n_a * c2_q_b;
  c2_o_a = c2_o_y;
  for (c2_i46 = 0; c2_i46 < 6; c2_i46++) {
    c2_l_y[c2_i46] = c2_o_a * c2_r_b[c2_i46];
  }

  for (c2_i47 = 0; c2_i47 < 6; c2_i47++) {
    c2_o_b[c2_i47] = c2_omega_squared[c2_i47];
  }

  c2_c_eml_scalar_eg(chartInstance);
  c2_c_eml_scalar_eg(chartInstance);
  c2_p_y = 0.0;
  for (c2_d_k = 1; c2_d_k < 7; c2_d_k++) {
    c2_e_k = c2_d_k;
    c2_p_y += c2_l_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c2_e_k), 1, 6, 1, 0) - 1] *
      c2_o_b[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c2_e_k), 1, 6, 1, 0) - 1];
  }

  c2_p_a = c2_b;
  for (c2_i48 = 0; c2_i48 < 6; c2_i48++) {
    c2_l_y[c2_i48] = c2_p_a * c2_s_b[c2_i48];
  }

  for (c2_i49 = 0; c2_i49 < 6; c2_i49++) {
    c2_o_b[c2_i49] = c2_omega_squared[c2_i49];
  }

  c2_c_eml_scalar_eg(chartInstance);
  c2_c_eml_scalar_eg(chartInstance);
  c2_q_y = 0.0;
  for (c2_f_k = 1; c2_f_k < 7; c2_f_k++) {
    c2_g_k = c2_f_k;
    c2_q_y += c2_l_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK
      ("", (real_T)c2_g_k), 1, 6, 1, 0) - 1] *
      c2_o_b[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c2_g_k), 1, 6, 1, 0) - 1];
  }

  c2_tau_b[0] = c2_c_y + c2_wind_b[0];
  c2_tau_b[1] = c2_f_y + c2_wind_b[1];
  c2_tau_b[2] = c2_i_y - c2_j_y;
  c2_tau_b[3] = c2_m_y;
  c2_tau_b[4] = c2_p_y;
  c2_tau_b[5] = c2_q_y;
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, 25);
  for (c2_i50 = 0; c2_i50 < 36; c2_i50++) {
    c2_q_a[c2_i50] = c2_C_rb[c2_i50];
  }

  for (c2_i51 = 0; c2_i51 < 6; c2_i51++) {
    c2_o_b[c2_i51] = c2_nu[c2_i51];
  }

  c2_d_eml_scalar_eg(chartInstance);
  c2_d_eml_scalar_eg(chartInstance);
  for (c2_i52 = 0; c2_i52 < 6; c2_i52++) {
    c2_r_y[c2_i52] = 0.0;
    c2_i53 = 0;
    for (c2_i54 = 0; c2_i54 < 6; c2_i54++) {
      c2_r_y[c2_i52] += c2_q_a[c2_i53 + c2_i52] * c2_o_b[c2_i54];
      c2_i53 += 6;
    }
  }

  for (c2_i55 = 0; c2_i55 < 36; c2_i55++) {
    c2_q_a[c2_i55] = c2_M_rb[c2_i55];
  }

  for (c2_i56 = 0; c2_i56 < 6; c2_i56++) {
    c2_r_y[c2_i56] = c2_tau_b[c2_i56] - c2_r_y[c2_i56];
  }

  c2_b_eml_matlab_zgetrf(chartInstance, c2_q_a, c2_ipiv, &c2_info);
  c2_b_info = c2_info;
  c2_c_info = c2_b_info;
  c2_d_info = c2_c_info;
  if (c2_d_info > 0) {
    c2_eml_warning(chartInstance);
  }

  c2_d_eml_scalar_eg(chartInstance);
  for (c2_i57 = 0; c2_i57 < 6; c2_i57++) {
    c2_nu_dot[c2_i57] = c2_r_y[c2_i57];
  }

  for (c2_i = 1; c2_i < 7; c2_i++) {
    c2_b_i = c2_i;
    if (c2_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_b_i), 1, 6, 1, 0) - 1] != c2_b_i) {
      c2_ip = c2_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c2_b_i), 1, 6, 1, 0) - 1];
      c2_temp = c2_nu_dot[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
        _SFD_INTEGER_CHECK("", (real_T)c2_b_i), 1, 6, 1, 0) - 1];
      c2_nu_dot[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c2_b_i), 1, 6, 1, 0) - 1] =
        c2_nu_dot[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c2_ip), 1, 6, 1, 0) - 1];
      c2_nu_dot[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c2_ip), 1, 6, 1, 0) - 1] = c2_temp;
    }
  }

  for (c2_i58 = 0; c2_i58 < 36; c2_i58++) {
    c2_dv2[c2_i58] = c2_q_a[c2_i58];
  }

  for (c2_i59 = 0; c2_i59 < 36; c2_i59++) {
    c2_dv3[c2_i59] = c2_dv2[c2_i59];
  }

  c2_c_eml_xtrsm(chartInstance, c2_dv3, c2_nu_dot);
  for (c2_i60 = 0; c2_i60 < 36; c2_i60++) {
    c2_dv4[c2_i60] = c2_q_a[c2_i60];
  }

  for (c2_i61 = 0; c2_i61 < 36; c2_i61++) {
    c2_dv5[c2_i61] = c2_dv4[c2_i61];
  }

  c2_d_eml_xtrsm(chartInstance, c2_dv5, c2_nu_dot);
  _SFD_EML_CALL(0U, chartInstance->c2_sfEvent, -25);
  _SFD_SYMBOL_SCOPE_POP();
  for (c2_i62 = 0; c2_i62 < 6; c2_i62++) {
    (*c2_b_nu_dot)[c2_i62] = c2_nu_dot[c2_i62];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c2_sfEvent);
}

static void initSimStructsc2_test(SFc2_testInstanceStruct *chartInstance)
{
}

static void registerMessagesc2_test(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_m2c(SFc2_testInstanceStruct *chartInstance, real_T c2_M[36],
                   real_T c2_nu[6], real_T c2_C[36])
{
  uint32_T c2_debug_family_var_map[14];
  real_T c2_nu1[3];
  real_T c2_nu2[3];
  real_T c2_Msym[36];
  real_T c2_M11[9];
  real_T c2_M12[9];
  real_T c2_M21[9];
  real_T c2_M22[9];
  real_T c2_dt_dnu1[3];
  real_T c2_dt_dnu2[3];
  real_T c2_nargin = 2.0;
  real_T c2_nargout = 1.0;
  int32_T c2_i63;
  int32_T c2_i64;
  int32_T c2_i65;
  int32_T c2_i66;
  int32_T c2_i67;
  int32_T c2_i68;
  real_T c2_b[36];
  int32_T c2_i69;
  int32_T c2_i70;
  int32_T c2_i71;
  int32_T c2_i72;
  int32_T c2_i73;
  int32_T c2_i74;
  int32_T c2_i75;
  int32_T c2_i76;
  int32_T c2_i77;
  int32_T c2_i78;
  int32_T c2_i79;
  int32_T c2_i80;
  int32_T c2_i81;
  int32_T c2_i82;
  int32_T c2_i83;
  int32_T c2_i84;
  int32_T c2_i85;
  int32_T c2_i86;
  int32_T c2_i87;
  real_T c2_a[9];
  int32_T c2_i88;
  real_T c2_b_b[3];
  int32_T c2_i89;
  real_T c2_y[3];
  int32_T c2_i90;
  int32_T c2_i91;
  int32_T c2_i92;
  int32_T c2_i93;
  int32_T c2_i94;
  real_T c2_b_y[3];
  int32_T c2_i95;
  int32_T c2_i96;
  int32_T c2_i97;
  int32_T c2_i98;
  int32_T c2_i99;
  int32_T c2_i100;
  int32_T c2_i101;
  int32_T c2_i102;
  int32_T c2_i103;
  int32_T c2_i104;
  int32_T c2_i105;
  int32_T c2_i106;
  int32_T c2_i107;
  int32_T c2_i108;
  int32_T c2_i109;
  real_T c2_r[3];
  uint32_T c2_b_debug_family_var_map[4];
  real_T c2_b_nargin = 1.0;
  real_T c2_b_nargout = 1.0;
  real_T c2_S[9];
  int32_T c2_i110;
  real_T c2_b_r[3];
  real_T c2_c_nargin = 1.0;
  real_T c2_c_nargout = 1.0;
  real_T c2_b_S[9];
  int32_T c2_i111;
  real_T c2_c_r[3];
  real_T c2_d_nargin = 1.0;
  real_T c2_d_nargout = 1.0;
  real_T c2_c_S[9];
  int32_T c2_i112;
  int32_T c2_i113;
  int32_T c2_i114;
  int32_T c2_i115;
  int32_T c2_i116;
  int32_T c2_i117;
  int32_T c2_i118;
  int32_T c2_i119;
  int32_T c2_i120;
  int32_T c2_i121;
  int32_T c2_i122;
  int32_T c2_i123;
  int32_T c2_i124;
  int32_T c2_i125;
  int32_T c2_i126;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 14U, 14U, c2_c_debug_family_names,
    c2_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_nu1, 0U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_nu2, 1U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_Msym, 2U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_M11, 3U, c2_d_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_M12, 4U, c2_d_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_M21, 5U, c2_d_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_M22, 6U, c2_d_sf_marshallOut,
    c2_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dt_dnu1, 7U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_dt_dnu2, 8U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargin, 9U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_nargout, 10U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_M, 11U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_nu, 12U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_C, 13U, c2_f_sf_marshallOut,
    c2_e_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 33);
  for (c2_i63 = 0; c2_i63 < 3; c2_i63++) {
    c2_nu1[c2_i63] = c2_nu[c2_i63];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 34);
  for (c2_i64 = 0; c2_i64 < 3; c2_i64++) {
    c2_nu2[c2_i64] = c2_nu[c2_i64 + 3];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 39);
  c2_i65 = 0;
  for (c2_i66 = 0; c2_i66 < 6; c2_i66++) {
    c2_i67 = 0;
    for (c2_i68 = 0; c2_i68 < 6; c2_i68++) {
      c2_b[c2_i68 + c2_i65] = c2_M[c2_i68 + c2_i65] + c2_M[c2_i67 + c2_i66];
      c2_i67 += 6;
    }

    c2_i65 += 6;
  }

  for (c2_i69 = 0; c2_i69 < 36; c2_i69++) {
    c2_Msym[c2_i69] = 0.5 * c2_b[c2_i69];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 40);
  for (c2_i70 = 0; c2_i70 < 36; c2_i70++) {
    c2_M[c2_i70] = c2_Msym[c2_i70];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 42);
  c2_i71 = 0;
  c2_i72 = 0;
  for (c2_i73 = 0; c2_i73 < 3; c2_i73++) {
    for (c2_i74 = 0; c2_i74 < 3; c2_i74++) {
      c2_M11[c2_i74 + c2_i71] = c2_M[c2_i74 + c2_i72];
    }

    c2_i71 += 3;
    c2_i72 += 6;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 43);
  c2_i75 = 0;
  c2_i76 = 0;
  for (c2_i77 = 0; c2_i77 < 3; c2_i77++) {
    for (c2_i78 = 0; c2_i78 < 3; c2_i78++) {
      c2_M12[c2_i78 + c2_i75] = c2_M[(c2_i78 + c2_i76) + 18];
    }

    c2_i75 += 3;
    c2_i76 += 6;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 44);
  c2_i79 = 0;
  for (c2_i80 = 0; c2_i80 < 3; c2_i80++) {
    c2_i81 = 0;
    for (c2_i82 = 0; c2_i82 < 3; c2_i82++) {
      c2_M21[c2_i82 + c2_i79] = c2_M12[c2_i81 + c2_i80];
      c2_i81 += 3;
    }

    c2_i79 += 3;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 45);
  c2_i83 = 0;
  c2_i84 = 0;
  for (c2_i85 = 0; c2_i85 < 3; c2_i85++) {
    for (c2_i86 = 0; c2_i86 < 3; c2_i86++) {
      c2_M22[c2_i86 + c2_i83] = c2_M[(c2_i86 + c2_i84) + 21];
    }

    c2_i83 += 3;
    c2_i84 += 6;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 47);
  for (c2_i87 = 0; c2_i87 < 9; c2_i87++) {
    c2_a[c2_i87] = c2_M11[c2_i87];
  }

  for (c2_i88 = 0; c2_i88 < 3; c2_i88++) {
    c2_b_b[c2_i88] = c2_nu1[c2_i88];
  }

  c2_eml_scalar_eg(chartInstance);
  c2_eml_scalar_eg(chartInstance);
  for (c2_i89 = 0; c2_i89 < 3; c2_i89++) {
    c2_y[c2_i89] = 0.0;
    c2_i90 = 0;
    for (c2_i91 = 0; c2_i91 < 3; c2_i91++) {
      c2_y[c2_i89] += c2_a[c2_i90 + c2_i89] * c2_b_b[c2_i91];
      c2_i90 += 3;
    }
  }

  for (c2_i92 = 0; c2_i92 < 9; c2_i92++) {
    c2_a[c2_i92] = c2_M12[c2_i92];
  }

  for (c2_i93 = 0; c2_i93 < 3; c2_i93++) {
    c2_b_b[c2_i93] = c2_nu2[c2_i93];
  }

  c2_eml_scalar_eg(chartInstance);
  c2_eml_scalar_eg(chartInstance);
  for (c2_i94 = 0; c2_i94 < 3; c2_i94++) {
    c2_b_y[c2_i94] = 0.0;
    c2_i95 = 0;
    for (c2_i96 = 0; c2_i96 < 3; c2_i96++) {
      c2_b_y[c2_i94] += c2_a[c2_i95 + c2_i94] * c2_b_b[c2_i96];
      c2_i95 += 3;
    }
  }

  for (c2_i97 = 0; c2_i97 < 3; c2_i97++) {
    c2_dt_dnu1[c2_i97] = c2_y[c2_i97] + c2_b_y[c2_i97];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 48);
  for (c2_i98 = 0; c2_i98 < 9; c2_i98++) {
    c2_a[c2_i98] = c2_M21[c2_i98];
  }

  for (c2_i99 = 0; c2_i99 < 3; c2_i99++) {
    c2_b_b[c2_i99] = c2_nu1[c2_i99];
  }

  c2_eml_scalar_eg(chartInstance);
  c2_eml_scalar_eg(chartInstance);
  for (c2_i100 = 0; c2_i100 < 3; c2_i100++) {
    c2_y[c2_i100] = 0.0;
    c2_i101 = 0;
    for (c2_i102 = 0; c2_i102 < 3; c2_i102++) {
      c2_y[c2_i100] += c2_a[c2_i101 + c2_i100] * c2_b_b[c2_i102];
      c2_i101 += 3;
    }
  }

  for (c2_i103 = 0; c2_i103 < 9; c2_i103++) {
    c2_a[c2_i103] = c2_M22[c2_i103];
  }

  for (c2_i104 = 0; c2_i104 < 3; c2_i104++) {
    c2_b_b[c2_i104] = c2_nu2[c2_i104];
  }

  c2_eml_scalar_eg(chartInstance);
  c2_eml_scalar_eg(chartInstance);
  for (c2_i105 = 0; c2_i105 < 3; c2_i105++) {
    c2_b_y[c2_i105] = 0.0;
    c2_i106 = 0;
    for (c2_i107 = 0; c2_i107 < 3; c2_i107++) {
      c2_b_y[c2_i105] += c2_a[c2_i106 + c2_i105] * c2_b_b[c2_i107];
      c2_i106 += 3;
    }
  }

  for (c2_i108 = 0; c2_i108 < 3; c2_i108++) {
    c2_dt_dnu2[c2_i108] = c2_y[c2_i108] + c2_b_y[c2_i108];
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, 50);
  for (c2_i109 = 0; c2_i109 < 3; c2_i109++) {
    c2_r[c2_i109] = c2_dt_dnu1[c2_i109];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c2_b_debug_family_names,
    c2_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_nargin, 0U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_b_nargout, 1U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_r, 2U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_S, 3U, c2_d_sf_marshallOut,
    c2_f_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, 31);
  c2_S[0] = 0.0;
  c2_S[3] = -c2_r[2];
  c2_S[6] = c2_r[1];
  c2_S[1] = c2_r[2];
  c2_S[4] = 0.0;
  c2_S[7] = -c2_r[0];
  c2_S[2] = -c2_r[1];
  c2_S[5] = c2_r[0];
  c2_S[8] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  for (c2_i110 = 0; c2_i110 < 3; c2_i110++) {
    c2_b_r[c2_i110] = c2_dt_dnu1[c2_i110];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c2_b_debug_family_names,
    c2_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c_nargin, 0U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_c_nargout, 1U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_b_r, 2U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_b_S, 3U, c2_d_sf_marshallOut,
    c2_f_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, 31);
  c2_b_S[0] = 0.0;
  c2_b_S[3] = -c2_b_r[2];
  c2_b_S[6] = c2_b_r[1];
  c2_b_S[1] = c2_b_r[2];
  c2_b_S[4] = 0.0;
  c2_b_S[7] = -c2_b_r[0];
  c2_b_S[2] = -c2_b_r[1];
  c2_b_S[5] = c2_b_r[0];
  c2_b_S[8] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  for (c2_i111 = 0; c2_i111 < 3; c2_i111++) {
    c2_c_r[c2_i111] = c2_dt_dnu2[c2_i111];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 4U, 4U, c2_b_debug_family_names,
    c2_b_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_d_nargin, 0U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c2_d_nargout, 1U, c2_b_sf_marshallOut,
    c2_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_c_r, 2U, c2_g_sf_marshallOut,
    c2_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c2_c_S, 3U, c2_d_sf_marshallOut,
    c2_f_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, 31);
  c2_c_S[0] = 0.0;
  c2_c_S[3] = -c2_c_r[2];
  c2_c_S[6] = c2_c_r[1];
  c2_c_S[1] = c2_c_r[2];
  c2_c_S[4] = 0.0;
  c2_c_S[7] = -c2_c_r[0];
  c2_c_S[2] = -c2_c_r[1];
  c2_c_S[5] = c2_c_r[0];
  c2_c_S[8] = 0.0;
  _SFD_SCRIPT_CALL(1U, chartInstance->c2_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  c2_i112 = 0;
  for (c2_i113 = 0; c2_i113 < 3; c2_i113++) {
    for (c2_i114 = 0; c2_i114 < 3; c2_i114++) {
      c2_C[c2_i114 + c2_i112] = 0.0;
    }

    c2_i112 += 6;
  }

  c2_i115 = 0;
  c2_i116 = 0;
  for (c2_i117 = 0; c2_i117 < 3; c2_i117++) {
    for (c2_i118 = 0; c2_i118 < 3; c2_i118++) {
      c2_C[(c2_i118 + c2_i115) + 18] = -c2_S[c2_i118 + c2_i116];
    }

    c2_i115 += 6;
    c2_i116 += 3;
  }

  c2_i119 = 0;
  c2_i120 = 0;
  for (c2_i121 = 0; c2_i121 < 3; c2_i121++) {
    for (c2_i122 = 0; c2_i122 < 3; c2_i122++) {
      c2_C[(c2_i122 + c2_i119) + 3] = -c2_b_S[c2_i122 + c2_i120];
    }

    c2_i119 += 6;
    c2_i120 += 3;
  }

  c2_i123 = 0;
  c2_i124 = 0;
  for (c2_i125 = 0; c2_i125 < 3; c2_i125++) {
    for (c2_i126 = 0; c2_i126 < 3; c2_i126++) {
      c2_C[(c2_i126 + c2_i123) + 21] = -c2_c_S[c2_i126 + c2_i124];
    }

    c2_i123 += 6;
    c2_i124 += 3;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c2_sfEvent, -50);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber)
{
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, 0U, sf_debug_get_script_id(
    "C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/m2c.m"));
  _SFD_SCRIPT_TRANSLATION(c2_chartNumber, 1U, sf_debug_get_script_id(
    "C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/Smtrx.m"));
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i127;
  real_T c2_b_inData[6];
  int32_T c2_i128;
  real_T c2_u[6];
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i127 = 0; c2_i127 < 6; c2_i127++) {
    c2_b_inData[c2_i127] = (*(real_T (*)[6])c2_inData)[c2_i127];
  }

  for (c2_i128 = 0; c2_i128 < 6; c2_i128++) {
    c2_u[c2_i128] = c2_b_inData[c2_i128];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static void c2_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_nu_dot, const char_T *c2_identifier, real_T c2_y[6])
{
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_nu_dot), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_nu_dot);
}

static void c2_b_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[6])
{
  real_T c2_dv6[6];
  int32_T c2_i129;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv6, 1, 0, 0U, 1, 0U, 1, 6);
  for (c2_i129 = 0; c2_i129 < 6; c2_i129++) {
    c2_y[c2_i129] = c2_dv6[c2_i129];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_nu_dot;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[6];
  int32_T c2_i130;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_nu_dot = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_nu_dot), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_nu_dot);
  for (c2_i130 = 0; c2_i130 < 6; c2_i130++) {
    (*(real_T (*)[6])c2_outData)[c2_i130] = c2_y[c2_i130];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static const mxArray *c2_c_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i131;
  real_T c2_b_inData[2];
  int32_T c2_i132;
  real_T c2_u[2];
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i131 = 0; c2_i131 < 2; c2_i131++) {
    c2_b_inData[c2_i131] = (*(real_T (*)[2])c2_inData)[c2_i131];
  }

  for (c2_i132 = 0; c2_i132 < 2; c2_i132++) {
    c2_u[c2_i132] = c2_b_inData[c2_i132];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static const mxArray *c2_d_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i133;
  int32_T c2_i134;
  int32_T c2_i135;
  real_T c2_b_inData[9];
  int32_T c2_i136;
  int32_T c2_i137;
  int32_T c2_i138;
  real_T c2_u[9];
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_i133 = 0;
  for (c2_i134 = 0; c2_i134 < 3; c2_i134++) {
    for (c2_i135 = 0; c2_i135 < 3; c2_i135++) {
      c2_b_inData[c2_i135 + c2_i133] = (*(real_T (*)[9])c2_inData)[c2_i135 +
        c2_i133];
    }

    c2_i133 += 3;
  }

  c2_i136 = 0;
  for (c2_i137 = 0; c2_i137 < 3; c2_i137++) {
    for (c2_i138 = 0; c2_i138 < 3; c2_i138++) {
      c2_u[c2_i138 + c2_i136] = c2_b_inData[c2_i138 + c2_i136];
    }

    c2_i136 += 3;
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static real_T c2_c_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d0;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_nargout;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_nargout = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_nargout), &c2_thisId);
  sf_mex_destroy(&c2_nargout);
  *(real_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static void c2_d_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[2])
{
  real_T c2_dv7[2];
  int32_T c2_i139;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv7, 1, 0, 0U, 1, 0U, 1, 2);
  for (c2_i139 = 0; c2_i139 < 2; c2_i139++) {
    c2_y[c2_i139] = c2_dv7[c2_i139];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_wind_b;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[2];
  int32_T c2_i140;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_wind_b = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_wind_b), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_wind_b);
  for (c2_i140 = 0; c2_i140 < 2; c2_i140++) {
    (*(real_T (*)[2])c2_outData)[c2_i140] = c2_y[c2_i140];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_e_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i141;
  int32_T c2_i142;
  int32_T c2_i143;
  real_T c2_b_inData[4];
  int32_T c2_i144;
  int32_T c2_i145;
  int32_T c2_i146;
  real_T c2_u[4];
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_i141 = 0;
  for (c2_i142 = 0; c2_i142 < 2; c2_i142++) {
    for (c2_i143 = 0; c2_i143 < 2; c2_i143++) {
      c2_b_inData[c2_i143 + c2_i141] = (*(real_T (*)[4])c2_inData)[c2_i143 +
        c2_i141];
    }

    c2_i141 += 2;
  }

  c2_i144 = 0;
  for (c2_i145 = 0; c2_i145 < 2; c2_i145++) {
    for (c2_i146 = 0; c2_i146 < 2; c2_i146++) {
      c2_u[c2_i146 + c2_i144] = c2_b_inData[c2_i146 + c2_i144];
    }

    c2_i144 += 2;
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 2, 2), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static void c2_e_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[4])
{
  real_T c2_dv8[4];
  int32_T c2_i147;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv8, 1, 0, 0U, 1, 0U, 2, 2, 2);
  for (c2_i147 = 0; c2_i147 < 4; c2_i147++) {
    c2_y[c2_i147] = c2_dv8[c2_i147];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_RotWind;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[4];
  int32_T c2_i148;
  int32_T c2_i149;
  int32_T c2_i150;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_RotWind = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_RotWind), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_RotWind);
  c2_i148 = 0;
  for (c2_i149 = 0; c2_i149 < 2; c2_i149++) {
    for (c2_i150 = 0; c2_i150 < 2; c2_i150++) {
      (*(real_T (*)[4])c2_outData)[c2_i150 + c2_i148] = c2_y[c2_i150 + c2_i148];
    }

    c2_i148 += 2;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_f_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i151;
  int32_T c2_i152;
  int32_T c2_i153;
  real_T c2_b_inData[36];
  int32_T c2_i154;
  int32_T c2_i155;
  int32_T c2_i156;
  real_T c2_u[36];
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_i151 = 0;
  for (c2_i152 = 0; c2_i152 < 6; c2_i152++) {
    for (c2_i153 = 0; c2_i153 < 6; c2_i153++) {
      c2_b_inData[c2_i153 + c2_i151] = (*(real_T (*)[36])c2_inData)[c2_i153 +
        c2_i151];
    }

    c2_i151 += 6;
  }

  c2_i154 = 0;
  for (c2_i155 = 0; c2_i155 < 6; c2_i155++) {
    for (c2_i156 = 0; c2_i156 < 6; c2_i156++) {
      c2_u[c2_i156 + c2_i154] = c2_b_inData[c2_i156 + c2_i154];
    }

    c2_i154 += 6;
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 2, 6, 6), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static void c2_f_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[36])
{
  real_T c2_dv9[36];
  int32_T c2_i157;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv9, 1, 0, 0U, 1, 0U, 2, 6, 6);
  for (c2_i157 = 0; c2_i157 < 36; c2_i157++) {
    c2_y[c2_i157] = c2_dv9[c2_i157];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_C_rb;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[36];
  int32_T c2_i158;
  int32_T c2_i159;
  int32_T c2_i160;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_C_rb = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_C_rb), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_C_rb);
  c2_i158 = 0;
  for (c2_i159 = 0; c2_i159 < 6; c2_i159++) {
    for (c2_i160 = 0; c2_i160 < 6; c2_i160++) {
      (*(real_T (*)[36])c2_outData)[c2_i160 + c2_i158] = c2_y[c2_i160 + c2_i158];
    }

    c2_i158 += 6;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static void c2_g_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[9])
{
  real_T c2_dv10[9];
  int32_T c2_i161;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv10, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c2_i161 = 0; c2_i161 < 9; c2_i161++) {
    c2_y[c2_i161] = c2_dv10[c2_i161];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_S;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[9];
  int32_T c2_i162;
  int32_T c2_i163;
  int32_T c2_i164;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_S = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_S), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_S);
  c2_i162 = 0;
  for (c2_i163 = 0; c2_i163 < 3; c2_i163++) {
    for (c2_i164 = 0; c2_i164 < 3; c2_i164++) {
      (*(real_T (*)[9])c2_outData)[c2_i164 + c2_i162] = c2_y[c2_i164 + c2_i162];
    }

    c2_i162 += 3;
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

static const mxArray *c2_g_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_i165;
  real_T c2_b_inData[3];
  int32_T c2_i166;
  real_T c2_u[3];
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  for (c2_i165 = 0; c2_i165 < 3; c2_i165++) {
    c2_b_inData[c2_i165] = (*(real_T (*)[3])c2_inData)[c2_i165];
  }

  for (c2_i166 = 0; c2_i166 < 3; c2_i166++) {
    c2_u[c2_i166] = c2_b_inData[c2_i166];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static void c2_h_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance, const
  mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId, real_T c2_y[3])
{
  real_T c2_dv11[3];
  int32_T c2_i167;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), c2_dv11, 1, 0, 0U, 1, 0U, 1, 3);
  for (c2_i167 = 0; c2_i167 < 3; c2_i167++) {
    c2_y[c2_i167] = c2_dv11[c2_i167];
  }

  sf_mex_destroy(&c2_u);
}

static void c2_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_r;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y[3];
  int32_T c2_i168;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_r = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_r), &c2_thisId, c2_y);
  sf_mex_destroy(&c2_r);
  for (c2_i168 = 0; c2_i168 < 3; c2_i168++) {
    (*(real_T (*)[3])c2_outData)[c2_i168] = c2_y[c2_i168];
  }

  sf_mex_destroy(&c2_mxArrayInData);
}

const mxArray *sf_c2_test_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo;
  c2_ResolvedFunctionInfo c2_info[159];
  const mxArray *c2_m0 = NULL;
  int32_T c2_i169;
  c2_ResolvedFunctionInfo *c2_r0;
  c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  c2_info_helper(c2_info);
  c2_b_info_helper(c2_info);
  c2_c_info_helper(c2_info);
  sf_mex_assign(&c2_m0, sf_mex_createstruct("nameCaptureInfo", 1, 159), FALSE);
  for (c2_i169 = 0; c2_i169 < 159; c2_i169++) {
    c2_r0 = &c2_info[c2_i169];
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->context)), "context", "nameCaptureInfo",
                    c2_i169);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c2_r0->name)), "name", "nameCaptureInfo", c2_i169);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c2_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c2_i169);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->resolved)), "resolved", "nameCaptureInfo",
                    c2_i169);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c2_i169);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c2_i169);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c2_i169);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c2_i169);
  }

  sf_mex_assign(&c2_nameCaptureInfo, c2_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c2_nameCaptureInfo);
  return c2_nameCaptureInfo;
}

static void c2_info_helper(c2_ResolvedFunctionInfo c2_info[159])
{
  c2_info[0].context = "";
  c2_info[0].name = "eye";
  c2_info[0].dominantType = "double";
  c2_info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m";
  c2_info[0].fileTimeLo = 1286818688U;
  c2_info[0].fileTimeHi = 0U;
  c2_info[0].mFileTimeLo = 0U;
  c2_info[0].mFileTimeHi = 0U;
  c2_info[1].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[1].name = "eml_assert_valid_size_arg";
  c2_info[1].dominantType = "double";
  c2_info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c2_info[1].fileTimeLo = 1286818694U;
  c2_info[1].fileTimeHi = 0U;
  c2_info[1].mFileTimeLo = 0U;
  c2_info[1].mFileTimeHi = 0U;
  c2_info[2].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral";
  c2_info[2].name = "isinf";
  c2_info[2].dominantType = "double";
  c2_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m";
  c2_info[2].fileTimeLo = 1286818760U;
  c2_info[2].fileTimeHi = 0U;
  c2_info[2].mFileTimeLo = 0U;
  c2_info[2].mFileTimeHi = 0U;
  c2_info[3].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!numel_for_size";
  c2_info[3].name = "mtimes";
  c2_info[3].dominantType = "double";
  c2_info[3].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[3].fileTimeLo = 1289519692U;
  c2_info[3].fileTimeHi = 0U;
  c2_info[3].mFileTimeLo = 0U;
  c2_info[3].mFileTimeHi = 0U;
  c2_info[4].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c2_info[4].name = "eml_index_class";
  c2_info[4].dominantType = "";
  c2_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[4].fileTimeLo = 1323170578U;
  c2_info[4].fileTimeHi = 0U;
  c2_info[4].mFileTimeLo = 0U;
  c2_info[4].mFileTimeHi = 0U;
  c2_info[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  c2_info[5].name = "intmax";
  c2_info[5].dominantType = "char";
  c2_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[5].fileTimeLo = 1311255316U;
  c2_info[5].fileTimeHi = 0U;
  c2_info[5].mFileTimeLo = 0U;
  c2_info[5].mFileTimeHi = 0U;
  c2_info[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[6].name = "eml_is_float_class";
  c2_info[6].dominantType = "char";
  c2_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  c2_info[6].fileTimeLo = 1286818782U;
  c2_info[6].fileTimeHi = 0U;
  c2_info[6].mFileTimeLo = 0U;
  c2_info[6].mFileTimeHi = 0U;
  c2_info[7].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[7].name = "min";
  c2_info[7].dominantType = "double";
  c2_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c2_info[7].fileTimeLo = 1311255318U;
  c2_info[7].fileTimeHi = 0U;
  c2_info[7].mFileTimeLo = 0U;
  c2_info[7].mFileTimeHi = 0U;
  c2_info[8].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c2_info[8].name = "eml_min_or_max";
  c2_info[8].dominantType = "char";
  c2_info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  c2_info[8].fileTimeLo = 1334071490U;
  c2_info[8].fileTimeHi = 0U;
  c2_info[8].mFileTimeLo = 0U;
  c2_info[8].mFileTimeHi = 0U;
  c2_info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c2_info[9].name = "eml_scalar_eg";
  c2_info[9].dominantType = "double";
  c2_info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[9].fileTimeLo = 1286818796U;
  c2_info[9].fileTimeHi = 0U;
  c2_info[9].mFileTimeLo = 0U;
  c2_info[9].mFileTimeHi = 0U;
  c2_info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c2_info[10].name = "eml_scalexp_alloc";
  c2_info[10].dominantType = "double";
  c2_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c2_info[10].fileTimeLo = 1352424860U;
  c2_info[10].fileTimeHi = 0U;
  c2_info[10].mFileTimeLo = 0U;
  c2_info[10].mFileTimeHi = 0U;
  c2_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c2_info[11].name = "eml_index_class";
  c2_info[11].dominantType = "";
  c2_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[11].fileTimeLo = 1323170578U;
  c2_info[11].fileTimeHi = 0U;
  c2_info[11].mFileTimeLo = 0U;
  c2_info[11].mFileTimeHi = 0U;
  c2_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum";
  c2_info[12].name = "eml_scalar_eg";
  c2_info[12].dominantType = "double";
  c2_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[12].fileTimeLo = 1286818796U;
  c2_info[12].fileTimeHi = 0U;
  c2_info[12].mFileTimeLo = 0U;
  c2_info[12].mFileTimeHi = 0U;
  c2_info[13].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[13].name = "eml_index_class";
  c2_info[13].dominantType = "";
  c2_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[13].fileTimeLo = 1323170578U;
  c2_info[13].fileTimeHi = 0U;
  c2_info[13].mFileTimeLo = 0U;
  c2_info[13].mFileTimeHi = 0U;
  c2_info[14].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eye.m!eye_internal";
  c2_info[14].name = "eml_int_forloop_overflow_check";
  c2_info[14].dominantType = "";
  c2_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[14].fileTimeLo = 1346510340U;
  c2_info[14].fileTimeHi = 0U;
  c2_info[14].mFileTimeLo = 0U;
  c2_info[14].mFileTimeHi = 0U;
  c2_info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  c2_info[15].name = "intmax";
  c2_info[15].dominantType = "char";
  c2_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[15].fileTimeLo = 1311255316U;
  c2_info[15].fileTimeHi = 0U;
  c2_info[15].mFileTimeLo = 0U;
  c2_info[15].mFileTimeHi = 0U;
  c2_info[16].context = "";
  c2_info[16].name = "mtimes";
  c2_info[16].dominantType = "double";
  c2_info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[16].fileTimeLo = 1289519692U;
  c2_info[16].fileTimeHi = 0U;
  c2_info[16].mFileTimeLo = 0U;
  c2_info[16].mFileTimeHi = 0U;
  c2_info[17].context = "";
  c2_info[17].name = "m2c";
  c2_info[17].dominantType = "double";
  c2_info[17].resolved = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/m2c.m";
  c2_info[17].fileTimeLo = 1206483598U;
  c2_info[17].fileTimeHi = 0U;
  c2_info[17].mFileTimeLo = 0U;
  c2_info[17].mFileTimeHi = 0U;
  c2_info[18].context = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/m2c.m";
  c2_info[18].name = "mtimes";
  c2_info[18].dominantType = "double";
  c2_info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[18].fileTimeLo = 1289519692U;
  c2_info[18].fileTimeHi = 0U;
  c2_info[18].mFileTimeLo = 0U;
  c2_info[18].mFileTimeHi = 0U;
  c2_info[19].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[19].name = "eml_index_class";
  c2_info[19].dominantType = "";
  c2_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[19].fileTimeLo = 1323170578U;
  c2_info[19].fileTimeHi = 0U;
  c2_info[19].mFileTimeLo = 0U;
  c2_info[19].mFileTimeHi = 0U;
  c2_info[20].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[20].name = "eml_scalar_eg";
  c2_info[20].dominantType = "double";
  c2_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[20].fileTimeLo = 1286818796U;
  c2_info[20].fileTimeHi = 0U;
  c2_info[20].mFileTimeLo = 0U;
  c2_info[20].mFileTimeHi = 0U;
  c2_info[21].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[21].name = "eml_xgemm";
  c2_info[21].dominantType = "char";
  c2_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c2_info[21].fileTimeLo = 1299076772U;
  c2_info[21].fileTimeHi = 0U;
  c2_info[21].mFileTimeLo = 0U;
  c2_info[21].mFileTimeHi = 0U;
  c2_info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c2_info[22].name = "eml_blas_inline";
  c2_info[22].dominantType = "";
  c2_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[22].fileTimeLo = 1299076768U;
  c2_info[22].fileTimeHi = 0U;
  c2_info[22].mFileTimeLo = 0U;
  c2_info[22].mFileTimeHi = 0U;
  c2_info[23].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold";
  c2_info[23].name = "mtimes";
  c2_info[23].dominantType = "double";
  c2_info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[23].fileTimeLo = 1289519692U;
  c2_info[23].fileTimeHi = 0U;
  c2_info[23].mFileTimeLo = 0U;
  c2_info[23].mFileTimeHi = 0U;
  c2_info[24].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c2_info[24].name = "eml_index_class";
  c2_info[24].dominantType = "";
  c2_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[24].fileTimeLo = 1323170578U;
  c2_info[24].fileTimeHi = 0U;
  c2_info[24].mFileTimeLo = 0U;
  c2_info[24].mFileTimeHi = 0U;
  c2_info[25].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c2_info[25].name = "eml_scalar_eg";
  c2_info[25].dominantType = "double";
  c2_info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[25].fileTimeLo = 1286818796U;
  c2_info[25].fileTimeHi = 0U;
  c2_info[25].mFileTimeLo = 0U;
  c2_info[25].mFileTimeHi = 0U;
  c2_info[26].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c2_info[26].name = "eml_refblas_xgemm";
  c2_info[26].dominantType = "char";
  c2_info[26].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c2_info[26].fileTimeLo = 1299076774U;
  c2_info[26].fileTimeHi = 0U;
  c2_info[26].mFileTimeLo = 0U;
  c2_info[26].mFileTimeHi = 0U;
  c2_info[27].context = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/m2c.m";
  c2_info[27].name = "Smtrx";
  c2_info[27].dominantType = "double";
  c2_info[27].resolved = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/Smtrx.m";
  c2_info[27].fileTimeLo = 1206483740U;
  c2_info[27].fileTimeHi = 0U;
  c2_info[27].mFileTimeLo = 0U;
  c2_info[27].mFileTimeHi = 0U;
  c2_info[28].context = "";
  c2_info[28].name = "cos";
  c2_info[28].dominantType = "double";
  c2_info[28].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c2_info[28].fileTimeLo = 1343830372U;
  c2_info[28].fileTimeHi = 0U;
  c2_info[28].mFileTimeLo = 0U;
  c2_info[28].mFileTimeHi = 0U;
  c2_info[29].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c2_info[29].name = "eml_scalar_cos";
  c2_info[29].dominantType = "double";
  c2_info[29].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  c2_info[29].fileTimeLo = 1286818722U;
  c2_info[29].fileTimeHi = 0U;
  c2_info[29].mFileTimeLo = 0U;
  c2_info[29].mFileTimeHi = 0U;
  c2_info[30].context = "";
  c2_info[30].name = "sin";
  c2_info[30].dominantType = "double";
  c2_info[30].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c2_info[30].fileTimeLo = 1343830386U;
  c2_info[30].fileTimeHi = 0U;
  c2_info[30].mFileTimeLo = 0U;
  c2_info[30].mFileTimeHi = 0U;
  c2_info[31].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c2_info[31].name = "eml_scalar_sin";
  c2_info[31].dominantType = "double";
  c2_info[31].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  c2_info[31].fileTimeLo = 1286818736U;
  c2_info[31].fileTimeHi = 0U;
  c2_info[31].mFileTimeLo = 0U;
  c2_info[31].mFileTimeHi = 0U;
  c2_info[32].context = "";
  c2_info[32].name = "sum";
  c2_info[32].dominantType = "double";
  c2_info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c2_info[32].fileTimeLo = 1314736612U;
  c2_info[32].fileTimeHi = 0U;
  c2_info[32].mFileTimeLo = 0U;
  c2_info[32].mFileTimeHi = 0U;
  c2_info[33].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c2_info[33].name = "isequal";
  c2_info[33].dominantType = "double";
  c2_info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c2_info[33].fileTimeLo = 1286818758U;
  c2_info[33].fileTimeHi = 0U;
  c2_info[33].mFileTimeLo = 0U;
  c2_info[33].mFileTimeHi = 0U;
  c2_info[34].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c2_info[34].name = "eml_isequal_core";
  c2_info[34].dominantType = "double";
  c2_info[34].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  c2_info[34].fileTimeLo = 1286818786U;
  c2_info[34].fileTimeHi = 0U;
  c2_info[34].mFileTimeLo = 0U;
  c2_info[34].mFileTimeHi = 0U;
  c2_info[35].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c2_info[35].name = "eml_const_nonsingleton_dim";
  c2_info[35].dominantType = "double";
  c2_info[35].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m";
  c2_info[35].fileTimeLo = 1286818696U;
  c2_info[35].fileTimeHi = 0U;
  c2_info[35].mFileTimeLo = 0U;
  c2_info[35].mFileTimeHi = 0U;
  c2_info[36].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c2_info[36].name = "eml_scalar_eg";
  c2_info[36].dominantType = "double";
  c2_info[36].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[36].fileTimeLo = 1286818796U;
  c2_info[36].fileTimeHi = 0U;
  c2_info[36].mFileTimeLo = 0U;
  c2_info[36].mFileTimeHi = 0U;
  c2_info[37].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c2_info[37].name = "eml_index_class";
  c2_info[37].dominantType = "";
  c2_info[37].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[37].fileTimeLo = 1323170578U;
  c2_info[37].fileTimeHi = 0U;
  c2_info[37].mFileTimeLo = 0U;
  c2_info[37].mFileTimeHi = 0U;
  c2_info[38].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c2_info[38].name = "eml_int_forloop_overflow_check";
  c2_info[38].dominantType = "";
  c2_info[38].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[38].fileTimeLo = 1346510340U;
  c2_info[38].fileTimeHi = 0U;
  c2_info[38].mFileTimeLo = 0U;
  c2_info[38].mFileTimeHi = 0U;
  c2_info[39].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[39].name = "eml_xdotu";
  c2_info[39].dominantType = "double";
  c2_info[39].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m";
  c2_info[39].fileTimeLo = 1299076772U;
  c2_info[39].fileTimeHi = 0U;
  c2_info[39].mFileTimeLo = 0U;
  c2_info[39].mFileTimeHi = 0U;
  c2_info[40].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m";
  c2_info[40].name = "eml_blas_inline";
  c2_info[40].dominantType = "";
  c2_info[40].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[40].fileTimeLo = 1299076768U;
  c2_info[40].fileTimeHi = 0U;
  c2_info[40].mFileTimeLo = 0U;
  c2_info[40].mFileTimeHi = 0U;
  c2_info[41].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m";
  c2_info[41].name = "eml_xdot";
  c2_info[41].dominantType = "double";
  c2_info[41].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m";
  c2_info[41].fileTimeLo = 1299076772U;
  c2_info[41].fileTimeHi = 0U;
  c2_info[41].mFileTimeLo = 0U;
  c2_info[41].mFileTimeHi = 0U;
  c2_info[42].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m";
  c2_info[42].name = "eml_blas_inline";
  c2_info[42].dominantType = "";
  c2_info[42].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[42].fileTimeLo = 1299076768U;
  c2_info[42].fileTimeHi = 0U;
  c2_info[42].mFileTimeLo = 0U;
  c2_info[42].mFileTimeHi = 0U;
  c2_info[43].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m";
  c2_info[43].name = "eml_index_class";
  c2_info[43].dominantType = "";
  c2_info[43].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[43].fileTimeLo = 1323170578U;
  c2_info[43].fileTimeHi = 0U;
  c2_info[43].mFileTimeLo = 0U;
  c2_info[43].mFileTimeHi = 0U;
  c2_info[44].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m";
  c2_info[44].name = "eml_refblas_xdot";
  c2_info[44].dominantType = "double";
  c2_info[44].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c2_info[44].fileTimeLo = 1299076772U;
  c2_info[44].fileTimeHi = 0U;
  c2_info[44].mFileTimeLo = 0U;
  c2_info[44].mFileTimeHi = 0U;
  c2_info[45].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m";
  c2_info[45].name = "eml_refblas_xdotx";
  c2_info[45].dominantType = "char";
  c2_info[45].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[45].fileTimeLo = 1299076774U;
  c2_info[45].fileTimeHi = 0U;
  c2_info[45].mFileTimeLo = 0U;
  c2_info[45].mFileTimeHi = 0U;
  c2_info[46].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[46].name = "eml_scalar_eg";
  c2_info[46].dominantType = "double";
  c2_info[46].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[46].fileTimeLo = 1286818796U;
  c2_info[46].fileTimeHi = 0U;
  c2_info[46].mFileTimeLo = 0U;
  c2_info[46].mFileTimeHi = 0U;
  c2_info[47].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[47].name = "eml_index_class";
  c2_info[47].dominantType = "";
  c2_info[47].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[47].fileTimeLo = 1323170578U;
  c2_info[47].fileTimeHi = 0U;
  c2_info[47].mFileTimeLo = 0U;
  c2_info[47].mFileTimeHi = 0U;
  c2_info[48].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[48].name = "eml_index_minus";
  c2_info[48].dominantType = "double";
  c2_info[48].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[48].fileTimeLo = 1286818778U;
  c2_info[48].fileTimeHi = 0U;
  c2_info[48].mFileTimeLo = 0U;
  c2_info[48].mFileTimeHi = 0U;
  c2_info[49].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[49].name = "eml_index_class";
  c2_info[49].dominantType = "";
  c2_info[49].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[49].fileTimeLo = 1323170578U;
  c2_info[49].fileTimeHi = 0U;
  c2_info[49].mFileTimeLo = 0U;
  c2_info[49].mFileTimeHi = 0U;
  c2_info[50].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[50].name = "eml_index_times";
  c2_info[50].dominantType = "coder.internal.indexInt";
  c2_info[50].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c2_info[50].fileTimeLo = 1286818780U;
  c2_info[50].fileTimeHi = 0U;
  c2_info[50].mFileTimeLo = 0U;
  c2_info[50].mFileTimeHi = 0U;
  c2_info[51].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c2_info[51].name = "eml_index_class";
  c2_info[51].dominantType = "";
  c2_info[51].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[51].fileTimeLo = 1323170578U;
  c2_info[51].fileTimeHi = 0U;
  c2_info[51].mFileTimeLo = 0U;
  c2_info[51].mFileTimeHi = 0U;
  c2_info[52].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[52].name = "eml_index_plus";
  c2_info[52].dominantType = "coder.internal.indexInt";
  c2_info[52].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[52].fileTimeLo = 1286818778U;
  c2_info[52].fileTimeHi = 0U;
  c2_info[52].mFileTimeLo = 0U;
  c2_info[52].mFileTimeHi = 0U;
  c2_info[53].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[53].name = "eml_index_class";
  c2_info[53].dominantType = "";
  c2_info[53].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[53].fileTimeLo = 1323170578U;
  c2_info[53].fileTimeHi = 0U;
  c2_info[53].mFileTimeLo = 0U;
  c2_info[53].mFileTimeHi = 0U;
  c2_info[54].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m";
  c2_info[54].name = "eml_int_forloop_overflow_check";
  c2_info[54].dominantType = "";
  c2_info[54].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[54].fileTimeLo = 1346510340U;
  c2_info[54].fileTimeHi = 0U;
  c2_info[54].mFileTimeLo = 0U;
  c2_info[54].mFileTimeHi = 0U;
  c2_info[55].context = "";
  c2_info[55].name = "sqrt";
  c2_info[55].dominantType = "double";
  c2_info[55].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[55].fileTimeLo = 1343830386U;
  c2_info[55].fileTimeHi = 0U;
  c2_info[55].mFileTimeLo = 0U;
  c2_info[55].mFileTimeHi = 0U;
  c2_info[56].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[56].name = "eml_error";
  c2_info[56].dominantType = "char";
  c2_info[56].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  c2_info[56].fileTimeLo = 1343830358U;
  c2_info[56].fileTimeHi = 0U;
  c2_info[56].mFileTimeLo = 0U;
  c2_info[56].mFileTimeHi = 0U;
  c2_info[57].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c2_info[57].name = "eml_scalar_sqrt";
  c2_info[57].dominantType = "double";
  c2_info[57].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  c2_info[57].fileTimeLo = 1286818738U;
  c2_info[57].fileTimeHi = 0U;
  c2_info[57].mFileTimeLo = 0U;
  c2_info[57].mFileTimeHi = 0U;
  c2_info[58].context = "";
  c2_info[58].name = "mrdivide";
  c2_info[58].dominantType = "double";
  c2_info[58].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c2_info[58].fileTimeLo = 1357951548U;
  c2_info[58].fileTimeHi = 0U;
  c2_info[58].mFileTimeLo = 1319729966U;
  c2_info[58].mFileTimeHi = 0U;
  c2_info[59].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c2_info[59].name = "rdivide";
  c2_info[59].dominantType = "double";
  c2_info[59].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[59].fileTimeLo = 1346510388U;
  c2_info[59].fileTimeHi = 0U;
  c2_info[59].mFileTimeLo = 0U;
  c2_info[59].mFileTimeHi = 0U;
  c2_info[60].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[60].name = "eml_scalexp_compatible";
  c2_info[60].dominantType = "double";
  c2_info[60].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  c2_info[60].fileTimeLo = 1286818796U;
  c2_info[60].fileTimeHi = 0U;
  c2_info[60].mFileTimeLo = 0U;
  c2_info[60].mFileTimeHi = 0U;
  c2_info[61].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c2_info[61].name = "eml_div";
  c2_info[61].dominantType = "double";
  c2_info[61].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c2_info[61].fileTimeLo = 1313347810U;
  c2_info[61].fileTimeHi = 0U;
  c2_info[61].mFileTimeLo = 0U;
  c2_info[61].mFileTimeHi = 0U;
  c2_info[62].context = "";
  c2_info[62].name = "mldivide";
  c2_info[62].dominantType = "double";
  c2_info[62].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c2_info[62].fileTimeLo = 1357951548U;
  c2_info[62].fileTimeHi = 0U;
  c2_info[62].mFileTimeLo = 1319729966U;
  c2_info[62].mFileTimeHi = 0U;
  c2_info[63].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mldivide.p";
  c2_info[63].name = "eml_lusolve";
  c2_info[63].dominantType = "double";
  c2_info[63].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c2_info[63].fileTimeLo = 1309451196U;
  c2_info[63].fileTimeHi = 0U;
  c2_info[63].mFileTimeLo = 0U;
  c2_info[63].mFileTimeHi = 0U;
}

static void c2_b_info_helper(c2_ResolvedFunctionInfo c2_info[159])
{
  c2_info[64].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m";
  c2_info[64].name = "eml_index_class";
  c2_info[64].dominantType = "";
  c2_info[64].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[64].fileTimeLo = 1323170578U;
  c2_info[64].fileTimeHi = 0U;
  c2_info[64].mFileTimeLo = 0U;
  c2_info[64].mFileTimeHi = 0U;
  c2_info[65].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c2_info[65].name = "eml_index_class";
  c2_info[65].dominantType = "";
  c2_info[65].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[65].fileTimeLo = 1323170578U;
  c2_info[65].fileTimeHi = 0U;
  c2_info[65].mFileTimeLo = 0U;
  c2_info[65].mFileTimeHi = 0U;
  c2_info[66].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c2_info[66].name = "eml_xgetrf";
  c2_info[66].dominantType = "double";
  c2_info[66].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c2_info[66].fileTimeLo = 1286818806U;
  c2_info[66].fileTimeHi = 0U;
  c2_info[66].mFileTimeLo = 0U;
  c2_info[66].mFileTimeHi = 0U;
  c2_info[67].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m";
  c2_info[67].name = "eml_lapack_xgetrf";
  c2_info[67].dominantType = "double";
  c2_info[67].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c2_info[67].fileTimeLo = 1286818810U;
  c2_info[67].fileTimeHi = 0U;
  c2_info[67].mFileTimeLo = 0U;
  c2_info[67].mFileTimeHi = 0U;
  c2_info[68].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m";
  c2_info[68].name = "eml_matlab_zgetrf";
  c2_info[68].dominantType = "double";
  c2_info[68].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[68].fileTimeLo = 1302688994U;
  c2_info[68].fileTimeHi = 0U;
  c2_info[68].mFileTimeLo = 0U;
  c2_info[68].mFileTimeHi = 0U;
  c2_info[69].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[69].name = "realmin";
  c2_info[69].dominantType = "char";
  c2_info[69].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m";
  c2_info[69].fileTimeLo = 1307651242U;
  c2_info[69].fileTimeHi = 0U;
  c2_info[69].mFileTimeLo = 0U;
  c2_info[69].mFileTimeHi = 0U;
  c2_info[70].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m";
  c2_info[70].name = "eml_realmin";
  c2_info[70].dominantType = "char";
  c2_info[70].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m";
  c2_info[70].fileTimeLo = 1307651244U;
  c2_info[70].fileTimeHi = 0U;
  c2_info[70].mFileTimeLo = 0U;
  c2_info[70].mFileTimeHi = 0U;
  c2_info[71].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m";
  c2_info[71].name = "eml_float_model";
  c2_info[71].dominantType = "char";
  c2_info[71].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m";
  c2_info[71].fileTimeLo = 1326727996U;
  c2_info[71].fileTimeHi = 0U;
  c2_info[71].mFileTimeLo = 0U;
  c2_info[71].mFileTimeHi = 0U;
  c2_info[72].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[72].name = "eps";
  c2_info[72].dominantType = "char";
  c2_info[72].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c2_info[72].fileTimeLo = 1326727996U;
  c2_info[72].fileTimeHi = 0U;
  c2_info[72].mFileTimeLo = 0U;
  c2_info[72].mFileTimeHi = 0U;
  c2_info[73].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c2_info[73].name = "eml_is_float_class";
  c2_info[73].dominantType = "char";
  c2_info[73].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  c2_info[73].fileTimeLo = 1286818782U;
  c2_info[73].fileTimeHi = 0U;
  c2_info[73].mFileTimeLo = 0U;
  c2_info[73].mFileTimeHi = 0U;
  c2_info[74].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  c2_info[74].name = "eml_eps";
  c2_info[74].dominantType = "char";
  c2_info[74].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m";
  c2_info[74].fileTimeLo = 1326727996U;
  c2_info[74].fileTimeHi = 0U;
  c2_info[74].mFileTimeLo = 0U;
  c2_info[74].mFileTimeHi = 0U;
  c2_info[75].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m";
  c2_info[75].name = "eml_float_model";
  c2_info[75].dominantType = "char";
  c2_info[75].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m";
  c2_info[75].fileTimeLo = 1326727996U;
  c2_info[75].fileTimeHi = 0U;
  c2_info[75].mFileTimeLo = 0U;
  c2_info[75].mFileTimeHi = 0U;
  c2_info[76].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[76].name = "min";
  c2_info[76].dominantType = "coder.internal.indexInt";
  c2_info[76].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c2_info[76].fileTimeLo = 1311255318U;
  c2_info[76].fileTimeHi = 0U;
  c2_info[76].mFileTimeLo = 0U;
  c2_info[76].mFileTimeHi = 0U;
  c2_info[77].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c2_info[77].name = "eml_scalar_eg";
  c2_info[77].dominantType = "coder.internal.indexInt";
  c2_info[77].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[77].fileTimeLo = 1286818796U;
  c2_info[77].fileTimeHi = 0U;
  c2_info[77].mFileTimeLo = 0U;
  c2_info[77].mFileTimeHi = 0U;
  c2_info[78].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  c2_info[78].name = "eml_scalexp_alloc";
  c2_info[78].dominantType = "coder.internal.indexInt";
  c2_info[78].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c2_info[78].fileTimeLo = 1352424860U;
  c2_info[78].fileTimeHi = 0U;
  c2_info[78].mFileTimeLo = 0U;
  c2_info[78].mFileTimeHi = 0U;
  c2_info[79].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum";
  c2_info[79].name = "eml_scalar_eg";
  c2_info[79].dominantType = "coder.internal.indexInt";
  c2_info[79].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[79].fileTimeLo = 1286818796U;
  c2_info[79].fileTimeHi = 0U;
  c2_info[79].mFileTimeLo = 0U;
  c2_info[79].mFileTimeHi = 0U;
  c2_info[80].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[80].name = "colon";
  c2_info[80].dominantType = "double";
  c2_info[80].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c2_info[80].fileTimeLo = 1348191928U;
  c2_info[80].fileTimeHi = 0U;
  c2_info[80].mFileTimeLo = 0U;
  c2_info[80].mFileTimeHi = 0U;
  c2_info[81].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c2_info[81].name = "colon";
  c2_info[81].dominantType = "double";
  c2_info[81].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c2_info[81].fileTimeLo = 1348191928U;
  c2_info[81].fileTimeHi = 0U;
  c2_info[81].mFileTimeLo = 0U;
  c2_info[81].mFileTimeHi = 0U;
  c2_info[82].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  c2_info[82].name = "floor";
  c2_info[82].dominantType = "double";
  c2_info[82].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c2_info[82].fileTimeLo = 1343830380U;
  c2_info[82].fileTimeHi = 0U;
  c2_info[82].mFileTimeLo = 0U;
  c2_info[82].mFileTimeHi = 0U;
  c2_info[83].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c2_info[83].name = "eml_scalar_floor";
  c2_info[83].dominantType = "double";
  c2_info[83].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  c2_info[83].fileTimeLo = 1286818726U;
  c2_info[83].fileTimeHi = 0U;
  c2_info[83].mFileTimeLo = 0U;
  c2_info[83].mFileTimeHi = 0U;
  c2_info[84].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  c2_info[84].name = "intmin";
  c2_info[84].dominantType = "char";
  c2_info[84].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c2_info[84].fileTimeLo = 1311255318U;
  c2_info[84].fileTimeHi = 0U;
  c2_info[84].mFileTimeLo = 0U;
  c2_info[84].mFileTimeHi = 0U;
  c2_info[85].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  c2_info[85].name = "intmax";
  c2_info[85].dominantType = "char";
  c2_info[85].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[85].fileTimeLo = 1311255316U;
  c2_info[85].fileTimeHi = 0U;
  c2_info[85].mFileTimeLo = 0U;
  c2_info[85].mFileTimeHi = 0U;
  c2_info[86].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  c2_info[86].name = "intmin";
  c2_info[86].dominantType = "char";
  c2_info[86].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c2_info[86].fileTimeLo = 1311255318U;
  c2_info[86].fileTimeHi = 0U;
  c2_info[86].mFileTimeLo = 0U;
  c2_info[86].mFileTimeHi = 0U;
  c2_info[87].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  c2_info[87].name = "intmax";
  c2_info[87].dominantType = "char";
  c2_info[87].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[87].fileTimeLo = 1311255316U;
  c2_info[87].fileTimeHi = 0U;
  c2_info[87].mFileTimeLo = 0U;
  c2_info[87].mFileTimeHi = 0U;
  c2_info[88].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  c2_info[88].name = "eml_isa_uint";
  c2_info[88].dominantType = "coder.internal.indexInt";
  c2_info[88].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m";
  c2_info[88].fileTimeLo = 1286818784U;
  c2_info[88].fileTimeHi = 0U;
  c2_info[88].mFileTimeLo = 0U;
  c2_info[88].mFileTimeHi = 0U;
  c2_info[89].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c2_info[89].name = "eml_unsigned_class";
  c2_info[89].dominantType = "char";
  c2_info[89].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  c2_info[89].fileTimeLo = 1323170580U;
  c2_info[89].fileTimeHi = 0U;
  c2_info[89].mFileTimeLo = 0U;
  c2_info[89].mFileTimeHi = 0U;
  c2_info[90].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  c2_info[90].name = "eml_index_class";
  c2_info[90].dominantType = "";
  c2_info[90].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[90].fileTimeLo = 1323170578U;
  c2_info[90].fileTimeHi = 0U;
  c2_info[90].mFileTimeLo = 0U;
  c2_info[90].mFileTimeHi = 0U;
  c2_info[91].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c2_info[91].name = "eml_index_class";
  c2_info[91].dominantType = "";
  c2_info[91].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[91].fileTimeLo = 1323170578U;
  c2_info[91].fileTimeHi = 0U;
  c2_info[91].mFileTimeLo = 0U;
  c2_info[91].mFileTimeHi = 0U;
  c2_info[92].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c2_info[92].name = "intmax";
  c2_info[92].dominantType = "char";
  c2_info[92].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[92].fileTimeLo = 1311255316U;
  c2_info[92].fileTimeHi = 0U;
  c2_info[92].mFileTimeLo = 0U;
  c2_info[92].mFileTimeHi = 0U;
  c2_info[93].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c2_info[93].name = "eml_isa_uint";
  c2_info[93].dominantType = "coder.internal.indexInt";
  c2_info[93].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m";
  c2_info[93].fileTimeLo = 1286818784U;
  c2_info[93].fileTimeHi = 0U;
  c2_info[93].mFileTimeLo = 0U;
  c2_info[93].mFileTimeHi = 0U;
  c2_info[94].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  c2_info[94].name = "eml_index_plus";
  c2_info[94].dominantType = "double";
  c2_info[94].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[94].fileTimeLo = 1286818778U;
  c2_info[94].fileTimeHi = 0U;
  c2_info[94].mFileTimeLo = 0U;
  c2_info[94].mFileTimeHi = 0U;
  c2_info[95].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon";
  c2_info[95].name = "eml_int_forloop_overflow_check";
  c2_info[95].dominantType = "";
  c2_info[95].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[95].fileTimeLo = 1346510340U;
  c2_info[95].fileTimeHi = 0U;
  c2_info[95].mFileTimeLo = 0U;
  c2_info[95].mFileTimeHi = 0U;
  c2_info[96].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[96].name = "eml_index_class";
  c2_info[96].dominantType = "";
  c2_info[96].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[96].fileTimeLo = 1323170578U;
  c2_info[96].fileTimeHi = 0U;
  c2_info[96].mFileTimeLo = 0U;
  c2_info[96].mFileTimeHi = 0U;
  c2_info[97].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[97].name = "eml_index_plus";
  c2_info[97].dominantType = "double";
  c2_info[97].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[97].fileTimeLo = 1286818778U;
  c2_info[97].fileTimeHi = 0U;
  c2_info[97].mFileTimeLo = 0U;
  c2_info[97].mFileTimeHi = 0U;
  c2_info[98].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[98].name = "eml_int_forloop_overflow_check";
  c2_info[98].dominantType = "";
  c2_info[98].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[98].fileTimeLo = 1346510340U;
  c2_info[98].fileTimeHi = 0U;
  c2_info[98].mFileTimeLo = 0U;
  c2_info[98].mFileTimeHi = 0U;
  c2_info[99].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[99].name = "eml_index_minus";
  c2_info[99].dominantType = "double";
  c2_info[99].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[99].fileTimeLo = 1286818778U;
  c2_info[99].fileTimeHi = 0U;
  c2_info[99].mFileTimeLo = 0U;
  c2_info[99].mFileTimeHi = 0U;
  c2_info[100].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[100].name = "eml_index_minus";
  c2_info[100].dominantType = "coder.internal.indexInt";
  c2_info[100].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[100].fileTimeLo = 1286818778U;
  c2_info[100].fileTimeHi = 0U;
  c2_info[100].mFileTimeLo = 0U;
  c2_info[100].mFileTimeHi = 0U;
  c2_info[101].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[101].name = "eml_index_times";
  c2_info[101].dominantType = "coder.internal.indexInt";
  c2_info[101].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c2_info[101].fileTimeLo = 1286818780U;
  c2_info[101].fileTimeHi = 0U;
  c2_info[101].mFileTimeLo = 0U;
  c2_info[101].mFileTimeHi = 0U;
  c2_info[102].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[102].name = "eml_index_plus";
  c2_info[102].dominantType = "coder.internal.indexInt";
  c2_info[102].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[102].fileTimeLo = 1286818778U;
  c2_info[102].fileTimeHi = 0U;
  c2_info[102].mFileTimeLo = 0U;
  c2_info[102].mFileTimeHi = 0U;
  c2_info[103].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[103].name = "eml_ixamax";
  c2_info[103].dominantType = "double";
  c2_info[103].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m";
  c2_info[103].fileTimeLo = 1299076770U;
  c2_info[103].fileTimeHi = 0U;
  c2_info[103].mFileTimeLo = 0U;
  c2_info[103].mFileTimeHi = 0U;
  c2_info[104].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m";
  c2_info[104].name = "eml_blas_inline";
  c2_info[104].dominantType = "";
  c2_info[104].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[104].fileTimeLo = 1299076768U;
  c2_info[104].fileTimeHi = 0U;
  c2_info[104].mFileTimeLo = 0U;
  c2_info[104].mFileTimeHi = 0U;
  c2_info[105].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m!below_threshold";
  c2_info[105].name = "length";
  c2_info[105].dominantType = "double";
  c2_info[105].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  c2_info[105].fileTimeLo = 1303146206U;
  c2_info[105].fileTimeHi = 0U;
  c2_info[105].mFileTimeLo = 0U;
  c2_info[105].mFileTimeHi = 0U;
  c2_info[106].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength";
  c2_info[106].name = "eml_index_class";
  c2_info[106].dominantType = "";
  c2_info[106].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[106].fileTimeLo = 1323170578U;
  c2_info[106].fileTimeHi = 0U;
  c2_info[106].mFileTimeLo = 0U;
  c2_info[106].mFileTimeHi = 0U;
  c2_info[107].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m";
  c2_info[107].name = "eml_index_class";
  c2_info[107].dominantType = "";
  c2_info[107].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[107].fileTimeLo = 1323170578U;
  c2_info[107].fileTimeHi = 0U;
  c2_info[107].mFileTimeLo = 0U;
  c2_info[107].mFileTimeHi = 0U;
  c2_info[108].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m";
  c2_info[108].name = "eml_refblas_ixamax";
  c2_info[108].dominantType = "double";
  c2_info[108].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c2_info[108].fileTimeLo = 1299076770U;
  c2_info[108].fileTimeHi = 0U;
  c2_info[108].mFileTimeLo = 0U;
  c2_info[108].mFileTimeHi = 0U;
  c2_info[109].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c2_info[109].name = "eml_index_class";
  c2_info[109].dominantType = "";
  c2_info[109].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[109].fileTimeLo = 1323170578U;
  c2_info[109].fileTimeHi = 0U;
  c2_info[109].mFileTimeLo = 0U;
  c2_info[109].mFileTimeHi = 0U;
  c2_info[110].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c2_info[110].name = "eml_xcabs1";
  c2_info[110].dominantType = "double";
  c2_info[110].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c2_info[110].fileTimeLo = 1286818706U;
  c2_info[110].fileTimeHi = 0U;
  c2_info[110].mFileTimeLo = 0U;
  c2_info[110].mFileTimeHi = 0U;
  c2_info[111].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m";
  c2_info[111].name = "abs";
  c2_info[111].dominantType = "double";
  c2_info[111].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[111].fileTimeLo = 1343830366U;
  c2_info[111].fileTimeHi = 0U;
  c2_info[111].mFileTimeLo = 0U;
  c2_info[111].mFileTimeHi = 0U;
  c2_info[112].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[112].name = "eml_scalar_abs";
  c2_info[112].dominantType = "double";
  c2_info[112].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c2_info[112].fileTimeLo = 1286818712U;
  c2_info[112].fileTimeHi = 0U;
  c2_info[112].mFileTimeLo = 0U;
  c2_info[112].mFileTimeHi = 0U;
  c2_info[113].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c2_info[113].name = "eml_int_forloop_overflow_check";
  c2_info[113].dominantType = "";
  c2_info[113].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[113].fileTimeLo = 1346510340U;
  c2_info[113].fileTimeHi = 0U;
  c2_info[113].mFileTimeLo = 0U;
  c2_info[113].mFileTimeHi = 0U;
  c2_info[114].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_ixamax.m";
  c2_info[114].name = "eml_index_plus";
  c2_info[114].dominantType = "coder.internal.indexInt";
  c2_info[114].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[114].fileTimeLo = 1286818778U;
  c2_info[114].fileTimeHi = 0U;
  c2_info[114].mFileTimeLo = 0U;
  c2_info[114].mFileTimeHi = 0U;
  c2_info[115].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[115].name = "eml_xswap";
  c2_info[115].dominantType = "double";
  c2_info[115].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m";
  c2_info[115].fileTimeLo = 1299076778U;
  c2_info[115].fileTimeHi = 0U;
  c2_info[115].mFileTimeLo = 0U;
  c2_info[115].mFileTimeHi = 0U;
  c2_info[116].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m";
  c2_info[116].name = "eml_blas_inline";
  c2_info[116].dominantType = "";
  c2_info[116].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[116].fileTimeLo = 1299076768U;
  c2_info[116].fileTimeHi = 0U;
  c2_info[116].mFileTimeLo = 0U;
  c2_info[116].mFileTimeHi = 0U;
  c2_info[117].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m";
  c2_info[117].name = "eml_index_class";
  c2_info[117].dominantType = "";
  c2_info[117].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[117].fileTimeLo = 1323170578U;
  c2_info[117].fileTimeHi = 0U;
  c2_info[117].mFileTimeLo = 0U;
  c2_info[117].mFileTimeHi = 0U;
  c2_info[118].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m";
  c2_info[118].name = "eml_refblas_xswap";
  c2_info[118].dominantType = "double";
  c2_info[118].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c2_info[118].fileTimeLo = 1299076786U;
  c2_info[118].fileTimeHi = 0U;
  c2_info[118].mFileTimeLo = 0U;
  c2_info[118].mFileTimeHi = 0U;
  c2_info[119].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c2_info[119].name = "eml_index_class";
  c2_info[119].dominantType = "";
  c2_info[119].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[119].fileTimeLo = 1323170578U;
  c2_info[119].fileTimeHi = 0U;
  c2_info[119].mFileTimeLo = 0U;
  c2_info[119].mFileTimeHi = 0U;
  c2_info[120].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c2_info[120].name = "abs";
  c2_info[120].dominantType = "coder.internal.indexInt";
  c2_info[120].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[120].fileTimeLo = 1343830366U;
  c2_info[120].fileTimeHi = 0U;
  c2_info[120].mFileTimeLo = 0U;
  c2_info[120].mFileTimeHi = 0U;
  c2_info[121].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[121].name = "eml_scalar_abs";
  c2_info[121].dominantType = "coder.internal.indexInt";
  c2_info[121].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c2_info[121].fileTimeLo = 1286818712U;
  c2_info[121].fileTimeHi = 0U;
  c2_info[121].mFileTimeLo = 0U;
  c2_info[121].mFileTimeHi = 0U;
  c2_info[122].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c2_info[122].name = "eml_int_forloop_overflow_check";
  c2_info[122].dominantType = "";
  c2_info[122].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[122].fileTimeLo = 1346510340U;
  c2_info[122].fileTimeHi = 0U;
  c2_info[122].mFileTimeLo = 0U;
  c2_info[122].mFileTimeHi = 0U;
  c2_info[123].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m";
  c2_info[123].name = "eml_index_plus";
  c2_info[123].dominantType = "coder.internal.indexInt";
  c2_info[123].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[123].fileTimeLo = 1286818778U;
  c2_info[123].fileTimeHi = 0U;
  c2_info[123].mFileTimeLo = 0U;
  c2_info[123].mFileTimeHi = 0U;
  c2_info[124].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[124].name = "eml_div";
  c2_info[124].dominantType = "double";
  c2_info[124].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c2_info[124].fileTimeLo = 1313347810U;
  c2_info[124].fileTimeHi = 0U;
  c2_info[124].mFileTimeLo = 0U;
  c2_info[124].mFileTimeHi = 0U;
  c2_info[125].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m";
  c2_info[125].name = "eml_xgeru";
  c2_info[125].dominantType = "double";
  c2_info[125].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c2_info[125].fileTimeLo = 1299076774U;
  c2_info[125].fileTimeHi = 0U;
  c2_info[125].mFileTimeLo = 0U;
  c2_info[125].mFileTimeHi = 0U;
  c2_info[126].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c2_info[126].name = "eml_blas_inline";
  c2_info[126].dominantType = "";
  c2_info[126].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[126].fileTimeLo = 1299076768U;
  c2_info[126].fileTimeHi = 0U;
  c2_info[126].mFileTimeLo = 0U;
  c2_info[126].mFileTimeHi = 0U;
  c2_info[127].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgeru.m";
  c2_info[127].name = "eml_xger";
  c2_info[127].dominantType = "double";
  c2_info[127].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m";
  c2_info[127].fileTimeLo = 1299076774U;
  c2_info[127].fileTimeHi = 0U;
  c2_info[127].mFileTimeLo = 0U;
  c2_info[127].mFileTimeHi = 0U;
}

static void c2_c_info_helper(c2_ResolvedFunctionInfo c2_info[159])
{
  c2_info[128].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xger.m";
  c2_info[128].name = "eml_blas_inline";
  c2_info[128].dominantType = "";
  c2_info[128].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[128].fileTimeLo = 1299076768U;
  c2_info[128].fileTimeHi = 0U;
  c2_info[128].mFileTimeLo = 0U;
  c2_info[128].mFileTimeHi = 0U;
  c2_info[129].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold";
  c2_info[129].name = "intmax";
  c2_info[129].dominantType = "char";
  c2_info[129].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c2_info[129].fileTimeLo = 1311255316U;
  c2_info[129].fileTimeHi = 0U;
  c2_info[129].mFileTimeLo = 0U;
  c2_info[129].mFileTimeHi = 0U;
  c2_info[130].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold";
  c2_info[130].name = "min";
  c2_info[130].dominantType = "double";
  c2_info[130].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  c2_info[130].fileTimeLo = 1311255318U;
  c2_info[130].fileTimeHi = 0U;
  c2_info[130].mFileTimeLo = 0U;
  c2_info[130].mFileTimeHi = 0U;
  c2_info[131].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m!below_threshold";
  c2_info[131].name = "mtimes";
  c2_info[131].dominantType = "double";
  c2_info[131].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[131].fileTimeLo = 1289519692U;
  c2_info[131].fileTimeHi = 0U;
  c2_info[131].mFileTimeLo = 0U;
  c2_info[131].mFileTimeHi = 0U;
  c2_info[132].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m";
  c2_info[132].name = "eml_index_class";
  c2_info[132].dominantType = "";
  c2_info[132].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[132].fileTimeLo = 1323170578U;
  c2_info[132].fileTimeHi = 0U;
  c2_info[132].mFileTimeLo = 0U;
  c2_info[132].mFileTimeHi = 0U;
  c2_info[133].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xger.m";
  c2_info[133].name = "eml_refblas_xger";
  c2_info[133].dominantType = "double";
  c2_info[133].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m";
  c2_info[133].fileTimeLo = 1299076776U;
  c2_info[133].fileTimeHi = 0U;
  c2_info[133].mFileTimeLo = 0U;
  c2_info[133].mFileTimeHi = 0U;
  c2_info[134].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xger.m";
  c2_info[134].name = "eml_refblas_xgerx";
  c2_info[134].dominantType = "char";
  c2_info[134].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[134].fileTimeLo = 1299076778U;
  c2_info[134].fileTimeHi = 0U;
  c2_info[134].mFileTimeLo = 0U;
  c2_info[134].mFileTimeHi = 0U;
  c2_info[135].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[135].name = "eml_index_class";
  c2_info[135].dominantType = "";
  c2_info[135].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[135].fileTimeLo = 1323170578U;
  c2_info[135].fileTimeHi = 0U;
  c2_info[135].mFileTimeLo = 0U;
  c2_info[135].mFileTimeHi = 0U;
  c2_info[136].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[136].name = "abs";
  c2_info[136].dominantType = "coder.internal.indexInt";
  c2_info[136].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c2_info[136].fileTimeLo = 1343830366U;
  c2_info[136].fileTimeHi = 0U;
  c2_info[136].mFileTimeLo = 0U;
  c2_info[136].mFileTimeHi = 0U;
  c2_info[137].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[137].name = "eml_index_minus";
  c2_info[137].dominantType = "double";
  c2_info[137].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[137].fileTimeLo = 1286818778U;
  c2_info[137].fileTimeHi = 0U;
  c2_info[137].mFileTimeLo = 0U;
  c2_info[137].mFileTimeHi = 0U;
  c2_info[138].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[138].name = "eml_int_forloop_overflow_check";
  c2_info[138].dominantType = "";
  c2_info[138].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[138].fileTimeLo = 1346510340U;
  c2_info[138].fileTimeHi = 0U;
  c2_info[138].mFileTimeLo = 0U;
  c2_info[138].mFileTimeHi = 0U;
  c2_info[139].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[139].name = "eml_index_plus";
  c2_info[139].dominantType = "double";
  c2_info[139].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[139].fileTimeLo = 1286818778U;
  c2_info[139].fileTimeHi = 0U;
  c2_info[139].mFileTimeLo = 0U;
  c2_info[139].mFileTimeHi = 0U;
  c2_info[140].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgerx.m";
  c2_info[140].name = "eml_index_plus";
  c2_info[140].dominantType = "coder.internal.indexInt";
  c2_info[140].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[140].fileTimeLo = 1286818778U;
  c2_info[140].fileTimeHi = 0U;
  c2_info[140].mFileTimeLo = 0U;
  c2_info[140].mFileTimeHi = 0U;
  c2_info[141].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!warn_singular";
  c2_info[141].name = "eml_warning";
  c2_info[141].dominantType = "char";
  c2_info[141].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m";
  c2_info[141].fileTimeLo = 1286818802U;
  c2_info[141].fileTimeHi = 0U;
  c2_info[141].mFileTimeLo = 0U;
  c2_info[141].mFileTimeHi = 0U;
  c2_info[142].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c2_info[142].name = "eml_scalar_eg";
  c2_info[142].dominantType = "double";
  c2_info[142].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[142].fileTimeLo = 1286818796U;
  c2_info[142].fileTimeHi = 0U;
  c2_info[142].mFileTimeLo = 0U;
  c2_info[142].mFileTimeHi = 0U;
  c2_info[143].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c2_info[143].name = "eml_int_forloop_overflow_check";
  c2_info[143].dominantType = "";
  c2_info[143].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[143].fileTimeLo = 1346510340U;
  c2_info[143].fileTimeHi = 0U;
  c2_info[143].mFileTimeLo = 0U;
  c2_info[143].mFileTimeHi = 0U;
  c2_info[144].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_lusolve.m!lusolveNxN";
  c2_info[144].name = "eml_xtrsm";
  c2_info[144].dominantType = "char";
  c2_info[144].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m";
  c2_info[144].fileTimeLo = 1299076778U;
  c2_info[144].fileTimeHi = 0U;
  c2_info[144].mFileTimeLo = 0U;
  c2_info[144].mFileTimeHi = 0U;
  c2_info[145].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m";
  c2_info[145].name = "eml_blas_inline";
  c2_info[145].dominantType = "";
  c2_info[145].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c2_info[145].fileTimeLo = 1299076768U;
  c2_info[145].fileTimeHi = 0U;
  c2_info[145].mFileTimeLo = 0U;
  c2_info[145].mFileTimeHi = 0U;
  c2_info[146].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m!below_threshold";
  c2_info[146].name = "mtimes";
  c2_info[146].dominantType = "double";
  c2_info[146].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c2_info[146].fileTimeLo = 1289519692U;
  c2_info[146].fileTimeHi = 0U;
  c2_info[146].mFileTimeLo = 0U;
  c2_info[146].mFileTimeHi = 0U;
  c2_info[147].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c2_info[147].name = "eml_index_class";
  c2_info[147].dominantType = "";
  c2_info[147].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[147].fileTimeLo = 1323170578U;
  c2_info[147].fileTimeHi = 0U;
  c2_info[147].mFileTimeLo = 0U;
  c2_info[147].mFileTimeHi = 0U;
  c2_info[148].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c2_info[148].name = "eml_scalar_eg";
  c2_info[148].dominantType = "double";
  c2_info[148].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[148].fileTimeLo = 1286818796U;
  c2_info[148].fileTimeHi = 0U;
  c2_info[148].mFileTimeLo = 0U;
  c2_info[148].mFileTimeHi = 0U;
  c2_info[149].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m";
  c2_info[149].name = "eml_refblas_xtrsm";
  c2_info[149].dominantType = "char";
  c2_info[149].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[149].fileTimeLo = 1299076786U;
  c2_info[149].fileTimeHi = 0U;
  c2_info[149].mFileTimeLo = 0U;
  c2_info[149].mFileTimeHi = 0U;
  c2_info[150].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[150].name = "eml_scalar_eg";
  c2_info[150].dominantType = "double";
  c2_info[150].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c2_info[150].fileTimeLo = 1286818796U;
  c2_info[150].fileTimeHi = 0U;
  c2_info[150].mFileTimeLo = 0U;
  c2_info[150].mFileTimeHi = 0U;
  c2_info[151].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[151].name = "eml_index_minus";
  c2_info[151].dominantType = "double";
  c2_info[151].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m";
  c2_info[151].fileTimeLo = 1286818778U;
  c2_info[151].fileTimeHi = 0U;
  c2_info[151].mFileTimeLo = 0U;
  c2_info[151].mFileTimeHi = 0U;
  c2_info[152].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[152].name = "eml_index_class";
  c2_info[152].dominantType = "";
  c2_info[152].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c2_info[152].fileTimeLo = 1323170578U;
  c2_info[152].fileTimeHi = 0U;
  c2_info[152].mFileTimeLo = 0U;
  c2_info[152].mFileTimeHi = 0U;
  c2_info[153].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[153].name = "eml_index_times";
  c2_info[153].dominantType = "coder.internal.indexInt";
  c2_info[153].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  c2_info[153].fileTimeLo = 1286818780U;
  c2_info[153].fileTimeHi = 0U;
  c2_info[153].mFileTimeLo = 0U;
  c2_info[153].mFileTimeHi = 0U;
  c2_info[154].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[154].name = "eml_index_plus";
  c2_info[154].dominantType = "coder.internal.indexInt";
  c2_info[154].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[154].fileTimeLo = 1286818778U;
  c2_info[154].fileTimeHi = 0U;
  c2_info[154].mFileTimeLo = 0U;
  c2_info[154].mFileTimeHi = 0U;
  c2_info[155].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[155].name = "eml_int_forloop_overflow_check";
  c2_info[155].dominantType = "";
  c2_info[155].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c2_info[155].fileTimeLo = 1346510340U;
  c2_info[155].fileTimeHi = 0U;
  c2_info[155].mFileTimeLo = 0U;
  c2_info[155].mFileTimeHi = 0U;
  c2_info[156].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[156].name = "eml_index_plus";
  c2_info[156].dominantType = "double";
  c2_info[156].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  c2_info[156].fileTimeLo = 1286818778U;
  c2_info[156].fileTimeHi = 0U;
  c2_info[156].mFileTimeLo = 0U;
  c2_info[156].mFileTimeHi = 0U;
  c2_info[157].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  c2_info[157].name = "intmin";
  c2_info[157].dominantType = "char";
  c2_info[157].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  c2_info[157].fileTimeLo = 1311255318U;
  c2_info[157].fileTimeHi = 0U;
  c2_info[157].mFileTimeLo = 0U;
  c2_info[157].mFileTimeHi = 0U;
  c2_info[158].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xtrsm.m";
  c2_info[158].name = "eml_div";
  c2_info[158].dominantType = "double";
  c2_info[158].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c2_info[158].fileTimeLo = 1313347810U;
  c2_info[158].fileTimeHi = 0U;
  c2_info[158].mFileTimeLo = 0U;
  c2_info[158].mFileTimeHi = 0U;
}

static void c2_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_b_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance)
{
}

static real_T c2_sum(SFc2_testInstanceStruct *chartInstance, real_T c2_x[6])
{
  real_T c2_y;
  int32_T c2_k;
  int32_T c2_b_k;
  c2_y = c2_x[0];
  for (c2_k = 2; c2_k < 7; c2_k++) {
    c2_b_k = c2_k;
    c2_y += c2_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c2_b_k), 1, 6, 1, 0) - 1];
  }

  return c2_y;
}

static void c2_c_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_d_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_realmin(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_eps(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_eml_matlab_zgetrf(SFc2_testInstanceStruct *chartInstance, real_T
  c2_A[36], real_T c2_b_A[36], int32_T c2_ipiv[6], int32_T *c2_info)
{
  int32_T c2_i170;
  for (c2_i170 = 0; c2_i170 < 36; c2_i170++) {
    c2_b_A[c2_i170] = c2_A[c2_i170];
  }

  c2_b_eml_matlab_zgetrf(chartInstance, c2_b_A, c2_ipiv, c2_info);
}

static void c2_check_forloop_overflow_error(SFc2_testInstanceStruct
  *chartInstance, boolean_T c2_overflow)
{
  int32_T c2_i171;
  static char_T c2_cv0[34] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'i', 'n', 't', '_', 'f', 'o', 'r', 'l', 'o', 'o', 'p',
    '_', 'o', 'v', 'e', 'r', 'f', 'l', 'o', 'w' };

  char_T c2_u[34];
  const mxArray *c2_y = NULL;
  int32_T c2_i172;
  static char_T c2_cv1[23] = { 'c', 'o', 'd', 'e', 'r', '.', 'i', 'n', 't', 'e',
    'r', 'n', 'a', 'l', '.', 'i', 'n', 'd', 'e', 'x', 'I', 'n', 't' };

  char_T c2_b_u[23];
  const mxArray *c2_b_y = NULL;
  if (!c2_overflow) {
  } else {
    for (c2_i171 = 0; c2_i171 < 34; c2_i171++) {
      c2_u[c2_i171] = c2_cv0[c2_i171];
    }

    c2_y = NULL;
    sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 34),
                  FALSE);
    for (c2_i172 = 0; c2_i172 < 23; c2_i172++) {
      c2_b_u[c2_i172] = c2_cv1[c2_i172];
    }

    c2_b_y = NULL;
    sf_mex_assign(&c2_b_y, sf_mex_create("y", c2_b_u, 10, 0U, 1U, 0U, 2, 1, 23),
                  FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U,
      14, c2_y, 14, c2_b_y));
  }
}

static void c2_eml_xger(SFc2_testInstanceStruct *chartInstance, int32_T c2_m,
  int32_T c2_n, real_T c2_alpha1, int32_T c2_ix0, int32_T c2_iy0, real_T c2_A[36],
  int32_T c2_ia0, real_T c2_b_A[36])
{
  int32_T c2_i173;
  for (c2_i173 = 0; c2_i173 < 36; c2_i173++) {
    c2_b_A[c2_i173] = c2_A[c2_i173];
  }

  c2_b_eml_xger(chartInstance, c2_m, c2_n, c2_alpha1, c2_ix0, c2_iy0, c2_b_A,
                c2_ia0);
}

static void c2_eml_warning(SFc2_testInstanceStruct *chartInstance)
{
  int32_T c2_i174;
  static char_T c2_varargin_1[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A',
    'T', 'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a',
    't', 'r', 'i', 'x' };

  char_T c2_u[27];
  const mxArray *c2_y = NULL;
  for (c2_i174 = 0; c2_i174 < 27; c2_i174++) {
    c2_u[c2_i174] = c2_varargin_1[c2_i174];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 27), FALSE);
  sf_mex_call_debug("warning", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U,
    14, c2_y));
}

static void c2_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A[36],
  real_T c2_B[6], real_T c2_b_B[6])
{
  int32_T c2_i175;
  int32_T c2_i176;
  real_T c2_b_A[36];
  for (c2_i175 = 0; c2_i175 < 6; c2_i175++) {
    c2_b_B[c2_i175] = c2_B[c2_i175];
  }

  for (c2_i176 = 0; c2_i176 < 36; c2_i176++) {
    c2_b_A[c2_i176] = c2_A[c2_i176];
  }

  c2_c_eml_xtrsm(chartInstance, c2_b_A, c2_b_B);
}

static void c2_below_threshold(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_e_eml_scalar_eg(SFc2_testInstanceStruct *chartInstance)
{
}

static void c2_b_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A
  [36], real_T c2_B[6], real_T c2_b_B[6])
{
  int32_T c2_i177;
  int32_T c2_i178;
  real_T c2_b_A[36];
  for (c2_i177 = 0; c2_i177 < 6; c2_i177++) {
    c2_b_B[c2_i177] = c2_B[c2_i177];
  }

  for (c2_i178 = 0; c2_i178 < 36; c2_i178++) {
    c2_b_A[c2_i178] = c2_A[c2_i178];
  }

  c2_d_eml_xtrsm(chartInstance, c2_b_A, c2_b_B);
}

static const mxArray *c2_h_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(int32_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static int32_T c2_i_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_y;
  int32_T c2_i179;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_i179, 1, 6, 0U, 0, 0U, 0);
  c2_y = c2_i179;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_b_sfEvent;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_y;
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)chartInstanceVoid;
  c2_b_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_sfEvent),
    &c2_thisId);
  sf_mex_destroy(&c2_b_sfEvent);
  *(int32_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_j_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_b_is_active_c2_test, const char_T *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_b_is_active_c2_test),
    &c2_thisId);
  sf_mex_destroy(&c2_b_is_active_c2_test);
  return c2_y;
}

static uint8_T c2_k_emlrt_marshallIn(SFc2_testInstanceStruct *chartInstance,
  const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_b_eml_matlab_zgetrf(SFc2_testInstanceStruct *chartInstance,
  real_T c2_A[36], int32_T c2_ipiv[6], int32_T *c2_info)
{
  int32_T c2_i180;
  int32_T c2_j;
  int32_T c2_b_j;
  int32_T c2_a;
  int32_T c2_jm1;
  int32_T c2_b;
  int32_T c2_mmj;
  int32_T c2_b_a;
  int32_T c2_c;
  int32_T c2_b_b;
  int32_T c2_jj;
  int32_T c2_c_a;
  int32_T c2_jp1j;
  int32_T c2_d_a;
  int32_T c2_b_c;
  int32_T c2_n;
  int32_T c2_ix0;
  int32_T c2_b_n;
  int32_T c2_b_ix0;
  int32_T c2_c_n;
  int32_T c2_c_ix0;
  int32_T c2_idxmax;
  int32_T c2_ix;
  real_T c2_x;
  real_T c2_b_x;
  real_T c2_c_x;
  real_T c2_y;
  real_T c2_d_x;
  real_T c2_e_x;
  real_T c2_b_y;
  real_T c2_smax;
  int32_T c2_d_n;
  int32_T c2_c_b;
  int32_T c2_d_b;
  boolean_T c2_overflow;
  int32_T c2_k;
  int32_T c2_b_k;
  int32_T c2_e_a;
  real_T c2_f_x;
  real_T c2_g_x;
  real_T c2_h_x;
  real_T c2_c_y;
  real_T c2_i_x;
  real_T c2_j_x;
  real_T c2_d_y;
  real_T c2_s;
  int32_T c2_f_a;
  int32_T c2_jpiv_offset;
  int32_T c2_g_a;
  int32_T c2_e_b;
  int32_T c2_jpiv;
  int32_T c2_h_a;
  int32_T c2_f_b;
  int32_T c2_c_c;
  int32_T c2_g_b;
  int32_T c2_jrow;
  int32_T c2_i_a;
  int32_T c2_h_b;
  int32_T c2_jprow;
  int32_T c2_d_ix0;
  int32_T c2_iy0;
  int32_T c2_e_ix0;
  int32_T c2_b_iy0;
  int32_T c2_f_ix0;
  int32_T c2_c_iy0;
  int32_T c2_b_ix;
  int32_T c2_iy;
  int32_T c2_c_k;
  real_T c2_temp;
  int32_T c2_j_a;
  int32_T c2_k_a;
  int32_T c2_b_jp1j;
  int32_T c2_l_a;
  int32_T c2_d_c;
  int32_T c2_m_a;
  int32_T c2_i_b;
  int32_T c2_i181;
  int32_T c2_n_a;
  int32_T c2_j_b;
  int32_T c2_o_a;
  int32_T c2_k_b;
  boolean_T c2_b_overflow;
  int32_T c2_i;
  int32_T c2_b_i;
  real_T c2_k_x;
  real_T c2_e_y;
  real_T c2_z;
  int32_T c2_l_b;
  int32_T c2_e_c;
  int32_T c2_p_a;
  int32_T c2_f_c;
  int32_T c2_q_a;
  int32_T c2_g_c;
  int32_T c2_m;
  int32_T c2_e_n;
  int32_T c2_g_ix0;
  int32_T c2_d_iy0;
  int32_T c2_ia0;
  real_T c2_d1;
  c2_realmin(chartInstance);
  c2_eps(chartInstance);
  for (c2_i180 = 0; c2_i180 < 6; c2_i180++) {
    c2_ipiv[c2_i180] = 1 + c2_i180;
  }

  *c2_info = 0;
  for (c2_j = 1; c2_j < 6; c2_j++) {
    c2_b_j = c2_j;
    c2_a = c2_b_j - 1;
    c2_jm1 = c2_a;
    c2_b = c2_b_j;
    c2_mmj = 6 - c2_b;
    c2_b_a = c2_jm1;
    c2_c = c2_b_a * 7;
    c2_b_b = c2_c + 1;
    c2_jj = c2_b_b;
    c2_c_a = c2_jj + 1;
    c2_jp1j = c2_c_a;
    c2_d_a = c2_mmj;
    c2_b_c = c2_d_a;
    c2_n = c2_b_c + 1;
    c2_ix0 = c2_jj;
    c2_b_n = c2_n;
    c2_b_ix0 = c2_ix0;
    c2_c_n = c2_b_n;
    c2_c_ix0 = c2_b_ix0;
    if (c2_c_n < 1) {
      c2_idxmax = 0;
    } else {
      c2_idxmax = 1;
      if (c2_c_n > 1) {
        c2_ix = c2_c_ix0;
        c2_x = c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
          "", (real_T)c2_ix), 1, 36, 1, 0) - 1];
        c2_b_x = c2_x;
        c2_c_x = c2_b_x;
        c2_y = muDoubleScalarAbs(c2_c_x);
        c2_d_x = 0.0;
        c2_e_x = c2_d_x;
        c2_b_y = muDoubleScalarAbs(c2_e_x);
        c2_smax = c2_y + c2_b_y;
        c2_d_n = c2_c_n;
        c2_c_b = c2_d_n;
        c2_d_b = c2_c_b;
        if (2 > c2_d_b) {
          c2_overflow = FALSE;
        } else {
          c2_overflow = (c2_d_b > 2147483646);
        }

        if (c2_overflow) {
          c2_check_forloop_overflow_error(chartInstance, c2_overflow);
        }

        for (c2_k = 2; c2_k <= c2_d_n; c2_k++) {
          c2_b_k = c2_k;
          c2_e_a = c2_ix + 1;
          c2_ix = c2_e_a;
          c2_f_x = c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c2_ix), 1, 36, 1, 0) - 1];
          c2_g_x = c2_f_x;
          c2_h_x = c2_g_x;
          c2_c_y = muDoubleScalarAbs(c2_h_x);
          c2_i_x = 0.0;
          c2_j_x = c2_i_x;
          c2_d_y = muDoubleScalarAbs(c2_j_x);
          c2_s = c2_c_y + c2_d_y;
          if (c2_s > c2_smax) {
            c2_idxmax = c2_b_k;
            c2_smax = c2_s;
          }
        }
      }
    }

    c2_f_a = c2_idxmax - 1;
    c2_jpiv_offset = c2_f_a;
    c2_g_a = c2_jj;
    c2_e_b = c2_jpiv_offset;
    c2_jpiv = c2_g_a + c2_e_b;
    if (c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_jpiv), 1, 36, 1, 0) - 1] != 0.0) {
      if (c2_jpiv_offset != 0) {
        c2_h_a = c2_b_j;
        c2_f_b = c2_jpiv_offset;
        c2_c_c = c2_h_a + c2_f_b;
        c2_ipiv[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_b_j), 1, 6, 1, 0) - 1] = c2_c_c;
        c2_g_b = c2_jm1 + 1;
        c2_jrow = c2_g_b;
        c2_i_a = c2_jrow;
        c2_h_b = c2_jpiv_offset;
        c2_jprow = c2_i_a + c2_h_b;
        c2_d_ix0 = c2_jrow;
        c2_iy0 = c2_jprow;
        c2_e_ix0 = c2_d_ix0;
        c2_b_iy0 = c2_iy0;
        c2_f_ix0 = c2_e_ix0;
        c2_c_iy0 = c2_b_iy0;
        c2_b_ix = c2_f_ix0;
        c2_iy = c2_c_iy0;
        for (c2_c_k = 1; c2_c_k < 7; c2_c_k++) {
          c2_temp = c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
            _SFD_INTEGER_CHECK("", (real_T)c2_b_ix), 1, 36, 1, 0) - 1];
          c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c2_b_ix), 1, 36, 1, 0) - 1] =
            c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c2_iy), 1, 36, 1, 0) - 1];
          c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c2_iy), 1, 36, 1, 0) - 1] = c2_temp;
          c2_j_a = c2_b_ix + 6;
          c2_b_ix = c2_j_a;
          c2_k_a = c2_iy + 6;
          c2_iy = c2_k_a;
        }
      }

      c2_b_jp1j = c2_jp1j;
      c2_l_a = c2_mmj;
      c2_d_c = c2_l_a;
      c2_m_a = c2_jp1j;
      c2_i_b = c2_d_c - 1;
      c2_i181 = c2_m_a + c2_i_b;
      c2_n_a = c2_b_jp1j;
      c2_j_b = c2_i181;
      c2_o_a = c2_n_a;
      c2_k_b = c2_j_b;
      if (c2_o_a > c2_k_b) {
        c2_b_overflow = FALSE;
      } else {
        c2_b_overflow = (c2_k_b > 2147483646);
      }

      if (c2_b_overflow) {
        c2_check_forloop_overflow_error(chartInstance, c2_b_overflow);
      }

      for (c2_i = c2_b_jp1j; c2_i <= c2_i181; c2_i++) {
        c2_b_i = c2_i;
        c2_k_x = c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c2_b_i), 1, 36, 1, 0) - 1];
        c2_e_y = c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c2_jj), 1, 36, 1, 0) - 1];
        c2_z = c2_k_x / c2_e_y;
        c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_b_i), 1, 36, 1, 0) - 1] = c2_z;
      }
    } else {
      *c2_info = c2_b_j;
    }

    c2_l_b = c2_b_j;
    c2_e_c = 6 - c2_l_b;
    c2_p_a = c2_jj;
    c2_f_c = c2_p_a;
    c2_q_a = c2_jj;
    c2_g_c = c2_q_a;
    c2_m = c2_mmj;
    c2_e_n = c2_e_c;
    c2_g_ix0 = c2_jp1j;
    c2_d_iy0 = c2_f_c + 6;
    c2_ia0 = c2_g_c + 7;
    c2_d1 = -1.0;
    c2_b_eml_xger(chartInstance, c2_m, c2_e_n, c2_d1, c2_g_ix0, c2_d_iy0, c2_A,
                  c2_ia0);
  }

  if (*c2_info == 0) {
    if (!(c2_A[35] != 0.0)) {
      *c2_info = 6;
    }
  }
}

static void c2_b_eml_xger(SFc2_testInstanceStruct *chartInstance, int32_T c2_m,
  int32_T c2_n, real_T c2_alpha1, int32_T c2_ix0, int32_T c2_iy0, real_T c2_A[36],
  int32_T c2_ia0)
{
  int32_T c2_b_m;
  int32_T c2_b_n;
  real_T c2_b_alpha1;
  int32_T c2_b_ix0;
  int32_T c2_b_iy0;
  int32_T c2_b_ia0;
  int32_T c2_c_m;
  int32_T c2_c_n;
  real_T c2_c_alpha1;
  int32_T c2_c_ix0;
  int32_T c2_c_iy0;
  int32_T c2_c_ia0;
  int32_T c2_d_m;
  int32_T c2_d_n;
  real_T c2_d_alpha1;
  int32_T c2_d_ix0;
  int32_T c2_d_iy0;
  int32_T c2_d_ia0;
  int32_T c2_ixstart;
  int32_T c2_a;
  int32_T c2_jA;
  int32_T c2_jy;
  int32_T c2_e_n;
  int32_T c2_b;
  int32_T c2_b_b;
  boolean_T c2_overflow;
  int32_T c2_j;
  real_T c2_yjy;
  real_T c2_temp;
  int32_T c2_ix;
  int32_T c2_c_b;
  int32_T c2_i182;
  int32_T c2_b_a;
  int32_T c2_d_b;
  int32_T c2_i183;
  int32_T c2_c_a;
  int32_T c2_e_b;
  int32_T c2_d_a;
  int32_T c2_f_b;
  boolean_T c2_b_overflow;
  int32_T c2_ijA;
  int32_T c2_b_ijA;
  int32_T c2_e_a;
  int32_T c2_f_a;
  int32_T c2_g_a;
  c2_b_m = c2_m;
  c2_b_n = c2_n;
  c2_b_alpha1 = c2_alpha1;
  c2_b_ix0 = c2_ix0;
  c2_b_iy0 = c2_iy0;
  c2_b_ia0 = c2_ia0;
  c2_c_m = c2_b_m;
  c2_c_n = c2_b_n;
  c2_c_alpha1 = c2_b_alpha1;
  c2_c_ix0 = c2_b_ix0;
  c2_c_iy0 = c2_b_iy0;
  c2_c_ia0 = c2_b_ia0;
  c2_d_m = c2_c_m;
  c2_d_n = c2_c_n;
  c2_d_alpha1 = c2_c_alpha1;
  c2_d_ix0 = c2_c_ix0;
  c2_d_iy0 = c2_c_iy0;
  c2_d_ia0 = c2_c_ia0;
  if (c2_d_alpha1 == 0.0) {
  } else {
    c2_ixstart = c2_d_ix0;
    c2_a = c2_d_ia0 - 1;
    c2_jA = c2_a;
    c2_jy = c2_d_iy0;
    c2_e_n = c2_d_n;
    c2_b = c2_e_n;
    c2_b_b = c2_b;
    if (1 > c2_b_b) {
      c2_overflow = FALSE;
    } else {
      c2_overflow = (c2_b_b > 2147483646);
    }

    if (c2_overflow) {
      c2_check_forloop_overflow_error(chartInstance, c2_overflow);
    }

    for (c2_j = 1; c2_j <= c2_e_n; c2_j++) {
      c2_yjy = c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK(
        "", (real_T)c2_jy), 1, 36, 1, 0) - 1];
      if (c2_yjy != 0.0) {
        c2_temp = c2_yjy * c2_d_alpha1;
        c2_ix = c2_ixstart;
        c2_c_b = c2_jA + 1;
        c2_i182 = c2_c_b;
        c2_b_a = c2_d_m;
        c2_d_b = c2_jA;
        c2_i183 = c2_b_a + c2_d_b;
        c2_c_a = c2_i182;
        c2_e_b = c2_i183;
        c2_d_a = c2_c_a;
        c2_f_b = c2_e_b;
        if (c2_d_a > c2_f_b) {
          c2_b_overflow = FALSE;
        } else {
          c2_b_overflow = (c2_f_b > 2147483646);
        }

        if (c2_b_overflow) {
          c2_check_forloop_overflow_error(chartInstance, c2_b_overflow);
        }

        for (c2_ijA = c2_i182; c2_ijA <= c2_i183; c2_ijA++) {
          c2_b_ijA = c2_ijA;
          c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c2_b_ijA), 1, 36, 1, 0) - 1] =
            c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c2_b_ijA), 1, 36, 1, 0) - 1] +
            c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
            (real_T)c2_ix), 1, 36, 1, 0) - 1] * c2_temp;
          c2_e_a = c2_ix + 1;
          c2_ix = c2_e_a;
        }
      }

      c2_f_a = c2_jy + 6;
      c2_jy = c2_f_a;
      c2_g_a = c2_jA + 6;
      c2_jA = c2_g_a;
    }
  }
}

static void c2_c_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A
  [36], real_T c2_B[6])
{
  int32_T c2_k;
  int32_T c2_b_k;
  int32_T c2_a;
  int32_T c2_c;
  int32_T c2_b;
  int32_T c2_b_c;
  int32_T c2_b_b;
  int32_T c2_kAcol;
  int32_T c2_b_a;
  int32_T c2_c_c;
  int32_T c2_c_a;
  int32_T c2_i184;
  boolean_T c2_overflow;
  int32_T c2_i;
  int32_T c2_b_i;
  int32_T c2_d_a;
  int32_T c2_d_c;
  int32_T c2_e_a;
  int32_T c2_e_c;
  int32_T c2_f_a;
  int32_T c2_f_c;
  int32_T c2_g_a;
  int32_T c2_c_b;
  int32_T c2_g_c;
  c2_below_threshold(chartInstance);
  c2_e_eml_scalar_eg(chartInstance);
  for (c2_k = 1; c2_k < 7; c2_k++) {
    c2_b_k = c2_k;
    c2_a = c2_b_k;
    c2_c = c2_a;
    c2_b = c2_c - 1;
    c2_b_c = 6 * c2_b;
    c2_b_b = c2_b_c;
    c2_kAcol = c2_b_b;
    c2_b_a = c2_b_k;
    c2_c_c = c2_b_a;
    if (c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_c_c), 1, 6, 1, 0) - 1] != 0.0) {
      c2_c_a = c2_b_k;
      c2_i184 = c2_c_a;
      c2_overflow = FALSE;
      if (c2_overflow) {
        c2_check_forloop_overflow_error(chartInstance, c2_overflow);
      }

      for (c2_i = c2_i184 + 1; c2_i < 7; c2_i++) {
        c2_b_i = c2_i;
        c2_d_a = c2_b_i;
        c2_d_c = c2_d_a;
        c2_e_a = c2_b_i;
        c2_e_c = c2_e_a;
        c2_f_a = c2_b_k;
        c2_f_c = c2_f_a;
        c2_g_a = c2_b_i;
        c2_c_b = c2_kAcol;
        c2_g_c = c2_g_a + c2_c_b;
        c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_d_c), 1, 6, 1, 0) - 1] = c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c2_e_c), 1, 6, 1, 0) - 1]
          - c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_f_c), 1, 6, 1, 0) - 1] * c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c2_g_c), 1, 36, 1, 0) - 1];
      }
    }
  }
}

static void c2_d_eml_xtrsm(SFc2_testInstanceStruct *chartInstance, real_T c2_A
  [36], real_T c2_B[6])
{
  int32_T c2_k;
  int32_T c2_b_k;
  int32_T c2_a;
  int32_T c2_c;
  int32_T c2_b;
  int32_T c2_b_c;
  int32_T c2_b_b;
  int32_T c2_kAcol;
  int32_T c2_b_a;
  int32_T c2_c_c;
  int32_T c2_c_a;
  int32_T c2_d_c;
  int32_T c2_d_a;
  int32_T c2_e_c;
  int32_T c2_e_a;
  int32_T c2_c_b;
  int32_T c2_f_c;
  real_T c2_x;
  real_T c2_y;
  real_T c2_z;
  int32_T c2_f_a;
  int32_T c2_i185;
  int32_T c2_d_b;
  int32_T c2_e_b;
  boolean_T c2_overflow;
  int32_T c2_i;
  int32_T c2_b_i;
  int32_T c2_g_a;
  int32_T c2_g_c;
  int32_T c2_h_a;
  int32_T c2_h_c;
  int32_T c2_i_a;
  int32_T c2_i_c;
  int32_T c2_j_a;
  int32_T c2_f_b;
  int32_T c2_j_c;
  c2_below_threshold(chartInstance);
  c2_e_eml_scalar_eg(chartInstance);
  for (c2_k = 6; c2_k > 0; c2_k--) {
    c2_b_k = c2_k;
    c2_a = c2_b_k;
    c2_c = c2_a;
    c2_b = c2_c - 1;
    c2_b_c = 6 * c2_b;
    c2_b_b = c2_b_c;
    c2_kAcol = c2_b_b;
    c2_b_a = c2_b_k;
    c2_c_c = c2_b_a;
    if (c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_c_c), 1, 6, 1, 0) - 1] != 0.0) {
      c2_c_a = c2_b_k;
      c2_d_c = c2_c_a;
      c2_d_a = c2_b_k;
      c2_e_c = c2_d_a;
      c2_e_a = c2_b_k;
      c2_c_b = c2_kAcol;
      c2_f_c = c2_e_a + c2_c_b;
      c2_x = c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c2_e_c), 1, 6, 1, 0) - 1];
      c2_y = c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c2_f_c), 1, 36, 1, 0) - 1];
      c2_z = c2_x / c2_y;
      c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
        (real_T)c2_d_c), 1, 6, 1, 0) - 1] = c2_z;
      c2_f_a = c2_b_k - 1;
      c2_i185 = c2_f_a;
      c2_d_b = c2_i185;
      c2_e_b = c2_d_b;
      if (1 > c2_e_b) {
        c2_overflow = FALSE;
      } else {
        c2_overflow = (c2_e_b > 2147483646);
      }

      if (c2_overflow) {
        c2_check_forloop_overflow_error(chartInstance, c2_overflow);
      }

      for (c2_i = 1; c2_i <= c2_i185; c2_i++) {
        c2_b_i = c2_i;
        c2_g_a = c2_b_i;
        c2_g_c = c2_g_a;
        c2_h_a = c2_b_i;
        c2_h_c = c2_h_a;
        c2_i_a = c2_b_k;
        c2_i_c = c2_i_a;
        c2_j_a = c2_b_i;
        c2_f_b = c2_kAcol;
        c2_j_c = c2_j_a + c2_f_b;
        c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_g_c), 1, 6, 1, 0) - 1] = c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c2_h_c), 1, 6, 1, 0) - 1]
          - c2_B[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c2_i_c), 1, 6, 1, 0) - 1] * c2_A[_SFD_EML_ARRAY_BOUNDS_CHECK(
          "", (int32_T)_SFD_INTEGER_CHECK("", (real_T)c2_j_c), 1, 36, 1, 0) - 1];
      }
    }
  }
}

static void init_dsm_address_info(SFc2_testInstanceStruct *chartInstance)
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

void sf_c2_test_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2481613192U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1533748073U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3087715314U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2333820233U);
}

mxArray *sf_c2_test_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("JNEdJMP8LPuP1mhFamrCWH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,11,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(6);
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
      pr[0] = (double)(6);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(3);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,7,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,7,"type",mxType);
    }

    mxSetField(mxData,7,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,8,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,8,"type",mxType);
    }

    mxSetField(mxData,8,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,9,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,9,"type",mxType);
    }

    mxSetField(mxData,9,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,10,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,10,"type",mxType);
    }

    mxSetField(mxData,10,"complexity",mxCreateDoubleScalar(0));
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

mxArray *sf_c2_test_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c2_test(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"nu_dot\",},{M[8],M[0],T\"is_active_c2_test\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_test_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_testInstanceStruct *chartInstance;
    chartInstance = (SFc2_testInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _testMachineNumber_,
           2,
           1,
           1,
           12,
           0,
           0,
           0,
           0,
           2,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"nu");
          _SFD_SET_DATA_PROPS(1,2,0,1,"nu_dot");
          _SFD_SET_DATA_PROPS(2,1,1,0,"eta");
          _SFD_SET_DATA_PROPS(3,1,1,0,"omega");
          _SFD_SET_DATA_PROPS(4,1,1,0,"g");
          _SFD_SET_DATA_PROPS(5,1,1,0,"m");
          _SFD_SET_DATA_PROPS(6,1,1,0,"I");
          _SFD_SET_DATA_PROPS(7,1,1,0,"b");
          _SFD_SET_DATA_PROPS(8,1,1,0,"k");
          _SFD_SET_DATA_PROPS(9,1,1,0,"l");
          _SFD_SET_DATA_PROPS(10,1,1,0,"wind");
          _SFD_SET_DATA_PROPS(11,1,1,0,"windDir");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,713);
        _SFD_CV_INIT_SCRIPT(0,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"m2c",0,-1,1704);
        _SFD_CV_INIT_SCRIPT(1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(1,0,"Smtrx",0,-1,1271);
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
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)
            c2_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c2_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c2_g;
          real_T *c2_m;
          real_T *c2_b;
          real_T *c2_k;
          real_T *c2_l;
          real_T *c2_windDir;
          real_T (*c2_nu)[6];
          real_T (*c2_nu_dot)[6];
          real_T (*c2_eta)[6];
          real_T (*c2_omega)[6];
          real_T (*c2_I)[9];
          real_T (*c2_wind)[2];
          c2_windDir = (real_T *)ssGetInputPortSignal(chartInstance->S, 10);
          c2_wind = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 9);
          c2_l = (real_T *)ssGetInputPortSignal(chartInstance->S, 8);
          c2_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 7);
          c2_b = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
          c2_I = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 5);
          c2_m = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c2_g = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c2_omega = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 2);
          c2_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
          c2_nu_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
          c2_nu = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c2_nu);
          _SFD_SET_DATA_VALUE_PTR(1U, *c2_nu_dot);
          _SFD_SET_DATA_VALUE_PTR(2U, *c2_eta);
          _SFD_SET_DATA_VALUE_PTR(3U, *c2_omega);
          _SFD_SET_DATA_VALUE_PTR(4U, c2_g);
          _SFD_SET_DATA_VALUE_PTR(5U, c2_m);
          _SFD_SET_DATA_VALUE_PTR(6U, *c2_I);
          _SFD_SET_DATA_VALUE_PTR(7U, c2_b);
          _SFD_SET_DATA_VALUE_PTR(8U, c2_k);
          _SFD_SET_DATA_VALUE_PTR(9U, c2_l);
          _SFD_SET_DATA_VALUE_PTR(10U, *c2_wind);
          _SFD_SET_DATA_VALUE_PTR(11U, c2_windDir);
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
  return "oA5wojXtBtQH0p1JJKjStD";
}

static void sf_opaque_initialize_c2_test(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc2_testInstanceStruct*) chartInstanceVar)->S,0);
  initialize_params_c2_test((SFc2_testInstanceStruct*) chartInstanceVar);
  initialize_c2_test((SFc2_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c2_test(void *chartInstanceVar)
{
  enable_c2_test((SFc2_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c2_test(void *chartInstanceVar)
{
  disable_c2_test((SFc2_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c2_test(void *chartInstanceVar)
{
  sf_c2_test((SFc2_testInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c2_test(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c2_test((SFc2_testInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_test();/* state var info */
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

extern void sf_internal_set_sim_state_c2_test(SimStruct* S, const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_test();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c2_test((SFc2_testInstanceStruct*)chartInfo->chartInstance,
                        mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c2_test(SimStruct* S)
{
  return sf_internal_get_sim_state_c2_test(S);
}

static void sf_opaque_set_sim_state_c2_test(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c2_test(S, st);
}

static void sf_opaque_terminate_c2_test(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_testInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_optimization_info();
    }

    finalize_c2_test((SFc2_testInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_test((SFc2_testInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_test(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c2_test((SFc2_testInstanceStruct*)(((ChartInfoStruct *)
      ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c2_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,2,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,2);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 7, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 8, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 9, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 10, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,11);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 11; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1069826812U));
  ssSetChecksum1(S,(4203885011U));
  ssSetChecksum2(S,(3376131596U));
  ssSetChecksum3(S,(2909028122U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_test(SimStruct *S)
{
  SFc2_testInstanceStruct *chartInstance;
  chartInstance = (SFc2_testInstanceStruct *)utMalloc(sizeof
    (SFc2_testInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc2_testInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c2_test;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c2_test;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c2_test;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c2_test;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c2_test;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c2_test;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c2_test;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c2_test;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_test;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_test;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c2_test;
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

void c2_test_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_test(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_test(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_test(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_test_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
