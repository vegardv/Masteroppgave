/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_sfun.h"
#include "c4_test.h"
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
static const char * c4_debug_family_names[16] = { "psi", "T", "B_inv",
  "Theta_ref", "nargin", "nargout", "eta", "eta_dot", "pos_ref", "m", "K_p_pos",
  "K_d_pos", "omega", "k", "phi_ref", "theta_ref" };

/* Function Declarations */
static void initialize_c4_test(SFc4_testInstanceStruct *chartInstance);
static void initialize_params_c4_test(SFc4_testInstanceStruct *chartInstance);
static void enable_c4_test(SFc4_testInstanceStruct *chartInstance);
static void disable_c4_test(SFc4_testInstanceStruct *chartInstance);
static void c4_update_debugger_state_c4_test(SFc4_testInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c4_test(SFc4_testInstanceStruct
  *chartInstance);
static void set_sim_state_c4_test(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_st);
static void finalize_c4_test(SFc4_testInstanceStruct *chartInstance);
static void sf_c4_test(SFc4_testInstanceStruct *chartInstance);
static void c4_chartstep_c4_test(SFc4_testInstanceStruct *chartInstance);
static void initSimStructsc4_test(SFc4_testInstanceStruct *chartInstance);
static void registerMessagesc4_test(SFc4_testInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber);
static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData);
static real_T c4_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_theta_ref, const char_T *c4_identifier);
static real_T c4_b_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static const mxArray *c4_e_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static void c4_c_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[2]);
static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static void c4_d_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[4]);
static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static void c4_info_helper(c4_ResolvedFunctionInfo c4_info[25]);
static real_T c4_sum(SFc4_testInstanceStruct *chartInstance, real_T c4_x[6]);
static void c4_eml_scalar_eg(SFc4_testInstanceStruct *chartInstance);
static const mxArray *c4_f_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData);
static int32_T c4_e_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData);
static uint8_T c4_f_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_b_is_active_c4_test, const char_T *c4_identifier);
static uint8_T c4_g_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId);
static void init_dsm_address_info(SFc4_testInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c4_test(SFc4_testInstanceStruct *chartInstance)
{
  chartInstance->c4_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c4_is_active_c4_test = 0U;
}

static void initialize_params_c4_test(SFc4_testInstanceStruct *chartInstance)
{
}

static void enable_c4_test(SFc4_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c4_test(SFc4_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c4_update_debugger_state_c4_test(SFc4_testInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c4_test(SFc4_testInstanceStruct
  *chartInstance)
{
  const mxArray *c4_st;
  const mxArray *c4_y = NULL;
  real_T c4_hoistedGlobal;
  real_T c4_u;
  const mxArray *c4_b_y = NULL;
  real_T c4_b_hoistedGlobal;
  real_T c4_b_u;
  const mxArray *c4_c_y = NULL;
  uint8_T c4_c_hoistedGlobal;
  uint8_T c4_c_u;
  const mxArray *c4_d_y = NULL;
  real_T *c4_phi_ref;
  real_T *c4_theta_ref;
  c4_theta_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c4_phi_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c4_st = NULL;
  c4_st = NULL;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_createcellarray(3), FALSE);
  c4_hoistedGlobal = *c4_phi_ref;
  c4_u = c4_hoistedGlobal;
  c4_b_y = NULL;
  sf_mex_assign(&c4_b_y, sf_mex_create("y", &c4_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c4_y, 0, c4_b_y);
  c4_b_hoistedGlobal = *c4_theta_ref;
  c4_b_u = c4_b_hoistedGlobal;
  c4_c_y = NULL;
  sf_mex_assign(&c4_c_y, sf_mex_create("y", &c4_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c4_y, 1, c4_c_y);
  c4_c_hoistedGlobal = chartInstance->c4_is_active_c4_test;
  c4_c_u = c4_c_hoistedGlobal;
  c4_d_y = NULL;
  sf_mex_assign(&c4_d_y, sf_mex_create("y", &c4_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c4_y, 2, c4_d_y);
  sf_mex_assign(&c4_st, c4_y, FALSE);
  return c4_st;
}

static void set_sim_state_c4_test(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_st)
{
  const mxArray *c4_u;
  real_T *c4_phi_ref;
  real_T *c4_theta_ref;
  c4_theta_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c4_phi_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c4_doneDoubleBufferReInit = TRUE;
  c4_u = sf_mex_dup(c4_st);
  *c4_phi_ref = c4_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c4_u, 0)), "phi_ref");
  *c4_theta_ref = c4_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c4_u, 1)), "theta_ref");
  chartInstance->c4_is_active_c4_test = c4_f_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c4_u, 2)), "is_active_c4_test");
  sf_mex_destroy(&c4_u);
  c4_update_debugger_state_c4_test(chartInstance);
  sf_mex_destroy(&c4_st);
}

static void finalize_c4_test(SFc4_testInstanceStruct *chartInstance)
{
}

static void sf_c4_test(SFc4_testInstanceStruct *chartInstance)
{
  int32_T c4_i0;
  int32_T c4_i1;
  int32_T c4_i2;
  int32_T c4_i3;
  int32_T c4_i4;
  int32_T c4_i5;
  real_T *c4_phi_ref;
  real_T *c4_theta_ref;
  real_T *c4_m;
  real_T *c4_k;
  real_T (*c4_omega)[6];
  real_T (*c4_K_d_pos)[4];
  real_T (*c4_K_p_pos)[4];
  real_T (*c4_pos_ref)[3];
  real_T (*c4_eta_dot)[6];
  real_T (*c4_eta)[6];
  c4_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 7);
  c4_omega = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 6);
  c4_K_d_pos = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 5);
  c4_K_p_pos = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 4);
  c4_m = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c4_theta_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c4_pos_ref = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
  c4_phi_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c4_eta_dot = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c4_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
  for (c4_i0 = 0; c4_i0 < 6; c4_i0++) {
    _SFD_DATA_RANGE_CHECK((*c4_eta)[c4_i0], 0U);
  }

  for (c4_i1 = 0; c4_i1 < 6; c4_i1++) {
    _SFD_DATA_RANGE_CHECK((*c4_eta_dot)[c4_i1], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c4_phi_ref, 2U);
  for (c4_i2 = 0; c4_i2 < 3; c4_i2++) {
    _SFD_DATA_RANGE_CHECK((*c4_pos_ref)[c4_i2], 3U);
  }

  _SFD_DATA_RANGE_CHECK(*c4_theta_ref, 4U);
  _SFD_DATA_RANGE_CHECK(*c4_m, 5U);
  for (c4_i3 = 0; c4_i3 < 4; c4_i3++) {
    _SFD_DATA_RANGE_CHECK((*c4_K_p_pos)[c4_i3], 6U);
  }

  for (c4_i4 = 0; c4_i4 < 4; c4_i4++) {
    _SFD_DATA_RANGE_CHECK((*c4_K_d_pos)[c4_i4], 7U);
  }

  for (c4_i5 = 0; c4_i5 < 6; c4_i5++) {
    _SFD_DATA_RANGE_CHECK((*c4_omega)[c4_i5], 8U);
  }

  _SFD_DATA_RANGE_CHECK(*c4_k, 9U);
  chartInstance->c4_sfEvent = CALL_EVENT;
  c4_chartstep_c4_test(chartInstance);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_testMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c4_chartstep_c4_test(SFc4_testInstanceStruct *chartInstance)
{
  real_T c4_hoistedGlobal;
  real_T c4_b_hoistedGlobal;
  int32_T c4_i6;
  real_T c4_eta[6];
  int32_T c4_i7;
  real_T c4_eta_dot[6];
  int32_T c4_i8;
  real_T c4_pos_ref[3];
  real_T c4_m;
  int32_T c4_i9;
  real_T c4_K_p_pos[4];
  int32_T c4_i10;
  real_T c4_K_d_pos[4];
  int32_T c4_i11;
  real_T c4_omega[6];
  real_T c4_k;
  uint32_T c4_debug_family_var_map[16];
  real_T c4_psi;
  real_T c4_T;
  real_T c4_B_inv[4];
  real_T c4_Theta_ref[2];
  real_T c4_nargin = 8.0;
  real_T c4_nargout = 2.0;
  real_T c4_phi_ref;
  real_T c4_theta_ref;
  real_T c4_a;
  int32_T c4_i12;
  real_T c4_b_omega[6];
  real_T c4_b;
  real_T c4_A;
  real_T c4_B;
  real_T c4_x;
  real_T c4_y;
  real_T c4_b_x;
  real_T c4_b_y;
  real_T c4_c_y;
  real_T c4_c_x;
  real_T c4_d_x;
  real_T c4_e_x;
  real_T c4_f_x;
  real_T c4_g_x;
  real_T c4_h_x;
  real_T c4_i_x;
  real_T c4_j_x;
  real_T c4_b_a;
  real_T c4_b_b[4];
  int32_T c4_i13;
  int32_T c4_i14;
  int32_T c4_i15;
  real_T c4_c_b[2];
  int32_T c4_i16;
  real_T c4_d_y[2];
  int32_T c4_i17;
  int32_T c4_i18;
  int32_T c4_i19;
  int32_T c4_i20;
  int32_T c4_i21;
  real_T c4_e_y[2];
  int32_T c4_i22;
  int32_T c4_i23;
  int32_T c4_i24;
  int32_T c4_i25;
  int32_T c4_i26;
  real_T c4_f_y[2];
  int32_T c4_i27;
  int32_T c4_i28;
  int32_T c4_i29;
  int32_T c4_i30;
  int32_T c4_i31;
  int32_T c4_i32;
  int32_T c4_i33;
  int32_T c4_i34;
  int32_T c4_i35;
  int32_T c4_i36;
  int32_T c4_i37;
  int32_T c4_i38;
  int32_T c4_i39;
  real_T *c4_b_m;
  real_T *c4_b_k;
  real_T *c4_b_phi_ref;
  real_T *c4_b_theta_ref;
  real_T (*c4_c_omega)[6];
  real_T (*c4_b_K_d_pos)[4];
  real_T (*c4_b_K_p_pos)[4];
  real_T (*c4_b_pos_ref)[3];
  real_T (*c4_b_eta_dot)[6];
  real_T (*c4_b_eta)[6];
  c4_b_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 7);
  c4_c_omega = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 6);
  c4_b_K_d_pos = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 5);
  c4_b_K_p_pos = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 4);
  c4_b_m = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c4_b_theta_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c4_b_pos_ref = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
  c4_b_phi_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c4_b_eta_dot = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c4_b_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
  c4_hoistedGlobal = *c4_b_m;
  c4_b_hoistedGlobal = *c4_b_k;
  for (c4_i6 = 0; c4_i6 < 6; c4_i6++) {
    c4_eta[c4_i6] = (*c4_b_eta)[c4_i6];
  }

  for (c4_i7 = 0; c4_i7 < 6; c4_i7++) {
    c4_eta_dot[c4_i7] = (*c4_b_eta_dot)[c4_i7];
  }

  for (c4_i8 = 0; c4_i8 < 3; c4_i8++) {
    c4_pos_ref[c4_i8] = (*c4_b_pos_ref)[c4_i8];
  }

  c4_m = c4_hoistedGlobal;
  for (c4_i9 = 0; c4_i9 < 4; c4_i9++) {
    c4_K_p_pos[c4_i9] = (*c4_b_K_p_pos)[c4_i9];
  }

  for (c4_i10 = 0; c4_i10 < 4; c4_i10++) {
    c4_K_d_pos[c4_i10] = (*c4_b_K_d_pos)[c4_i10];
  }

  for (c4_i11 = 0; c4_i11 < 6; c4_i11++) {
    c4_omega[c4_i11] = (*c4_c_omega)[c4_i11];
  }

  c4_k = c4_b_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 16U, 16U, c4_debug_family_names,
    c4_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_psi, 0U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_T, 1U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_B_inv, 2U, c4_c_sf_marshallOut,
    c4_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c4_Theta_ref, 3U, c4_e_sf_marshallOut,
    c4_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargin, 4U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_nargout, 5U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_eta, 6U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_eta_dot, 7U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_pos_ref, 8U, c4_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_m, 9U, c4_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_K_p_pos, 10U, c4_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_K_d_pos, 11U, c4_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c4_omega, 12U, c4_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c4_k, 13U, c4_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_phi_ref, 14U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c4_theta_ref, 15U, c4_sf_marshallOut,
    c4_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 4);
  c4_psi = c4_eta[5];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 5);
  c4_a = c4_k;
  for (c4_i12 = 0; c4_i12 < 6; c4_i12++) {
    c4_b_omega[c4_i12] = c4_omega[c4_i12];
  }

  c4_b = c4_sum(chartInstance, c4_b_omega);
  c4_T = c4_a * c4_b;
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 6);
  if (CV_EML_IF(0, 1, 0, c4_T == 0.0)) {
    _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 7);
    c4_T = 0.1;
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 9);
  c4_A = c4_m;
  c4_B = c4_T;
  c4_x = c4_A;
  c4_y = c4_B;
  c4_b_x = c4_x;
  c4_b_y = c4_y;
  c4_c_y = c4_b_x / c4_b_y;
  c4_c_x = c4_psi;
  c4_d_x = c4_c_x;
  c4_d_x = muDoubleScalarSin(c4_d_x);
  c4_e_x = c4_psi;
  c4_f_x = c4_e_x;
  c4_f_x = muDoubleScalarCos(c4_f_x);
  c4_g_x = c4_psi;
  c4_h_x = c4_g_x;
  c4_h_x = muDoubleScalarCos(c4_h_x);
  c4_i_x = c4_psi;
  c4_j_x = c4_i_x;
  c4_j_x = muDoubleScalarSin(c4_j_x);
  c4_b_a = c4_c_y;
  c4_b_b[0] = -c4_d_x;
  c4_b_b[2] = c4_f_x;
  c4_b_b[1] = -c4_h_x;
  c4_b_b[3] = -c4_j_x;
  for (c4_i13 = 0; c4_i13 < 4; c4_i13++) {
    c4_B_inv[c4_i13] = c4_b_a * c4_b_b[c4_i13];
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 11);
  for (c4_i14 = 0; c4_i14 < 4; c4_i14++) {
    c4_b_b[c4_i14] = -c4_K_d_pos[c4_i14];
  }

  for (c4_i15 = 0; c4_i15 < 2; c4_i15++) {
    c4_c_b[c4_i15] = c4_eta_dot[c4_i15];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i16 = 0; c4_i16 < 2; c4_i16++) {
    c4_d_y[c4_i16] = 0.0;
    c4_i17 = 0;
    for (c4_i18 = 0; c4_i18 < 2; c4_i18++) {
      c4_d_y[c4_i16] += c4_b_b[c4_i17 + c4_i16] * c4_c_b[c4_i18];
      c4_i17 += 2;
    }
  }

  for (c4_i19 = 0; c4_i19 < 4; c4_i19++) {
    c4_b_b[c4_i19] = c4_K_p_pos[c4_i19];
  }

  for (c4_i20 = 0; c4_i20 < 2; c4_i20++) {
    c4_c_b[c4_i20] = c4_eta[c4_i20];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i21 = 0; c4_i21 < 2; c4_i21++) {
    c4_e_y[c4_i21] = 0.0;
    c4_i22 = 0;
    for (c4_i23 = 0; c4_i23 < 2; c4_i23++) {
      c4_e_y[c4_i21] += c4_b_b[c4_i22 + c4_i21] * c4_c_b[c4_i23];
      c4_i22 += 2;
    }
  }

  for (c4_i24 = 0; c4_i24 < 4; c4_i24++) {
    c4_b_b[c4_i24] = c4_K_p_pos[c4_i24];
  }

  for (c4_i25 = 0; c4_i25 < 2; c4_i25++) {
    c4_c_b[c4_i25] = c4_pos_ref[c4_i25];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i26 = 0; c4_i26 < 2; c4_i26++) {
    c4_f_y[c4_i26] = 0.0;
    c4_i27 = 0;
    for (c4_i28 = 0; c4_i28 < 2; c4_i28++) {
      c4_f_y[c4_i26] += c4_b_b[c4_i27 + c4_i26] * c4_c_b[c4_i28];
      c4_i27 += 2;
    }
  }

  for (c4_i29 = 0; c4_i29 < 4; c4_i29++) {
    c4_b_b[c4_i29] = c4_B_inv[c4_i29];
  }

  for (c4_i30 = 0; c4_i30 < 2; c4_i30++) {
    c4_d_y[c4_i30] = (c4_d_y[c4_i30] - c4_e_y[c4_i30]) + c4_f_y[c4_i30];
  }

  c4_eml_scalar_eg(chartInstance);
  c4_eml_scalar_eg(chartInstance);
  for (c4_i31 = 0; c4_i31 < 2; c4_i31++) {
    c4_Theta_ref[c4_i31] = 0.0;
  }

  for (c4_i32 = 0; c4_i32 < 2; c4_i32++) {
    c4_Theta_ref[c4_i32] = 0.0;
  }

  for (c4_i33 = 0; c4_i33 < 2; c4_i33++) {
    c4_c_b[c4_i33] = c4_Theta_ref[c4_i33];
  }

  for (c4_i34 = 0; c4_i34 < 2; c4_i34++) {
    c4_Theta_ref[c4_i34] = c4_c_b[c4_i34];
  }

  for (c4_i35 = 0; c4_i35 < 2; c4_i35++) {
    c4_c_b[c4_i35] = c4_Theta_ref[c4_i35];
  }

  for (c4_i36 = 0; c4_i36 < 2; c4_i36++) {
    c4_Theta_ref[c4_i36] = c4_c_b[c4_i36];
  }

  for (c4_i37 = 0; c4_i37 < 2; c4_i37++) {
    c4_Theta_ref[c4_i37] = 0.0;
    c4_i38 = 0;
    for (c4_i39 = 0; c4_i39 < 2; c4_i39++) {
      c4_Theta_ref[c4_i37] += c4_b_b[c4_i38 + c4_i37] * c4_d_y[c4_i39];
      c4_i38 += 2;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 12);
  c4_theta_ref = c4_Theta_ref[1];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, 13);
  c4_phi_ref = c4_Theta_ref[0];
  _SFD_EML_CALL(0U, chartInstance->c4_sfEvent, -13);
  _SFD_SYMBOL_SCOPE_POP();
  *c4_b_phi_ref = c4_phi_ref;
  *c4_b_theta_ref = c4_theta_ref;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U, chartInstance->c4_sfEvent);
}

static void initSimStructsc4_test(SFc4_testInstanceStruct *chartInstance)
{
}

static void registerMessagesc4_test(SFc4_testInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c4_machineNumber, uint32_T
  c4_chartNumber)
{
}

static const mxArray *c4_sf_marshallOut(void *chartInstanceVoid, void *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  real_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(real_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static real_T c4_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_theta_ref, const char_T *c4_identifier)
{
  real_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_theta_ref),
    &c4_thisId);
  sf_mex_destroy(&c4_theta_ref);
  return c4_y;
}

static real_T c4_b_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  real_T c4_y;
  real_T c4_d0;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_d0, 1, 0, 0U, 0, 0U, 0);
  c4_y = c4_d0;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_theta_ref;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_theta_ref = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_theta_ref),
    &c4_thisId);
  sf_mex_destroy(&c4_theta_ref);
  *(real_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static const mxArray *c4_b_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i40;
  real_T c4_b_inData[6];
  int32_T c4_i41;
  real_T c4_u[6];
  const mxArray *c4_y = NULL;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i40 = 0; c4_i40 < 6; c4_i40++) {
    c4_b_inData[c4_i40] = (*(real_T (*)[6])c4_inData)[c4_i40];
  }

  for (c4_i41 = 0; c4_i41 < 6; c4_i41++) {
    c4_u[c4_i41] = c4_b_inData[c4_i41];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static const mxArray *c4_c_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i42;
  int32_T c4_i43;
  int32_T c4_i44;
  real_T c4_b_inData[4];
  int32_T c4_i45;
  int32_T c4_i46;
  int32_T c4_i47;
  real_T c4_u[4];
  const mxArray *c4_y = NULL;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_i42 = 0;
  for (c4_i43 = 0; c4_i43 < 2; c4_i43++) {
    for (c4_i44 = 0; c4_i44 < 2; c4_i44++) {
      c4_b_inData[c4_i44 + c4_i42] = (*(real_T (*)[4])c4_inData)[c4_i44 + c4_i42];
    }

    c4_i42 += 2;
  }

  c4_i45 = 0;
  for (c4_i46 = 0; c4_i46 < 2; c4_i46++) {
    for (c4_i47 = 0; c4_i47 < 2; c4_i47++) {
      c4_u[c4_i47 + c4_i45] = c4_b_inData[c4_i47 + c4_i45];
    }

    c4_i45 += 2;
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 2, 2, 2), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static const mxArray *c4_d_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i48;
  real_T c4_b_inData[3];
  int32_T c4_i49;
  real_T c4_u[3];
  const mxArray *c4_y = NULL;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i48 = 0; c4_i48 < 3; c4_i48++) {
    c4_b_inData[c4_i48] = (*(real_T (*)[3])c4_inData)[c4_i48];
  }

  for (c4_i49 = 0; c4_i49 < 3; c4_i49++) {
    c4_u[c4_i49] = c4_b_inData[c4_i49];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static const mxArray *c4_e_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_i50;
  real_T c4_b_inData[2];
  int32_T c4_i51;
  real_T c4_u[2];
  const mxArray *c4_y = NULL;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  for (c4_i50 = 0; c4_i50 < 2; c4_i50++) {
    c4_b_inData[c4_i50] = (*(real_T (*)[2])c4_inData)[c4_i50];
  }

  for (c4_i51 = 0; c4_i51 < 2; c4_i51++) {
    c4_u[c4_i51] = c4_b_inData[c4_i51];
  }

  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", c4_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static void c4_c_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[2])
{
  real_T c4_dv0[2];
  int32_T c4_i52;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv0, 1, 0, 0U, 1, 0U, 1, 2);
  for (c4_i52 = 0; c4_i52 < 2; c4_i52++) {
    c4_y[c4_i52] = c4_dv0[c4_i52];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_Theta_ref;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[2];
  int32_T c4_i53;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_Theta_ref = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_Theta_ref), &c4_thisId,
                        c4_y);
  sf_mex_destroy(&c4_Theta_ref);
  for (c4_i53 = 0; c4_i53 < 2; c4_i53++) {
    (*(real_T (*)[2])c4_outData)[c4_i53] = c4_y[c4_i53];
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

static void c4_d_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance, const
  mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId, real_T c4_y[4])
{
  real_T c4_dv1[4];
  int32_T c4_i54;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), c4_dv1, 1, 0, 0U, 1, 0U, 2, 2, 2);
  for (c4_i54 = 0; c4_i54 < 4; c4_i54++) {
    c4_y[c4_i54] = c4_dv1[c4_i54];
  }

  sf_mex_destroy(&c4_u);
}

static void c4_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_B_inv;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  real_T c4_y[4];
  int32_T c4_i55;
  int32_T c4_i56;
  int32_T c4_i57;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_B_inv = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_B_inv), &c4_thisId, c4_y);
  sf_mex_destroy(&c4_B_inv);
  c4_i55 = 0;
  for (c4_i56 = 0; c4_i56 < 2; c4_i56++) {
    for (c4_i57 = 0; c4_i57 < 2; c4_i57++) {
      (*(real_T (*)[4])c4_outData)[c4_i57 + c4_i55] = c4_y[c4_i57 + c4_i55];
    }

    c4_i55 += 2;
  }

  sf_mex_destroy(&c4_mxArrayInData);
}

const mxArray *sf_c4_test_get_eml_resolved_functions_info(void)
{
  const mxArray *c4_nameCaptureInfo;
  c4_ResolvedFunctionInfo c4_info[25];
  const mxArray *c4_m0 = NULL;
  int32_T c4_i58;
  c4_ResolvedFunctionInfo *c4_r0;
  c4_nameCaptureInfo = NULL;
  c4_nameCaptureInfo = NULL;
  c4_info_helper(c4_info);
  sf_mex_assign(&c4_m0, sf_mex_createstruct("nameCaptureInfo", 1, 25), FALSE);
  for (c4_i58 = 0; c4_i58 < 25; c4_i58++) {
    c4_r0 = &c4_info[c4_i58];
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", c4_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c4_r0->context)), "context", "nameCaptureInfo",
                    c4_i58);
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", c4_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c4_r0->name)), "name", "nameCaptureInfo", c4_i58);
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", c4_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c4_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c4_i58);
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", c4_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c4_r0->resolved)), "resolved", "nameCaptureInfo",
                    c4_i58);
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", &c4_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c4_i58);
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", &c4_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c4_i58);
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", &c4_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c4_i58);
    sf_mex_addfield(c4_m0, sf_mex_create("nameCaptureInfo", &c4_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c4_i58);
  }

  sf_mex_assign(&c4_nameCaptureInfo, c4_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c4_nameCaptureInfo);
  return c4_nameCaptureInfo;
}

static void c4_info_helper(c4_ResolvedFunctionInfo c4_info[25])
{
  c4_info[0].context = "";
  c4_info[0].name = "sum";
  c4_info[0].dominantType = "double";
  c4_info[0].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c4_info[0].fileTimeLo = 1314736612U;
  c4_info[0].fileTimeHi = 0U;
  c4_info[0].mFileTimeLo = 0U;
  c4_info[0].mFileTimeHi = 0U;
  c4_info[1].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c4_info[1].name = "isequal";
  c4_info[1].dominantType = "double";
  c4_info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c4_info[1].fileTimeLo = 1286818758U;
  c4_info[1].fileTimeHi = 0U;
  c4_info[1].mFileTimeLo = 0U;
  c4_info[1].mFileTimeHi = 0U;
  c4_info[2].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  c4_info[2].name = "eml_isequal_core";
  c4_info[2].dominantType = "double";
  c4_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  c4_info[2].fileTimeLo = 1286818786U;
  c4_info[2].fileTimeHi = 0U;
  c4_info[2].mFileTimeLo = 0U;
  c4_info[2].mFileTimeHi = 0U;
  c4_info[3].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c4_info[3].name = "eml_const_nonsingleton_dim";
  c4_info[3].dominantType = "double";
  c4_info[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_const_nonsingleton_dim.m";
  c4_info[3].fileTimeLo = 1286818696U;
  c4_info[3].fileTimeHi = 0U;
  c4_info[3].mFileTimeLo = 0U;
  c4_info[3].mFileTimeHi = 0U;
  c4_info[4].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c4_info[4].name = "eml_scalar_eg";
  c4_info[4].dominantType = "double";
  c4_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c4_info[4].fileTimeLo = 1286818796U;
  c4_info[4].fileTimeHi = 0U;
  c4_info[4].mFileTimeLo = 0U;
  c4_info[4].mFileTimeHi = 0U;
  c4_info[5].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c4_info[5].name = "eml_index_class";
  c4_info[5].dominantType = "";
  c4_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c4_info[5].fileTimeLo = 1323170578U;
  c4_info[5].fileTimeHi = 0U;
  c4_info[5].mFileTimeLo = 0U;
  c4_info[5].mFileTimeHi = 0U;
  c4_info[6].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/sum.m";
  c4_info[6].name = "eml_int_forloop_overflow_check";
  c4_info[6].dominantType = "";
  c4_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  c4_info[6].fileTimeLo = 1346510340U;
  c4_info[6].fileTimeHi = 0U;
  c4_info[6].mFileTimeLo = 0U;
  c4_info[6].mFileTimeHi = 0U;
  c4_info[7].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  c4_info[7].name = "intmax";
  c4_info[7].dominantType = "char";
  c4_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  c4_info[7].fileTimeLo = 1311255316U;
  c4_info[7].fileTimeHi = 0U;
  c4_info[7].mFileTimeLo = 0U;
  c4_info[7].mFileTimeHi = 0U;
  c4_info[8].context = "";
  c4_info[8].name = "mtimes";
  c4_info[8].dominantType = "double";
  c4_info[8].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c4_info[8].fileTimeLo = 1289519692U;
  c4_info[8].fileTimeHi = 0U;
  c4_info[8].mFileTimeLo = 0U;
  c4_info[8].mFileTimeHi = 0U;
  c4_info[9].context = "";
  c4_info[9].name = "mrdivide";
  c4_info[9].dominantType = "double";
  c4_info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c4_info[9].fileTimeLo = 1357951548U;
  c4_info[9].fileTimeHi = 0U;
  c4_info[9].mFileTimeLo = 1319729966U;
  c4_info[9].mFileTimeHi = 0U;
  c4_info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c4_info[10].name = "rdivide";
  c4_info[10].dominantType = "double";
  c4_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c4_info[10].fileTimeLo = 1346510388U;
  c4_info[10].fileTimeHi = 0U;
  c4_info[10].mFileTimeLo = 0U;
  c4_info[10].mFileTimeHi = 0U;
  c4_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c4_info[11].name = "eml_scalexp_compatible";
  c4_info[11].dominantType = "double";
  c4_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  c4_info[11].fileTimeLo = 1286818796U;
  c4_info[11].fileTimeHi = 0U;
  c4_info[11].mFileTimeLo = 0U;
  c4_info[11].mFileTimeHi = 0U;
  c4_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c4_info[12].name = "eml_div";
  c4_info[12].dominantType = "double";
  c4_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c4_info[12].fileTimeLo = 1313347810U;
  c4_info[12].fileTimeHi = 0U;
  c4_info[12].mFileTimeLo = 0U;
  c4_info[12].mFileTimeHi = 0U;
  c4_info[13].context = "";
  c4_info[13].name = "sin";
  c4_info[13].dominantType = "double";
  c4_info[13].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c4_info[13].fileTimeLo = 1343830386U;
  c4_info[13].fileTimeHi = 0U;
  c4_info[13].mFileTimeLo = 0U;
  c4_info[13].mFileTimeHi = 0U;
  c4_info[14].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c4_info[14].name = "eml_scalar_sin";
  c4_info[14].dominantType = "double";
  c4_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  c4_info[14].fileTimeLo = 1286818736U;
  c4_info[14].fileTimeHi = 0U;
  c4_info[14].mFileTimeLo = 0U;
  c4_info[14].mFileTimeHi = 0U;
  c4_info[15].context = "";
  c4_info[15].name = "cos";
  c4_info[15].dominantType = "double";
  c4_info[15].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c4_info[15].fileTimeLo = 1343830372U;
  c4_info[15].fileTimeHi = 0U;
  c4_info[15].mFileTimeLo = 0U;
  c4_info[15].mFileTimeHi = 0U;
  c4_info[16].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c4_info[16].name = "eml_scalar_cos";
  c4_info[16].dominantType = "double";
  c4_info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  c4_info[16].fileTimeLo = 1286818722U;
  c4_info[16].fileTimeHi = 0U;
  c4_info[16].mFileTimeLo = 0U;
  c4_info[16].mFileTimeHi = 0U;
  c4_info[17].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c4_info[17].name = "eml_index_class";
  c4_info[17].dominantType = "";
  c4_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c4_info[17].fileTimeLo = 1323170578U;
  c4_info[17].fileTimeHi = 0U;
  c4_info[17].mFileTimeLo = 0U;
  c4_info[17].mFileTimeHi = 0U;
  c4_info[18].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c4_info[18].name = "eml_scalar_eg";
  c4_info[18].dominantType = "double";
  c4_info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c4_info[18].fileTimeLo = 1286818796U;
  c4_info[18].fileTimeHi = 0U;
  c4_info[18].mFileTimeLo = 0U;
  c4_info[18].mFileTimeHi = 0U;
  c4_info[19].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c4_info[19].name = "eml_xgemm";
  c4_info[19].dominantType = "char";
  c4_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c4_info[19].fileTimeLo = 1299076772U;
  c4_info[19].fileTimeHi = 0U;
  c4_info[19].mFileTimeLo = 0U;
  c4_info[19].mFileTimeHi = 0U;
  c4_info[20].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c4_info[20].name = "eml_blas_inline";
  c4_info[20].dominantType = "";
  c4_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c4_info[20].fileTimeLo = 1299076768U;
  c4_info[20].fileTimeHi = 0U;
  c4_info[20].mFileTimeLo = 0U;
  c4_info[20].mFileTimeHi = 0U;
  c4_info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold";
  c4_info[21].name = "mtimes";
  c4_info[21].dominantType = "double";
  c4_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c4_info[21].fileTimeLo = 1289519692U;
  c4_info[21].fileTimeHi = 0U;
  c4_info[21].mFileTimeLo = 0U;
  c4_info[21].mFileTimeHi = 0U;
  c4_info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c4_info[22].name = "eml_index_class";
  c4_info[22].dominantType = "";
  c4_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c4_info[22].fileTimeLo = 1323170578U;
  c4_info[22].fileTimeHi = 0U;
  c4_info[22].mFileTimeLo = 0U;
  c4_info[22].mFileTimeHi = 0U;
  c4_info[23].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c4_info[23].name = "eml_scalar_eg";
  c4_info[23].dominantType = "double";
  c4_info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c4_info[23].fileTimeLo = 1286818796U;
  c4_info[23].fileTimeHi = 0U;
  c4_info[23].mFileTimeLo = 0U;
  c4_info[23].mFileTimeHi = 0U;
  c4_info[24].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c4_info[24].name = "eml_refblas_xgemm";
  c4_info[24].dominantType = "char";
  c4_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c4_info[24].fileTimeLo = 1299076774U;
  c4_info[24].fileTimeHi = 0U;
  c4_info[24].mFileTimeLo = 0U;
  c4_info[24].mFileTimeHi = 0U;
}

static real_T c4_sum(SFc4_testInstanceStruct *chartInstance, real_T c4_x[6])
{
  real_T c4_y;
  int32_T c4_k;
  int32_T c4_b_k;
  c4_y = c4_x[0];
  for (c4_k = 2; c4_k < 7; c4_k++) {
    c4_b_k = c4_k;
    c4_y += c4_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c4_b_k), 1, 6, 1, 0) - 1];
  }

  return c4_y;
}

static void c4_eml_scalar_eg(SFc4_testInstanceStruct *chartInstance)
{
}

static const mxArray *c4_f_sf_marshallOut(void *chartInstanceVoid, void
  *c4_inData)
{
  const mxArray *c4_mxArrayOutData = NULL;
  int32_T c4_u;
  const mxArray *c4_y = NULL;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_mxArrayOutData = NULL;
  c4_u = *(int32_T *)c4_inData;
  c4_y = NULL;
  sf_mex_assign(&c4_y, sf_mex_create("y", &c4_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c4_mxArrayOutData, c4_y, FALSE);
  return c4_mxArrayOutData;
}

static int32_T c4_e_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  int32_T c4_y;
  int32_T c4_i59;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_i59, 1, 6, 0U, 0, 0U, 0);
  c4_y = c4_i59;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void c4_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c4_mxArrayInData, const char_T *c4_varName, void *c4_outData)
{
  const mxArray *c4_b_sfEvent;
  const char_T *c4_identifier;
  emlrtMsgIdentifier c4_thisId;
  int32_T c4_y;
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)chartInstanceVoid;
  c4_b_sfEvent = sf_mex_dup(c4_mxArrayInData);
  c4_identifier = c4_varName;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_sfEvent),
    &c4_thisId);
  sf_mex_destroy(&c4_b_sfEvent);
  *(int32_T *)c4_outData = c4_y;
  sf_mex_destroy(&c4_mxArrayInData);
}

static uint8_T c4_f_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_b_is_active_c4_test, const char_T *c4_identifier)
{
  uint8_T c4_y;
  emlrtMsgIdentifier c4_thisId;
  c4_thisId.fIdentifier = c4_identifier;
  c4_thisId.fParent = NULL;
  c4_y = c4_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c4_b_is_active_c4_test),
    &c4_thisId);
  sf_mex_destroy(&c4_b_is_active_c4_test);
  return c4_y;
}

static uint8_T c4_g_emlrt_marshallIn(SFc4_testInstanceStruct *chartInstance,
  const mxArray *c4_u, const emlrtMsgIdentifier *c4_parentId)
{
  uint8_T c4_y;
  uint8_T c4_u0;
  sf_mex_import(c4_parentId, sf_mex_dup(c4_u), &c4_u0, 1, 3, 0U, 0, 0U, 0);
  c4_y = c4_u0;
  sf_mex_destroy(&c4_u);
  return c4_y;
}

static void init_dsm_address_info(SFc4_testInstanceStruct *chartInstance)
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

void sf_c4_test_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3642607870U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1541752355U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3189045559U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1072265335U);
}

mxArray *sf_c4_test_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("fGQowVdhPwVkXjeY2GFgRH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,8,3,dataFields);

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
      pr[0] = (double)(3);
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
      pr[0] = (double)(2);
      pr[1] = (double)(2);
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
      pr[0] = (double)(2);
      pr[1] = (double)(2);
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
      pr[0] = (double)(6);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c4_test_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c4_test(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[5],T\"phi_ref\",},{M[1],M[7],T\"theta_ref\",},{M[8],M[0],T\"is_active_c4_test\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c4_test_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc4_testInstanceStruct *chartInstance;
    chartInstance = (SFc4_testInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _testMachineNumber_,
           4,
           1,
           1,
           10,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"eta");
          _SFD_SET_DATA_PROPS(1,1,1,0,"eta_dot");
          _SFD_SET_DATA_PROPS(2,2,0,1,"phi_ref");
          _SFD_SET_DATA_PROPS(3,1,1,0,"pos_ref");
          _SFD_SET_DATA_PROPS(4,2,0,1,"theta_ref");
          _SFD_SET_DATA_PROPS(5,1,1,0,"m");
          _SFD_SET_DATA_PROPS(6,1,1,0,"K_p_pos");
          _SFD_SET_DATA_PROPS(7,1,1,0,"K_d_pos");
          _SFD_SET_DATA_PROPS(8,1,1,0,"omega");
          _SFD_SET_DATA_PROPS(9,1,1,0,"k");
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
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,457);
        _SFD_CV_INIT_EML_IF(0,1,0,133,140,-1,157);
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
            1.0,0,0,(MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)c4_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)c4_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[2];
          dimVector[0]= 2;
          dimVector[1]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 2;
          dimVector[1]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c4_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c4_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c4_phi_ref;
          real_T *c4_theta_ref;
          real_T *c4_m;
          real_T *c4_k;
          real_T (*c4_eta)[6];
          real_T (*c4_eta_dot)[6];
          real_T (*c4_pos_ref)[3];
          real_T (*c4_K_p_pos)[4];
          real_T (*c4_K_d_pos)[4];
          real_T (*c4_omega)[6];
          c4_k = (real_T *)ssGetInputPortSignal(chartInstance->S, 7);
          c4_omega = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 6);
          c4_K_d_pos = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 5);
          c4_K_p_pos = (real_T (*)[4])ssGetInputPortSignal(chartInstance->S, 4);
          c4_m = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c4_theta_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c4_pos_ref = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 2);
          c4_phi_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          c4_eta_dot = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
          c4_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c4_eta);
          _SFD_SET_DATA_VALUE_PTR(1U, *c4_eta_dot);
          _SFD_SET_DATA_VALUE_PTR(2U, c4_phi_ref);
          _SFD_SET_DATA_VALUE_PTR(3U, *c4_pos_ref);
          _SFD_SET_DATA_VALUE_PTR(4U, c4_theta_ref);
          _SFD_SET_DATA_VALUE_PTR(5U, c4_m);
          _SFD_SET_DATA_VALUE_PTR(6U, *c4_K_p_pos);
          _SFD_SET_DATA_VALUE_PTR(7U, *c4_K_d_pos);
          _SFD_SET_DATA_VALUE_PTR(8U, *c4_omega);
          _SFD_SET_DATA_VALUE_PTR(9U, c4_k);
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
  return "vGt0XTKU6Zzq0E0GCv18QH";
}

static void sf_opaque_initialize_c4_test(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc4_testInstanceStruct*) chartInstanceVar)->S,0);
  initialize_params_c4_test((SFc4_testInstanceStruct*) chartInstanceVar);
  initialize_c4_test((SFc4_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c4_test(void *chartInstanceVar)
{
  enable_c4_test((SFc4_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c4_test(void *chartInstanceVar)
{
  disable_c4_test((SFc4_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c4_test(void *chartInstanceVar)
{
  sf_c4_test((SFc4_testInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c4_test(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c4_test((SFc4_testInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c4_test();/* state var info */
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

extern void sf_internal_set_sim_state_c4_test(SimStruct* S, const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c4_test();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c4_test((SFc4_testInstanceStruct*)chartInfo->chartInstance,
                        mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c4_test(SimStruct* S)
{
  return sf_internal_get_sim_state_c4_test(S);
}

static void sf_opaque_set_sim_state_c4_test(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c4_test(S, st);
}

static void sf_opaque_terminate_c4_test(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc4_testInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_optimization_info();
    }

    finalize_c4_test((SFc4_testInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc4_test((SFc4_testInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c4_test(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c4_test((SFc4_testInstanceStruct*)(((ChartInfoStruct *)
      ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c4_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      4);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,4,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,4,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,4);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 7, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,4,8);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,4,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 8; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,4);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3793238913U));
  ssSetChecksum1(S,(743881160U));
  ssSetChecksum2(S,(3865548663U));
  ssSetChecksum3(S,(1467882396U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c4_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c4_test(SimStruct *S)
{
  SFc4_testInstanceStruct *chartInstance;
  chartInstance = (SFc4_testInstanceStruct *)utMalloc(sizeof
    (SFc4_testInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc4_testInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c4_test;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c4_test;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c4_test;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c4_test;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c4_test;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c4_test;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c4_test;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c4_test;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c4_test;
  chartInstance->chartInfo.mdlStart = mdlStart_c4_test;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c4_test;
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

void c4_test_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c4_test(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c4_test(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c4_test(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c4_test_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
