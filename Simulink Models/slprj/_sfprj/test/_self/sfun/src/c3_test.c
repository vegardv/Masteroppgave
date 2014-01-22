/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_sfun.h"
#include "c3_test.h"
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
static const char * c3_debug_family_names[8] = { "J", "J11", "J22", "nargin",
  "nargout", "eta", "nu", "eta_dot" };

static const char * c3_b_debug_family_names[12] = { "cphi", "sphi", "cth", "sth",
  "cpsi", "spsi", "nargin", "nargout", "phi", "theta", "psi", "R" };

static const char * c3_c_debug_family_names[12] = { "cphi", "sphi", "cth", "sth",
  "nargin", "nargout", "phi", "theta", "psi", "J", "J1", "J2" };

/* Function Declarations */
static void initialize_c3_test(SFc3_testInstanceStruct *chartInstance);
static void initialize_params_c3_test(SFc3_testInstanceStruct *chartInstance);
static void enable_c3_test(SFc3_testInstanceStruct *chartInstance);
static void disable_c3_test(SFc3_testInstanceStruct *chartInstance);
static void c3_update_debugger_state_c3_test(SFc3_testInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c3_test(SFc3_testInstanceStruct
  *chartInstance);
static void set_sim_state_c3_test(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_st);
static void finalize_c3_test(SFc3_testInstanceStruct *chartInstance);
static void sf_c3_test(SFc3_testInstanceStruct *chartInstance);
static void c3_chartstep_c3_test(SFc3_testInstanceStruct *chartInstance);
static void initSimStructsc3_test(SFc3_testInstanceStruct *chartInstance);
static void registerMessagesc3_test(SFc3_testInstanceStruct *chartInstance);
static void c3_eulerang(SFc3_testInstanceStruct *chartInstance, real_T c3_phi,
  real_T c3_theta, real_T c3_psi, real_T c3_J[36], real_T c3_J1[9], real_T
  c3_J2[9]);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static void c3_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_eta_dot, const char_T *c3_identifier, real_T c3_y[6]);
static void c3_b_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[6]);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static real_T c3_c_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_d_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[9]);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static void c3_e_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[36]);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(c3_ResolvedFunctionInfo c3_info[24]);
static void c3_eml_scalar_eg(SFc3_testInstanceStruct *chartInstance);
static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_f_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static uint8_T c3_g_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_b_is_active_c3_test, const char_T *c3_identifier);
static uint8_T c3_h_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void init_dsm_address_info(SFc3_testInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c3_test(SFc3_testInstanceStruct *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c3_is_active_c3_test = 0U;
}

static void initialize_params_c3_test(SFc3_testInstanceStruct *chartInstance)
{
}

static void enable_c3_test(SFc3_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c3_test(SFc3_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c3_update_debugger_state_c3_test(SFc3_testInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c3_test(SFc3_testInstanceStruct
  *chartInstance)
{
  const mxArray *c3_st;
  const mxArray *c3_y = NULL;
  int32_T c3_i0;
  real_T c3_u[6];
  const mxArray *c3_b_y = NULL;
  uint8_T c3_hoistedGlobal;
  uint8_T c3_b_u;
  const mxArray *c3_c_y = NULL;
  real_T (*c3_eta_dot)[6];
  c3_eta_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_st = NULL;
  c3_st = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellarray(2), FALSE);
  for (c3_i0 = 0; c3_i0 < 6; c3_i0++) {
    c3_u[c3_i0] = (*c3_eta_dot)[c3_i0];
  }

  c3_b_y = NULL;
  sf_mex_assign(&c3_b_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_setcell(c3_y, 0, c3_b_y);
  c3_hoistedGlobal = chartInstance->c3_is_active_c3_test;
  c3_b_u = c3_hoistedGlobal;
  c3_c_y = NULL;
  sf_mex_assign(&c3_c_y, sf_mex_create("y", &c3_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c3_y, 1, c3_c_y);
  sf_mex_assign(&c3_st, c3_y, FALSE);
  return c3_st;
}

static void set_sim_state_c3_test(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_st)
{
  const mxArray *c3_u;
  real_T c3_dv0[6];
  int32_T c3_i1;
  real_T (*c3_eta_dot)[6];
  c3_eta_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c3_doneDoubleBufferReInit = TRUE;
  c3_u = sf_mex_dup(c3_st);
  c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)),
                      "eta_dot", c3_dv0);
  for (c3_i1 = 0; c3_i1 < 6; c3_i1++) {
    (*c3_eta_dot)[c3_i1] = c3_dv0[c3_i1];
  }

  chartInstance->c3_is_active_c3_test = c3_g_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c3_u, 1)), "is_active_c3_test");
  sf_mex_destroy(&c3_u);
  c3_update_debugger_state_c3_test(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void finalize_c3_test(SFc3_testInstanceStruct *chartInstance)
{
}

static void sf_c3_test(SFc3_testInstanceStruct *chartInstance)
{
  int32_T c3_i2;
  int32_T c3_i3;
  int32_T c3_i4;
  real_T (*c3_nu)[6];
  real_T (*c3_eta_dot)[6];
  real_T (*c3_eta)[6];
  c3_nu = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c3_eta_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  for (c3_i2 = 0; c3_i2 < 6; c3_i2++) {
    _SFD_DATA_RANGE_CHECK((*c3_eta)[c3_i2], 0U);
  }

  for (c3_i3 = 0; c3_i3 < 6; c3_i3++) {
    _SFD_DATA_RANGE_CHECK((*c3_eta_dot)[c3_i3], 1U);
  }

  for (c3_i4 = 0; c3_i4 < 6; c3_i4++) {
    _SFD_DATA_RANGE_CHECK((*c3_nu)[c3_i4], 2U);
  }

  chartInstance->c3_sfEvent = CALL_EVENT;
  c3_chartstep_c3_test(chartInstance);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_testMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c3_chartstep_c3_test(SFc3_testInstanceStruct *chartInstance)
{
  int32_T c3_i5;
  real_T c3_eta[6];
  int32_T c3_i6;
  real_T c3_nu[6];
  uint32_T c3_debug_family_var_map[8];
  real_T c3_J[36];
  real_T c3_J11[9];
  real_T c3_J22[9];
  real_T c3_nargin = 2.0;
  real_T c3_nargout = 1.0;
  real_T c3_eta_dot[6];
  real_T c3_b_J22[9];
  real_T c3_b_J11[9];
  real_T c3_b_J[36];
  int32_T c3_i7;
  int32_T c3_i8;
  int32_T c3_i9;
  int32_T c3_i10;
  int32_T c3_i11;
  real_T c3_b[6];
  int32_T c3_i12;
  int32_T c3_i13;
  int32_T c3_i14;
  real_T c3_C[6];
  int32_T c3_i15;
  int32_T c3_i16;
  int32_T c3_i17;
  int32_T c3_i18;
  int32_T c3_i19;
  int32_T c3_i20;
  int32_T c3_i21;
  real_T (*c3_b_eta_dot)[6];
  real_T (*c3_b_nu)[6];
  real_T (*c3_b_eta)[6];
  c3_b_nu = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c3_b_eta_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
  c3_b_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
  for (c3_i5 = 0; c3_i5 < 6; c3_i5++) {
    c3_eta[c3_i5] = (*c3_b_eta)[c3_i5];
  }

  for (c3_i6 = 0; c3_i6 < 6; c3_i6++) {
    c3_nu[c3_i6] = (*c3_b_nu)[c3_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 8U, 8U, c3_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_J, 0U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_J11, 1U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_J22, 2U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 3U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 4U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_eta, 5U, c3_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c3_nu, 6U, c3_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_eta_dot, 7U, c3_sf_marshallOut,
    c3_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 4);
  c3_eulerang(chartInstance, c3_eta[3], c3_eta[4], c3_eta[5], c3_b_J, c3_b_J11,
              c3_b_J22);
  for (c3_i7 = 0; c3_i7 < 36; c3_i7++) {
    c3_J[c3_i7] = c3_b_J[c3_i7];
  }

  for (c3_i8 = 0; c3_i8 < 9; c3_i8++) {
    c3_J11[c3_i8] = c3_b_J11[c3_i8];
  }

  for (c3_i9 = 0; c3_i9 < 9; c3_i9++) {
    c3_J22[c3_i9] = c3_b_J22[c3_i9];
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, 6);
  for (c3_i10 = 0; c3_i10 < 36; c3_i10++) {
    c3_b_J[c3_i10] = c3_J[c3_i10];
  }

  for (c3_i11 = 0; c3_i11 < 6; c3_i11++) {
    c3_b[c3_i11] = c3_nu[c3_i11];
  }

  c3_eml_scalar_eg(chartInstance);
  c3_eml_scalar_eg(chartInstance);
  for (c3_i12 = 0; c3_i12 < 6; c3_i12++) {
    c3_eta_dot[c3_i12] = 0.0;
  }

  for (c3_i13 = 0; c3_i13 < 6; c3_i13++) {
    c3_eta_dot[c3_i13] = 0.0;
  }

  for (c3_i14 = 0; c3_i14 < 6; c3_i14++) {
    c3_C[c3_i14] = c3_eta_dot[c3_i14];
  }

  for (c3_i15 = 0; c3_i15 < 6; c3_i15++) {
    c3_eta_dot[c3_i15] = c3_C[c3_i15];
  }

  for (c3_i16 = 0; c3_i16 < 6; c3_i16++) {
    c3_C[c3_i16] = c3_eta_dot[c3_i16];
  }

  for (c3_i17 = 0; c3_i17 < 6; c3_i17++) {
    c3_eta_dot[c3_i17] = c3_C[c3_i17];
  }

  for (c3_i18 = 0; c3_i18 < 6; c3_i18++) {
    c3_eta_dot[c3_i18] = 0.0;
    c3_i19 = 0;
    for (c3_i20 = 0; c3_i20 < 6; c3_i20++) {
      c3_eta_dot[c3_i18] += c3_b_J[c3_i19 + c3_i18] * c3_b[c3_i20];
      c3_i19 += 6;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c3_sfEvent, -6);
  _SFD_SYMBOL_SCOPE_POP();
  for (c3_i21 = 0; c3_i21 < 6; c3_i21++) {
    (*c3_b_eta_dot)[c3_i21] = c3_eta_dot[c3_i21];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
}

static void initSimStructsc3_test(SFc3_testInstanceStruct *chartInstance)
{
}

static void registerMessagesc3_test(SFc3_testInstanceStruct *chartInstance)
{
}

static void c3_eulerang(SFc3_testInstanceStruct *chartInstance, real_T c3_phi,
  real_T c3_theta, real_T c3_psi, real_T c3_J[36], real_T c3_J1[9], real_T
  c3_J2[9])
{
  uint32_T c3_debug_family_var_map[12];
  real_T c3_cphi;
  real_T c3_sphi;
  real_T c3_cth;
  real_T c3_sth;
  real_T c3_nargin = 3.0;
  real_T c3_nargout = 3.0;
  real_T c3_x;
  real_T c3_b_x;
  real_T c3_c_x;
  real_T c3_d_x;
  real_T c3_e_x;
  real_T c3_f_x;
  real_T c3_g_x;
  real_T c3_h_x;
  real_T c3_b_phi;
  real_T c3_b_theta;
  real_T c3_b_psi;
  real_T c3_b_cphi;
  real_T c3_b_sphi;
  real_T c3_b_cth;
  real_T c3_b_sth;
  real_T c3_cpsi;
  real_T c3_spsi;
  real_T c3_b_nargin = 3.0;
  real_T c3_b_nargout = 1.0;
  real_T c3_i_x;
  real_T c3_j_x;
  real_T c3_k_x;
  real_T c3_l_x;
  real_T c3_m_x;
  real_T c3_n_x;
  real_T c3_o_x;
  real_T c3_p_x;
  real_T c3_q_x;
  real_T c3_r_x;
  real_T c3_s_x;
  real_T c3_t_x;
  real_T c3_a;
  real_T c3_b;
  real_T c3_y;
  real_T c3_b_a;
  real_T c3_b_b;
  real_T c3_b_y;
  real_T c3_c_a;
  real_T c3_c_b;
  real_T c3_c_y;
  real_T c3_d_a;
  real_T c3_d_b;
  real_T c3_d_y;
  real_T c3_e_a;
  real_T c3_e_b;
  real_T c3_e_y;
  real_T c3_f_a;
  real_T c3_f_b;
  real_T c3_f_y;
  real_T c3_g_a;
  real_T c3_g_b;
  real_T c3_g_y;
  real_T c3_h_a;
  real_T c3_h_b;
  real_T c3_h_y;
  real_T c3_i_a;
  real_T c3_i_b;
  real_T c3_i_y;
  real_T c3_j_a;
  real_T c3_j_b;
  real_T c3_j_y;
  real_T c3_k_a;
  real_T c3_k_b;
  real_T c3_k_y;
  real_T c3_l_a;
  real_T c3_l_b;
  real_T c3_l_y;
  real_T c3_m_a;
  real_T c3_m_b;
  real_T c3_m_y;
  real_T c3_n_a;
  real_T c3_n_b;
  real_T c3_n_y;
  real_T c3_o_a;
  real_T c3_o_b;
  real_T c3_o_y;
  real_T c3_p_a;
  real_T c3_p_b;
  real_T c3_p_y;
  int32_T c3_i22;
  static char_T c3_varargin_1[39] = { 'J', '2', ' ', 'i', 's', ' ', 's', 'i',
    'n', 'g', 'u', 'l', 'a', 'r', ' ', 'f', 'o', 'r', ' ', 't', 'h', 'e', 't',
    'a', ' ', '=', ' ', '+', '-', '9', '0', ' ', 'd', 'e', 'g', 'r', 'e', 'e',
    's' };

  char_T c3_u[39];
  const mxArray *c3_q_y = NULL;
  real_T c3_q_a;
  real_T c3_q_b;
  real_T c3_r_y;
  real_T c3_A;
  real_T c3_B;
  real_T c3_u_x;
  real_T c3_s_y;
  real_T c3_v_x;
  real_T c3_t_y;
  real_T c3_u_y;
  real_T c3_r_a;
  real_T c3_r_b;
  real_T c3_v_y;
  real_T c3_b_A;
  real_T c3_b_B;
  real_T c3_w_x;
  real_T c3_w_y;
  real_T c3_x_x;
  real_T c3_x_y;
  real_T c3_y_y;
  real_T c3_c_A;
  real_T c3_c_B;
  real_T c3_y_x;
  real_T c3_ab_y;
  real_T c3_ab_x;
  real_T c3_bb_y;
  real_T c3_cb_y;
  real_T c3_d_A;
  real_T c3_d_B;
  real_T c3_bb_x;
  real_T c3_db_y;
  real_T c3_cb_x;
  real_T c3_eb_y;
  real_T c3_fb_y;
  int32_T c3_i23;
  int32_T c3_i24;
  int32_T c3_i25;
  int32_T c3_i26;
  int32_T c3_i27;
  int32_T c3_i28;
  int32_T c3_i29;
  int32_T c3_i30;
  int32_T c3_i31;
  int32_T c3_i32;
  int32_T c3_i33;
  int32_T c3_i34;
  int32_T c3_i35;
  int32_T c3_i36;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 12U, 12U, c3_c_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cphi, 0U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_sphi, 1U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cth, 2U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_sth, 3U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 4U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 5U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_phi, 6U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_theta, 7U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_psi, 8U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_J, 9U, c3_d_sf_marshallOut,
    c3_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_J1, 10U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_J2, 11U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  CV_SCRIPT_FCN(0, 0);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 32);
  c3_x = c3_phi;
  c3_cphi = c3_x;
  c3_b_x = c3_cphi;
  c3_cphi = c3_b_x;
  c3_cphi = muDoubleScalarCos(c3_cphi);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 33);
  c3_c_x = c3_phi;
  c3_sphi = c3_c_x;
  c3_d_x = c3_sphi;
  c3_sphi = c3_d_x;
  c3_sphi = muDoubleScalarSin(c3_sphi);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 34);
  c3_e_x = c3_theta;
  c3_cth = c3_e_x;
  c3_f_x = c3_cth;
  c3_cth = c3_f_x;
  c3_cth = muDoubleScalarCos(c3_cth);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 35);
  c3_g_x = c3_theta;
  c3_sth = c3_g_x;
  c3_h_x = c3_sth;
  c3_sth = c3_h_x;
  c3_sth = muDoubleScalarSin(c3_sth);
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 37);
  c3_b_phi = c3_phi;
  c3_b_theta = c3_theta;
  c3_b_psi = c3_psi;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 12U, 12U, c3_b_debug_family_names,
    c3_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_cphi, 0U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_sphi, 1U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_cth, 2U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_sth, 3U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cpsi, 4U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_spsi, 5U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargin, 6U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargout, 7U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_phi, 8U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_theta, 9U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_psi, 10U, c3_b_sf_marshallOut,
    c3_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c3_J1, 11U, c3_c_sf_marshallOut,
    c3_c_sf_marshallIn);
  CV_SCRIPT_FCN(1, 0);
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, 31);
  c3_i_x = c3_b_phi;
  c3_b_cphi = c3_i_x;
  c3_j_x = c3_b_cphi;
  c3_b_cphi = c3_j_x;
  c3_b_cphi = muDoubleScalarCos(c3_b_cphi);
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, 32);
  c3_k_x = c3_b_phi;
  c3_b_sphi = c3_k_x;
  c3_l_x = c3_b_sphi;
  c3_b_sphi = c3_l_x;
  c3_b_sphi = muDoubleScalarSin(c3_b_sphi);
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, 33);
  c3_m_x = c3_b_theta;
  c3_b_cth = c3_m_x;
  c3_n_x = c3_b_cth;
  c3_b_cth = c3_n_x;
  c3_b_cth = muDoubleScalarCos(c3_b_cth);
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, 34);
  c3_o_x = c3_b_theta;
  c3_b_sth = c3_o_x;
  c3_p_x = c3_b_sth;
  c3_b_sth = c3_p_x;
  c3_b_sth = muDoubleScalarSin(c3_b_sth);
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, 35);
  c3_q_x = c3_b_psi;
  c3_cpsi = c3_q_x;
  c3_r_x = c3_cpsi;
  c3_cpsi = c3_r_x;
  c3_cpsi = muDoubleScalarCos(c3_cpsi);
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, 36);
  c3_s_x = c3_b_psi;
  c3_spsi = c3_s_x;
  c3_t_x = c3_spsi;
  c3_spsi = c3_t_x;
  c3_spsi = muDoubleScalarSin(c3_spsi);
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, 38);
  c3_a = c3_cpsi;
  c3_b = c3_b_cth;
  c3_y = c3_a * c3_b;
  c3_b_a = -c3_spsi;
  c3_b_b = c3_b_cphi;
  c3_b_y = c3_b_a * c3_b_b;
  c3_c_a = c3_cpsi;
  c3_c_b = c3_b_sth;
  c3_c_y = c3_c_a * c3_c_b;
  c3_d_a = c3_c_y;
  c3_d_b = c3_b_sphi;
  c3_d_y = c3_d_a * c3_d_b;
  c3_e_a = c3_spsi;
  c3_e_b = c3_b_sphi;
  c3_e_y = c3_e_a * c3_e_b;
  c3_f_a = c3_cpsi;
  c3_f_b = c3_b_cphi;
  c3_f_y = c3_f_a * c3_f_b;
  c3_g_a = c3_f_y;
  c3_g_b = c3_b_sth;
  c3_g_y = c3_g_a * c3_g_b;
  c3_h_a = c3_spsi;
  c3_h_b = c3_b_cth;
  c3_h_y = c3_h_a * c3_h_b;
  c3_i_a = c3_cpsi;
  c3_i_b = c3_b_cphi;
  c3_i_y = c3_i_a * c3_i_b;
  c3_j_a = c3_b_sphi;
  c3_j_b = c3_b_sth;
  c3_j_y = c3_j_a * c3_j_b;
  c3_k_a = c3_j_y;
  c3_k_b = c3_spsi;
  c3_k_y = c3_k_a * c3_k_b;
  c3_l_a = -c3_cpsi;
  c3_l_b = c3_b_sphi;
  c3_l_y = c3_l_a * c3_l_b;
  c3_m_a = c3_b_sth;
  c3_m_b = c3_spsi;
  c3_m_y = c3_m_a * c3_m_b;
  c3_n_a = c3_m_y;
  c3_n_b = c3_b_cphi;
  c3_n_y = c3_n_a * c3_n_b;
  c3_o_a = c3_b_cth;
  c3_o_b = c3_b_sphi;
  c3_o_y = c3_o_a * c3_o_b;
  c3_p_a = c3_b_cth;
  c3_p_b = c3_b_cphi;
  c3_p_y = c3_p_a * c3_p_b;
  c3_J1[0] = c3_y;
  c3_J1[3] = c3_b_y + c3_d_y;
  c3_J1[6] = c3_e_y + c3_g_y;
  c3_J1[1] = c3_h_y;
  c3_J1[4] = c3_i_y + c3_k_y;
  c3_J1[7] = c3_l_y + c3_n_y;
  c3_J1[2] = -c3_b_sth;
  c3_J1[5] = c3_o_y;
  c3_J1[8] = c3_p_y;
  _SFD_SCRIPT_CALL(1U, chartInstance->c3_sfEvent, -38);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 39);
  if (CV_SCRIPT_IF(0, 0, c3_cth == 0.0)) {
    _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 39);
    for (c3_i22 = 0; c3_i22 < 39; c3_i22++) {
      c3_u[c3_i22] = c3_varargin_1[c3_i22];
    }

    c3_q_y = NULL;
    sf_mex_assign(&c3_q_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 39),
                  FALSE);
    sf_mex_call_debug("error", 0U, 1U, 14, c3_q_y);
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 41);
  c3_q_a = c3_sphi;
  c3_q_b = c3_sth;
  c3_r_y = c3_q_a * c3_q_b;
  c3_A = c3_r_y;
  c3_B = c3_cth;
  c3_u_x = c3_A;
  c3_s_y = c3_B;
  c3_v_x = c3_u_x;
  c3_t_y = c3_s_y;
  c3_u_y = c3_v_x / c3_t_y;
  c3_r_a = c3_cphi;
  c3_r_b = c3_sth;
  c3_v_y = c3_r_a * c3_r_b;
  c3_b_A = c3_v_y;
  c3_b_B = c3_cth;
  c3_w_x = c3_b_A;
  c3_w_y = c3_b_B;
  c3_x_x = c3_w_x;
  c3_x_y = c3_w_y;
  c3_y_y = c3_x_x / c3_x_y;
  c3_c_A = c3_sphi;
  c3_c_B = c3_cth;
  c3_y_x = c3_c_A;
  c3_ab_y = c3_c_B;
  c3_ab_x = c3_y_x;
  c3_bb_y = c3_ab_y;
  c3_cb_y = c3_ab_x / c3_bb_y;
  c3_d_A = c3_cphi;
  c3_d_B = c3_cth;
  c3_bb_x = c3_d_A;
  c3_db_y = c3_d_B;
  c3_cb_x = c3_bb_x;
  c3_eb_y = c3_db_y;
  c3_fb_y = c3_cb_x / c3_eb_y;
  c3_J2[0] = 1.0;
  c3_J2[3] = c3_u_y;
  c3_J2[6] = c3_y_y;
  c3_J2[1] = 0.0;
  c3_J2[4] = c3_cphi;
  c3_J2[7] = -c3_sphi;
  c3_J2[2] = 0.0;
  c3_J2[5] = c3_cb_y;
  c3_J2[8] = c3_fb_y;
  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, 46);
  c3_i23 = 0;
  c3_i24 = 0;
  for (c3_i25 = 0; c3_i25 < 3; c3_i25++) {
    for (c3_i26 = 0; c3_i26 < 3; c3_i26++) {
      c3_J[c3_i26 + c3_i23] = c3_J1[c3_i26 + c3_i24];
    }

    c3_i23 += 6;
    c3_i24 += 3;
  }

  c3_i27 = 0;
  for (c3_i28 = 0; c3_i28 < 3; c3_i28++) {
    for (c3_i29 = 0; c3_i29 < 3; c3_i29++) {
      c3_J[(c3_i29 + c3_i27) + 18] = 0.0;
    }

    c3_i27 += 6;
  }

  c3_i30 = 0;
  for (c3_i31 = 0; c3_i31 < 3; c3_i31++) {
    for (c3_i32 = 0; c3_i32 < 3; c3_i32++) {
      c3_J[(c3_i32 + c3_i30) + 3] = 0.0;
    }

    c3_i30 += 6;
  }

  c3_i33 = 0;
  c3_i34 = 0;
  for (c3_i35 = 0; c3_i35 < 3; c3_i35++) {
    for (c3_i36 = 0; c3_i36 < 3; c3_i36++) {
      c3_J[(c3_i36 + c3_i33) + 21] = c3_J2[c3_i36 + c3_i34];
    }

    c3_i33 += 6;
    c3_i34 += 3;
  }

  _SFD_SCRIPT_CALL(0U, chartInstance->c3_sfEvent, -46);
  _SFD_SYMBOL_SCOPE_POP();
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber)
{
  _SFD_SCRIPT_TRANSLATION(c3_chartNumber, 0U, sf_debug_get_script_id(
    "C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m"));
  _SFD_SCRIPT_TRANSLATION(c3_chartNumber, 1U, sf_debug_get_script_id(
    "C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/Rzyx.m"));
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i37;
  real_T c3_b_inData[6];
  int32_T c3_i38;
  real_T c3_u[6];
  const mxArray *c3_y = NULL;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  for (c3_i37 = 0; c3_i37 < 6; c3_i37++) {
    c3_b_inData[c3_i37] = (*(real_T (*)[6])c3_inData)[c3_i37];
  }

  for (c3_i38 = 0; c3_i38 < 6; c3_i38++) {
    c3_u[c3_i38] = c3_b_inData[c3_i38];
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static void c3_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_eta_dot, const char_T *c3_identifier, real_T c3_y[6])
{
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_eta_dot), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_eta_dot);
}

static void c3_b_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[6])
{
  real_T c3_dv1[6];
  int32_T c3_i39;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv1, 1, 0, 0U, 1, 0U, 1, 6);
  for (c3_i39 = 0; c3_i39 < 6; c3_i39++) {
    c3_y[c3_i39] = c3_dv1[c3_i39];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_eta_dot;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[6];
  int32_T c3_i40;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_eta_dot = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_eta_dot), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_eta_dot);
  for (c3_i40 = 0; c3_i40 < 6; c3_i40++) {
    (*(real_T (*)[6])c3_outData)[c3_i40] = c3_y[c3_i40];
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  real_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(real_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static real_T c3_c_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d0;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d0, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_nargout;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_nargout = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  *(real_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i41;
  int32_T c3_i42;
  int32_T c3_i43;
  real_T c3_b_inData[9];
  int32_T c3_i44;
  int32_T c3_i45;
  int32_T c3_i46;
  real_T c3_u[9];
  const mxArray *c3_y = NULL;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i41 = 0;
  for (c3_i42 = 0; c3_i42 < 3; c3_i42++) {
    for (c3_i43 = 0; c3_i43 < 3; c3_i43++) {
      c3_b_inData[c3_i43 + c3_i41] = (*(real_T (*)[9])c3_inData)[c3_i43 + c3_i41];
    }

    c3_i41 += 3;
  }

  c3_i44 = 0;
  for (c3_i45 = 0; c3_i45 < 3; c3_i45++) {
    for (c3_i46 = 0; c3_i46 < 3; c3_i46++) {
      c3_u[c3_i46 + c3_i44] = c3_b_inData[c3_i46 + c3_i44];
    }

    c3_i44 += 3;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static void c3_d_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[9])
{
  real_T c3_dv2[9];
  int32_T c3_i47;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv2, 1, 0, 0U, 1, 0U, 2, 3, 3);
  for (c3_i47 = 0; c3_i47 < 9; c3_i47++) {
    c3_y[c3_i47] = c3_dv2[c3_i47];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_J22;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[9];
  int32_T c3_i48;
  int32_T c3_i49;
  int32_T c3_i50;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_J22 = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_J22), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_J22);
  c3_i48 = 0;
  for (c3_i49 = 0; c3_i49 < 3; c3_i49++) {
    for (c3_i50 = 0; c3_i50 < 3; c3_i50++) {
      (*(real_T (*)[9])c3_outData)[c3_i50 + c3_i48] = c3_y[c3_i50 + c3_i48];
    }

    c3_i48 += 3;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_i51;
  int32_T c3_i52;
  int32_T c3_i53;
  real_T c3_b_inData[36];
  int32_T c3_i54;
  int32_T c3_i55;
  int32_T c3_i56;
  real_T c3_u[36];
  const mxArray *c3_y = NULL;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_i51 = 0;
  for (c3_i52 = 0; c3_i52 < 6; c3_i52++) {
    for (c3_i53 = 0; c3_i53 < 6; c3_i53++) {
      c3_b_inData[c3_i53 + c3_i51] = (*(real_T (*)[36])c3_inData)[c3_i53 +
        c3_i51];
    }

    c3_i51 += 6;
  }

  c3_i54 = 0;
  for (c3_i55 = 0; c3_i55 < 6; c3_i55++) {
    for (c3_i56 = 0; c3_i56 < 6; c3_i56++) {
      c3_u[c3_i56 + c3_i54] = c3_b_inData[c3_i56 + c3_i54];
    }

    c3_i54 += 6;
  }

  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 0, 0U, 1U, 0U, 2, 6, 6), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static void c3_e_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance, const
  mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId, real_T c3_y[36])
{
  real_T c3_dv3[36];
  int32_T c3_i57;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), c3_dv3, 1, 0, 0U, 1, 0U, 2, 6, 6);
  for (c3_i57 = 0; c3_i57 < 36; c3_i57++) {
    c3_y[c3_i57] = c3_dv3[c3_i57];
  }

  sf_mex_destroy(&c3_u);
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_J;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  real_T c3_y[36];
  int32_T c3_i58;
  int32_T c3_i59;
  int32_T c3_i60;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_J = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_J), &c3_thisId, c3_y);
  sf_mex_destroy(&c3_J);
  c3_i58 = 0;
  for (c3_i59 = 0; c3_i59 < 6; c3_i59++) {
    for (c3_i60 = 0; c3_i60 < 6; c3_i60++) {
      (*(real_T (*)[36])c3_outData)[c3_i60 + c3_i58] = c3_y[c3_i60 + c3_i58];
    }

    c3_i58 += 6;
  }

  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray *sf_c3_test_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo;
  c3_ResolvedFunctionInfo c3_info[24];
  const mxArray *c3_m0 = NULL;
  int32_T c3_i61;
  c3_ResolvedFunctionInfo *c3_r0;
  c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  c3_info_helper(c3_info);
  sf_mex_assign(&c3_m0, sf_mex_createstruct("nameCaptureInfo", 1, 24), FALSE);
  for (c3_i61 = 0; c3_i61 < 24; c3_i61++) {
    c3_r0 = &c3_info[c3_i61];
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", c3_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c3_r0->context)), "context", "nameCaptureInfo",
                    c3_i61);
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", c3_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c3_r0->name)), "name", "nameCaptureInfo", c3_i61);
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", c3_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c3_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c3_i61);
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", c3_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c3_r0->resolved)), "resolved", "nameCaptureInfo",
                    c3_i61);
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", &c3_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c3_i61);
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", &c3_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c3_i61);
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", &c3_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c3_i61);
    sf_mex_addfield(c3_m0, sf_mex_create("nameCaptureInfo", &c3_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c3_i61);
  }

  sf_mex_assign(&c3_nameCaptureInfo, c3_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(c3_ResolvedFunctionInfo c3_info[24])
{
  c3_info[0].context = "";
  c3_info[0].name = "eulerang";
  c3_info[0].dominantType = "double";
  c3_info[0].resolved =
    "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m";
  c3_info[0].fileTimeLo = 1206483514U;
  c3_info[0].fileTimeHi = 0U;
  c3_info[0].mFileTimeLo = 0U;
  c3_info[0].mFileTimeHi = 0U;
  c3_info[1].context =
    "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m";
  c3_info[1].name = "cos";
  c3_info[1].dominantType = "double";
  c3_info[1].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c3_info[1].fileTimeLo = 1343830372U;
  c3_info[1].fileTimeHi = 0U;
  c3_info[1].mFileTimeLo = 0U;
  c3_info[1].mFileTimeHi = 0U;
  c3_info[2].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c3_info[2].name = "eml_scalar_cos";
  c3_info[2].dominantType = "double";
  c3_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  c3_info[2].fileTimeLo = 1286818722U;
  c3_info[2].fileTimeHi = 0U;
  c3_info[2].mFileTimeLo = 0U;
  c3_info[2].mFileTimeHi = 0U;
  c3_info[3].context =
    "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m";
  c3_info[3].name = "sin";
  c3_info[3].dominantType = "double";
  c3_info[3].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c3_info[3].fileTimeLo = 1343830386U;
  c3_info[3].fileTimeHi = 0U;
  c3_info[3].mFileTimeLo = 0U;
  c3_info[3].mFileTimeHi = 0U;
  c3_info[4].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c3_info[4].name = "eml_scalar_sin";
  c3_info[4].dominantType = "double";
  c3_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  c3_info[4].fileTimeLo = 1286818736U;
  c3_info[4].fileTimeHi = 0U;
  c3_info[4].mFileTimeLo = 0U;
  c3_info[4].mFileTimeHi = 0U;
  c3_info[5].context =
    "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m";
  c3_info[5].name = "Rzyx";
  c3_info[5].dominantType = "double";
  c3_info[5].resolved = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/Rzyx.m";
  c3_info[5].fileTimeLo = 1206483734U;
  c3_info[5].fileTimeHi = 0U;
  c3_info[5].mFileTimeLo = 0U;
  c3_info[5].mFileTimeHi = 0U;
  c3_info[6].context = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/Rzyx.m";
  c3_info[6].name = "cos";
  c3_info[6].dominantType = "double";
  c3_info[6].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c3_info[6].fileTimeLo = 1343830372U;
  c3_info[6].fileTimeHi = 0U;
  c3_info[6].mFileTimeLo = 0U;
  c3_info[6].mFileTimeHi = 0U;
  c3_info[7].context = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/Rzyx.m";
  c3_info[7].name = "sin";
  c3_info[7].dominantType = "double";
  c3_info[7].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c3_info[7].fileTimeLo = 1343830386U;
  c3_info[7].fileTimeHi = 0U;
  c3_info[7].mFileTimeLo = 0U;
  c3_info[7].mFileTimeHi = 0U;
  c3_info[8].context = "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/Rzyx.m";
  c3_info[8].name = "mtimes";
  c3_info[8].dominantType = "double";
  c3_info[8].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c3_info[8].fileTimeLo = 1289519692U;
  c3_info[8].fileTimeHi = 0U;
  c3_info[8].mFileTimeLo = 0U;
  c3_info[8].mFileTimeHi = 0U;
  c3_info[9].context =
    "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m";
  c3_info[9].name = "error";
  c3_info[9].dominantType = "char";
  c3_info[9].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/lang/error.m";
  c3_info[9].fileTimeLo = 1319729966U;
  c3_info[9].fileTimeHi = 0U;
  c3_info[9].mFileTimeLo = 0U;
  c3_info[9].mFileTimeHi = 0U;
  c3_info[10].context =
    "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m";
  c3_info[10].name = "mtimes";
  c3_info[10].dominantType = "double";
  c3_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c3_info[10].fileTimeLo = 1289519692U;
  c3_info[10].fileTimeHi = 0U;
  c3_info[10].mFileTimeLo = 0U;
  c3_info[10].mFileTimeHi = 0U;
  c3_info[11].context =
    "[E]C:/Program Files/MATLAB/mss/gnc/gnc_mfiles/eulerang.m";
  c3_info[11].name = "mrdivide";
  c3_info[11].dominantType = "double";
  c3_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c3_info[11].fileTimeLo = 1357951548U;
  c3_info[11].fileTimeHi = 0U;
  c3_info[11].mFileTimeLo = 1319729966U;
  c3_info[11].mFileTimeHi = 0U;
  c3_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c3_info[12].name = "rdivide";
  c3_info[12].dominantType = "double";
  c3_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c3_info[12].fileTimeLo = 1346510388U;
  c3_info[12].fileTimeHi = 0U;
  c3_info[12].mFileTimeLo = 0U;
  c3_info[12].mFileTimeHi = 0U;
  c3_info[13].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c3_info[13].name = "eml_scalexp_compatible";
  c3_info[13].dominantType = "double";
  c3_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  c3_info[13].fileTimeLo = 1286818796U;
  c3_info[13].fileTimeHi = 0U;
  c3_info[13].mFileTimeLo = 0U;
  c3_info[13].mFileTimeHi = 0U;
  c3_info[14].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c3_info[14].name = "eml_div";
  c3_info[14].dominantType = "double";
  c3_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c3_info[14].fileTimeLo = 1313347810U;
  c3_info[14].fileTimeHi = 0U;
  c3_info[14].mFileTimeLo = 0U;
  c3_info[14].mFileTimeHi = 0U;
  c3_info[15].context = "";
  c3_info[15].name = "mtimes";
  c3_info[15].dominantType = "double";
  c3_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c3_info[15].fileTimeLo = 1289519692U;
  c3_info[15].fileTimeHi = 0U;
  c3_info[15].mFileTimeLo = 0U;
  c3_info[15].mFileTimeHi = 0U;
  c3_info[16].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c3_info[16].name = "eml_index_class";
  c3_info[16].dominantType = "";
  c3_info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c3_info[16].fileTimeLo = 1323170578U;
  c3_info[16].fileTimeHi = 0U;
  c3_info[16].mFileTimeLo = 0U;
  c3_info[16].mFileTimeHi = 0U;
  c3_info[17].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c3_info[17].name = "eml_scalar_eg";
  c3_info[17].dominantType = "double";
  c3_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c3_info[17].fileTimeLo = 1286818796U;
  c3_info[17].fileTimeHi = 0U;
  c3_info[17].mFileTimeLo = 0U;
  c3_info[17].mFileTimeHi = 0U;
  c3_info[18].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c3_info[18].name = "eml_xgemm";
  c3_info[18].dominantType = "char";
  c3_info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c3_info[18].fileTimeLo = 1299076772U;
  c3_info[18].fileTimeHi = 0U;
  c3_info[18].mFileTimeLo = 0U;
  c3_info[18].mFileTimeHi = 0U;
  c3_info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c3_info[19].name = "eml_blas_inline";
  c3_info[19].dominantType = "";
  c3_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c3_info[19].fileTimeLo = 1299076768U;
  c3_info[19].fileTimeHi = 0U;
  c3_info[19].mFileTimeLo = 0U;
  c3_info[19].mFileTimeHi = 0U;
  c3_info[20].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold";
  c3_info[20].name = "mtimes";
  c3_info[20].dominantType = "double";
  c3_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c3_info[20].fileTimeLo = 1289519692U;
  c3_info[20].fileTimeHi = 0U;
  c3_info[20].mFileTimeLo = 0U;
  c3_info[20].mFileTimeHi = 0U;
  c3_info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c3_info[21].name = "eml_index_class";
  c3_info[21].dominantType = "";
  c3_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c3_info[21].fileTimeLo = 1323170578U;
  c3_info[21].fileTimeHi = 0U;
  c3_info[21].mFileTimeLo = 0U;
  c3_info[21].mFileTimeHi = 0U;
  c3_info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c3_info[22].name = "eml_scalar_eg";
  c3_info[22].dominantType = "double";
  c3_info[22].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c3_info[22].fileTimeLo = 1286818796U;
  c3_info[22].fileTimeHi = 0U;
  c3_info[22].mFileTimeLo = 0U;
  c3_info[22].mFileTimeHi = 0U;
  c3_info[23].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c3_info[23].name = "eml_refblas_xgemm";
  c3_info[23].dominantType = "char";
  c3_info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c3_info[23].fileTimeLo = 1299076774U;
  c3_info[23].fileTimeHi = 0U;
  c3_info[23].mFileTimeLo = 0U;
  c3_info[23].mFileTimeHi = 0U;
}

static void c3_eml_scalar_eg(SFc3_testInstanceStruct *chartInstance)
{
}

static const mxArray *c3_e_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  int32_T c3_u;
  const mxArray *c3_y = NULL;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  c3_u = *(int32_T *)c3_inData;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c3_mxArrayOutData, c3_y, FALSE);
  return c3_mxArrayOutData;
}

static int32_T c3_f_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i62;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i62, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i62;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  const mxArray *c3_b_sfEvent;
  const char_T *c3_identifier;
  emlrtMsgIdentifier c3_thisId;
  int32_T c3_y;
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)chartInstanceVoid;
  c3_b_sfEvent = sf_mex_dup(c3_mxArrayInData);
  c3_identifier = c3_varName;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  *(int32_T *)c3_outData = c3_y;
  sf_mex_destroy(&c3_mxArrayInData);
}

static uint8_T c3_g_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_b_is_active_c3_test, const char_T *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_is_active_c3_test),
    &c3_thisId);
  sf_mex_destroy(&c3_b_is_active_c3_test);
  return c3_y;
}

static uint8_T c3_h_emlrt_marshallIn(SFc3_testInstanceStruct *chartInstance,
  const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void init_dsm_address_info(SFc3_testInstanceStruct *chartInstance)
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

void sf_c3_test_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1549329948U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(824723978U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1170484214U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3786438667U);
}

mxArray *sf_c3_test_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("nAH5NhuBjQecKrJgRstPuH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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

mxArray *sf_c3_test_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c3_test(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"eta_dot\",},{M[8],M[0],T\"is_active_c3_test\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_test_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_testInstanceStruct *chartInstance;
    chartInstance = (SFc3_testInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _testMachineNumber_,
           3,
           1,
           1,
           3,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"eta");
          _SFD_SET_DATA_PROPS(1,2,0,1,"eta_dot");
          _SFD_SET_DATA_PROPS(2,1,1,0,"nu");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,109);
        _SFD_CV_INIT_SCRIPT(0,1,1,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(0,0,"eulerang",0,-1,1507);
        _SFD_CV_INIT_SCRIPT_IF(0,0,1282,1291,-1,1346);
        _SFD_CV_INIT_SCRIPT(1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_SCRIPT_FCN(1,0,"Rzyx",0,-1,1478);
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
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)
            c3_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T (*c3_eta)[6];
          real_T (*c3_eta_dot)[6];
          real_T (*c3_nu)[6];
          c3_nu = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
          c3_eta_dot = (real_T (*)[6])ssGetOutputPortSignal(chartInstance->S, 1);
          c3_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c3_eta);
          _SFD_SET_DATA_VALUE_PTR(1U, *c3_eta_dot);
          _SFD_SET_DATA_VALUE_PTR(2U, *c3_nu);
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
  return "hcSUisGJFpCYJOcvnWIGrE";
}

static void sf_opaque_initialize_c3_test(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc3_testInstanceStruct*) chartInstanceVar)->S,0);
  initialize_params_c3_test((SFc3_testInstanceStruct*) chartInstanceVar);
  initialize_c3_test((SFc3_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c3_test(void *chartInstanceVar)
{
  enable_c3_test((SFc3_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c3_test(void *chartInstanceVar)
{
  disable_c3_test((SFc3_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c3_test(void *chartInstanceVar)
{
  sf_c3_test((SFc3_testInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c3_test(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c3_test((SFc3_testInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_test();/* state var info */
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

extern void sf_internal_set_sim_state_c3_test(SimStruct* S, const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_test();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c3_test((SFc3_testInstanceStruct*)chartInfo->chartInstance,
                        mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c3_test(SimStruct* S)
{
  return sf_internal_get_sim_state_c3_test(S);
}

static void sf_opaque_set_sim_state_c3_test(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c3_test(S, st);
}

static void sf_opaque_terminate_c3_test(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_testInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_optimization_info();
    }

    finalize_c3_test((SFc3_testInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_test((SFc3_testInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_test(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c3_test((SFc3_testInstanceStruct*)(((ChartInfoStruct *)
      ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,3,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3888210082U));
  ssSetChecksum1(S,(3798156016U));
  ssSetChecksum2(S,(4171902106U));
  ssSetChecksum3(S,(3730112902U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c3_test(SimStruct *S)
{
  SFc3_testInstanceStruct *chartInstance;
  chartInstance = (SFc3_testInstanceStruct *)utMalloc(sizeof
    (SFc3_testInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc3_testInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c3_test;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c3_test;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c3_test;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c3_test;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c3_test;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c3_test;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c3_test;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c3_test;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_test;
  chartInstance->chartInfo.mdlStart = mdlStart_c3_test;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c3_test;
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

void c3_test_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_test(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_test(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_test(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_test_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
