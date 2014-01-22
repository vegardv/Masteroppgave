/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_sfun.h"
#include "c5_test.h"
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
static const char * c5_debug_family_names[42] = { "x_c", "y_c", "phi", "theta",
  "psi", "height", "d1", "lengthPerPixel", "sinPhi", "cosPhi", "sinTheta",
  "cosTheta", "sinPsi", "cosPsi", "corners", "x", "y", "i", "deltaX", "deltaY",
  "alpha", "d", "beta", "sinAlpha", "cosAlpha", "sinBeta", "cosBeta", "z1", "z2",
  "z3", "d2", "result", "j", "nargin", "nargout", "nodePos", "eta",
  "measuredNodePos", "measuredNodeHeading", "inFrame", "cornersX", "cornersY" };

/* Function Declarations */
static void initialize_c5_test(SFc5_testInstanceStruct *chartInstance);
static void initialize_params_c5_test(SFc5_testInstanceStruct *chartInstance);
static void enable_c5_test(SFc5_testInstanceStruct *chartInstance);
static void disable_c5_test(SFc5_testInstanceStruct *chartInstance);
static void c5_update_debugger_state_c5_test(SFc5_testInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c5_test(SFc5_testInstanceStruct
  *chartInstance);
static void set_sim_state_c5_test(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_st);
static void finalize_c5_test(SFc5_testInstanceStruct *chartInstance);
static void sf_c5_test(SFc5_testInstanceStruct *chartInstance);
static void c5_chartstep_c5_test(SFc5_testInstanceStruct *chartInstance);
static void initSimStructsc5_test(SFc5_testInstanceStruct *chartInstance);
static void registerMessagesc5_test(SFc5_testInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber);
static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData);
static void c5_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_cornersY, const char_T *c5_identifier, real_T c5_y[4]);
static void c5_b_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[4]);
static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static boolean_T c5_c_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_inFrame, const char_T *c5_identifier);
static boolean_T c5_d_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static real_T c5_e_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_measuredNodeHeading, const char_T *c5_identifier);
static real_T c5_f_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static const mxArray *c5_d_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static const mxArray *c5_e_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static const mxArray *c5_f_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static void c5_g_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[8]);
static void c5_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static void c5_info_helper(c5_ResolvedFunctionInfo c5_info[28]);
static real_T c5_atan2(SFc5_testInstanceStruct *chartInstance, real_T c5_y,
  real_T c5_x);
static void c5_eml_scalar_eg(SFc5_testInstanceStruct *chartInstance);
static real_T c5_mpower(SFc5_testInstanceStruct *chartInstance, real_T c5_a);
static real_T c5_sqrt(SFc5_testInstanceStruct *chartInstance, real_T c5_x);
static void c5_eml_error(SFc5_testInstanceStruct *chartInstance);
static const mxArray *c5_g_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData);
static int32_T c5_h_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData);
static uint8_T c5_i_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_b_is_active_c5_test, const char_T *c5_identifier);
static uint8_T c5_j_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId);
static void c5_b_sqrt(SFc5_testInstanceStruct *chartInstance, real_T *c5_x);
static void init_dsm_address_info(SFc5_testInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c5_test(SFc5_testInstanceStruct *chartInstance)
{
  chartInstance->c5_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c5_is_active_c5_test = 0U;
}

static void initialize_params_c5_test(SFc5_testInstanceStruct *chartInstance)
{
}

static void enable_c5_test(SFc5_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c5_test(SFc5_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c5_update_debugger_state_c5_test(SFc5_testInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c5_test(SFc5_testInstanceStruct
  *chartInstance)
{
  const mxArray *c5_st;
  const mxArray *c5_y = NULL;
  int32_T c5_i0;
  real_T c5_u[4];
  const mxArray *c5_b_y = NULL;
  int32_T c5_i1;
  real_T c5_b_u[4];
  const mxArray *c5_c_y = NULL;
  boolean_T c5_hoistedGlobal;
  boolean_T c5_c_u;
  const mxArray *c5_d_y = NULL;
  real_T c5_b_hoistedGlobal;
  real_T c5_d_u;
  const mxArray *c5_e_y = NULL;
  real_T c5_c_hoistedGlobal;
  real_T c5_e_u;
  const mxArray *c5_f_y = NULL;
  uint8_T c5_d_hoistedGlobal;
  uint8_T c5_f_u;
  const mxArray *c5_g_y = NULL;
  boolean_T *c5_inFrame;
  real_T *c5_measuredNodeHeading;
  real_T *c5_measuredNodePos;
  real_T (*c5_cornersY)[4];
  real_T (*c5_cornersX)[4];
  c5_cornersY = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 5);
  c5_cornersX = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c5_inFrame = (boolean_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c5_measuredNodeHeading = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c5_measuredNodePos = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c5_st = NULL;
  c5_st = NULL;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_createcellarray(6), FALSE);
  for (c5_i0 = 0; c5_i0 < 4; c5_i0++) {
    c5_u[c5_i0] = (*c5_cornersX)[c5_i0];
  }

  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_setcell(c5_y, 0, c5_b_y);
  for (c5_i1 = 0; c5_i1 < 4; c5_i1++) {
    c5_b_u[c5_i1] = (*c5_cornersY)[c5_i1];
  }

  c5_c_y = NULL;
  sf_mex_assign(&c5_c_y, sf_mex_create("y", c5_b_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_setcell(c5_y, 1, c5_c_y);
  c5_hoistedGlobal = *c5_inFrame;
  c5_c_u = c5_hoistedGlobal;
  c5_d_y = NULL;
  sf_mex_assign(&c5_d_y, sf_mex_create("y", &c5_c_u, 11, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 2, c5_d_y);
  c5_b_hoistedGlobal = *c5_measuredNodeHeading;
  c5_d_u = c5_b_hoistedGlobal;
  c5_e_y = NULL;
  sf_mex_assign(&c5_e_y, sf_mex_create("y", &c5_d_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 3, c5_e_y);
  c5_c_hoistedGlobal = *c5_measuredNodePos;
  c5_e_u = c5_c_hoistedGlobal;
  c5_f_y = NULL;
  sf_mex_assign(&c5_f_y, sf_mex_create("y", &c5_e_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 4, c5_f_y);
  c5_d_hoistedGlobal = chartInstance->c5_is_active_c5_test;
  c5_f_u = c5_d_hoistedGlobal;
  c5_g_y = NULL;
  sf_mex_assign(&c5_g_y, sf_mex_create("y", &c5_f_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c5_y, 5, c5_g_y);
  sf_mex_assign(&c5_st, c5_y, FALSE);
  return c5_st;
}

static void set_sim_state_c5_test(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_st)
{
  const mxArray *c5_u;
  real_T c5_dv0[4];
  int32_T c5_i2;
  real_T c5_dv1[4];
  int32_T c5_i3;
  boolean_T *c5_inFrame;
  real_T *c5_measuredNodeHeading;
  real_T *c5_measuredNodePos;
  real_T (*c5_cornersX)[4];
  real_T (*c5_cornersY)[4];
  c5_cornersY = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 5);
  c5_cornersX = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c5_inFrame = (boolean_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c5_measuredNodeHeading = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c5_measuredNodePos = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c5_doneDoubleBufferReInit = TRUE;
  c5_u = sf_mex_dup(c5_st);
  c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 0)),
                      "cornersX", c5_dv0);
  for (c5_i2 = 0; c5_i2 < 4; c5_i2++) {
    (*c5_cornersX)[c5_i2] = c5_dv0[c5_i2];
  }

  c5_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c5_u, 1)),
                      "cornersY", c5_dv1);
  for (c5_i3 = 0; c5_i3 < 4; c5_i3++) {
    (*c5_cornersY)[c5_i3] = c5_dv1[c5_i3];
  }

  *c5_inFrame = c5_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c5_u, 2)), "inFrame");
  *c5_measuredNodeHeading = c5_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 3)), "measuredNodeHeading");
  *c5_measuredNodePos = c5_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c5_u, 4)), "measuredNodePos");
  chartInstance->c5_is_active_c5_test = c5_i_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c5_u, 5)), "is_active_c5_test");
  sf_mex_destroy(&c5_u);
  c5_update_debugger_state_c5_test(chartInstance);
  sf_mex_destroy(&c5_st);
}

static void finalize_c5_test(SFc5_testInstanceStruct *chartInstance)
{
}

static void sf_c5_test(SFc5_testInstanceStruct *chartInstance)
{
  int32_T c5_i4;
  int32_T c5_i5;
  int32_T c5_i6;
  int32_T c5_i7;
  real_T *c5_measuredNodePos;
  real_T *c5_measuredNodeHeading;
  boolean_T *c5_inFrame;
  real_T (*c5_cornersY)[4];
  real_T (*c5_cornersX)[4];
  real_T (*c5_eta)[6];
  real_T (*c5_nodePos)[3];
  c5_cornersY = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 5);
  c5_cornersX = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c5_inFrame = (boolean_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c5_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c5_measuredNodeHeading = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c5_measuredNodePos = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c5_nodePos = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  for (c5_i4 = 0; c5_i4 < 3; c5_i4++) {
    _SFD_DATA_RANGE_CHECK((*c5_nodePos)[c5_i4], 0U);
  }

  _SFD_DATA_RANGE_CHECK(*c5_measuredNodePos, 1U);
  _SFD_DATA_RANGE_CHECK(*c5_measuredNodeHeading, 2U);
  for (c5_i5 = 0; c5_i5 < 6; c5_i5++) {
    _SFD_DATA_RANGE_CHECK((*c5_eta)[c5_i5], 3U);
  }

  _SFD_DATA_RANGE_CHECK((real_T)*c5_inFrame, 4U);
  for (c5_i6 = 0; c5_i6 < 4; c5_i6++) {
    _SFD_DATA_RANGE_CHECK((*c5_cornersX)[c5_i6], 5U);
  }

  for (c5_i7 = 0; c5_i7 < 4; c5_i7++) {
    _SFD_DATA_RANGE_CHECK((*c5_cornersY)[c5_i7], 6U);
  }

  chartInstance->c5_sfEvent = CALL_EVENT;
  c5_chartstep_c5_test(chartInstance);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_testMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c5_chartstep_c5_test(SFc5_testInstanceStruct *chartInstance)
{
  int32_T c5_i8;
  real_T c5_nodePos[3];
  int32_T c5_i9;
  real_T c5_eta[6];
  uint32_T c5_debug_family_var_map[42];
  real_T c5_x_c;
  real_T c5_y_c;
  real_T c5_phi;
  real_T c5_theta;
  real_T c5_psi;
  real_T c5_height;
  real_T c5_d1;
  real_T c5_lengthPerPixel;
  real_T c5_sinPhi;
  real_T c5_cosPhi;
  real_T c5_sinTheta;
  real_T c5_cosTheta;
  real_T c5_sinPsi;
  real_T c5_cosPsi;
  real_T c5_corners[8];
  real_T c5_x;
  real_T c5_y;
  real_T c5_i;
  real_T c5_deltaX;
  real_T c5_deltaY;
  real_T c5_alpha;
  real_T c5_d;
  real_T c5_beta;
  real_T c5_sinAlpha;
  real_T c5_cosAlpha;
  real_T c5_sinBeta;
  real_T c5_cosBeta;
  real_T c5_z1;
  real_T c5_z2;
  real_T c5_z3;
  real_T c5_d2;
  boolean_T c5_result;
  real_T c5_j;
  real_T c5_nargin = 2.0;
  real_T c5_nargout = 5.0;
  real_T c5_measuredNodePos;
  real_T c5_measuredNodeHeading;
  boolean_T c5_inFrame;
  real_T c5_cornersX[4];
  real_T c5_cornersY[4];
  real_T c5_b_x;
  real_T c5_c_x;
  real_T c5_d_x;
  real_T c5_e_x;
  real_T c5_f_x;
  real_T c5_g_x;
  real_T c5_h_x;
  real_T c5_i_x;
  real_T c5_j_x;
  real_T c5_k_x;
  real_T c5_l_x;
  real_T c5_m_x;
  int32_T c5_i10;
  int32_T c5_b_i;
  real_T c5_a;
  real_T c5_n_x;
  real_T c5_o_x;
  real_T c5_p_x;
  real_T c5_q_x;
  real_T c5_r_x;
  real_T c5_s_x;
  real_T c5_t_x;
  real_T c5_u_x;
  real_T c5_v_x;
  real_T c5_w_x;
  real_T c5_b_a;
  real_T c5_b;
  real_T c5_b_y;
  real_T c5_c_a;
  real_T c5_b_b;
  real_T c5_c_y;
  real_T c5_d_a;
  real_T c5_c_b;
  real_T c5_d_y;
  real_T c5_e_a;
  real_T c5_d_b;
  real_T c5_e_y;
  real_T c5_f_a;
  real_T c5_e_b;
  real_T c5_f_y;
  real_T c5_g_a;
  real_T c5_f_b;
  real_T c5_g_y;
  real_T c5_h_a;
  real_T c5_g_b;
  real_T c5_h_y;
  real_T c5_i_a;
  real_T c5_h_b;
  real_T c5_i_y;
  real_T c5_j_a;
  real_T c5_i_b;
  real_T c5_j_y;
  real_T c5_k_a;
  real_T c5_j_b;
  real_T c5_k_y;
  real_T c5_l_a;
  real_T c5_k_b;
  real_T c5_l_y;
  real_T c5_m_a;
  real_T c5_l_b;
  real_T c5_m_y;
  real_T c5_n_a;
  real_T c5_m_b;
  real_T c5_n_y;
  real_T c5_o_a;
  real_T c5_n_b;
  real_T c5_o_y;
  real_T c5_p_a;
  real_T c5_o_b;
  real_T c5_p_y;
  real_T c5_q_a;
  real_T c5_p_b;
  real_T c5_q_y;
  real_T c5_r_a;
  real_T c5_q_b;
  real_T c5_r_y;
  real_T c5_s_a;
  real_T c5_r_b;
  real_T c5_s_y;
  real_T c5_t_a;
  real_T c5_s_b;
  real_T c5_t_y;
  real_T c5_u_a;
  real_T c5_t_b;
  real_T c5_u_y;
  real_T c5_v_a;
  real_T c5_u_b;
  real_T c5_v_y;
  real_T c5_w_a;
  real_T c5_v_b;
  real_T c5_w_y;
  real_T c5_x_a;
  real_T c5_w_b;
  real_T c5_x_y;
  real_T c5_y_a;
  real_T c5_x_b;
  real_T c5_y_y;
  real_T c5_ab_a;
  real_T c5_y_b;
  real_T c5_ab_y;
  real_T c5_bb_a;
  real_T c5_ab_b;
  real_T c5_bb_y;
  real_T c5_cb_a;
  real_T c5_bb_b;
  real_T c5_cb_y;
  real_T c5_db_a;
  real_T c5_cb_b;
  real_T c5_db_y;
  real_T c5_eb_a;
  real_T c5_db_b;
  real_T c5_eb_y;
  real_T c5_fb_a;
  real_T c5_eb_b;
  real_T c5_fb_y;
  real_T c5_gb_a;
  real_T c5_fb_b;
  real_T c5_gb_y;
  real_T c5_x_x;
  real_T c5_y_x;
  real_T c5_hb_y;
  real_T c5_gb_b;
  real_T c5_ib_y;
  real_T c5_hb_a;
  real_T c5_hb_b;
  real_T c5_jb_y;
  real_T c5_A;
  real_T c5_B;
  real_T c5_ab_x;
  real_T c5_kb_y;
  real_T c5_bb_x;
  real_T c5_lb_y;
  real_T c5_ib_a;
  real_T c5_ib_b;
  real_T c5_mb_y;
  real_T c5_jb_a;
  real_T c5_jb_b;
  real_T c5_nb_y;
  real_T c5_kb_a;
  real_T c5_kb_b;
  real_T c5_ob_y;
  real_T c5_lb_a;
  real_T c5_lb_b;
  real_T c5_pb_y;
  real_T c5_mb_b;
  real_T c5_qb_y;
  real_T c5_mb_a;
  real_T c5_nb_b;
  real_T c5_rb_y;
  real_T c5_nb_a;
  real_T c5_ob_b;
  real_T c5_sb_y;
  real_T c5_ob_a;
  real_T c5_pb_b;
  real_T c5_tb_y;
  real_T c5_pb_a;
  real_T c5_qb_b;
  real_T c5_ub_y;
  real_T c5_rb_b;
  real_T c5_vb_y;
  int32_T c5_c_i;
  real_T c5_qb_a;
  real_T c5_sb_b;
  real_T c5_wb_y;
  real_T c5_b_A;
  real_T c5_b_B;
  real_T c5_cb_x;
  real_T c5_xb_y;
  real_T c5_db_x;
  real_T c5_yb_y;
  real_T c5_ac_y;
  int32_T c5_i11;
  int32_T c5_i12;
  int32_T c5_i13;
  int32_T c5_i14;
  boolean_T *c5_b_inFrame;
  real_T *c5_b_measuredNodeHeading;
  real_T *c5_b_measuredNodePos;
  real_T (*c5_b_cornersX)[4];
  real_T (*c5_b_cornersY)[4];
  real_T (*c5_b_eta)[6];
  real_T (*c5_b_nodePos)[3];
  boolean_T guard1 = FALSE;
  c5_b_cornersY = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 5);
  c5_b_cornersX = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c5_b_inFrame = (boolean_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c5_b_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
  c5_b_measuredNodeHeading = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c5_b_measuredNodePos = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c5_b_nodePos = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
  for (c5_i8 = 0; c5_i8 < 3; c5_i8++) {
    c5_nodePos[c5_i8] = (*c5_b_nodePos)[c5_i8];
  }

  for (c5_i9 = 0; c5_i9 < 6; c5_i9++) {
    c5_eta[c5_i9] = (*c5_b_eta)[c5_i9];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 42U, 42U, c5_debug_family_names,
    c5_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_x_c, 0U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_y_c, 1U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_phi, 2U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_theta, 3U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_psi, 4U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_height, 5U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_d1, 6U, c5_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c5_lengthPerPixel, 7U, c5_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_sinPhi, 8U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_cosPhi, 9U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_sinTheta, 10U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_cosTheta, 11U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_sinPsi, 12U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_cosPsi, 13U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_corners, 14U, c5_f_sf_marshallOut,
    c5_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_x, 15U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_y, 16U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_i, 17U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_deltaX, 18U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_deltaY, 19U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_alpha, 20U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_d, 21U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_beta, 22U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_sinAlpha, 23U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_cosAlpha, 24U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_sinBeta, 25U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_cosBeta, 26U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_z1, 27U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_z2, 28U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_z3, 29U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_d2, 30U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_result, 31U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_j, 32U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargin, 33U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_nargout, 34U, c5_c_sf_marshallOut,
    c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c5_nodePos, 35U, c5_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c5_eta, 36U, c5_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_measuredNodePos, 37U,
    c5_c_sf_marshallOut, c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_measuredNodeHeading, 38U,
    c5_c_sf_marshallOut, c5_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c5_inFrame, 39U, c5_b_sf_marshallOut,
    c5_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_cornersX, 40U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c5_cornersY, 41U, c5_sf_marshallOut,
    c5_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 7);
  c5_x_c = 320.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 8);
  c5_y_c = 240.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 10);
  c5_phi = c5_eta[3];
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 11);
  c5_theta = c5_eta[4];
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 12);
  c5_psi = c5_eta[5];
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 13);
  c5_height = c5_eta[2];
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 15);
  c5_d1 = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 16);
  c5_lengthPerPixel = 0.001328497723204;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 18);
  c5_b_x = c5_phi;
  c5_sinPhi = c5_b_x;
  c5_c_x = c5_sinPhi;
  c5_sinPhi = c5_c_x;
  c5_sinPhi = muDoubleScalarSin(c5_sinPhi);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 19);
  c5_d_x = c5_phi;
  c5_cosPhi = c5_d_x;
  c5_e_x = c5_cosPhi;
  c5_cosPhi = c5_e_x;
  c5_cosPhi = muDoubleScalarCos(c5_cosPhi);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 21);
  c5_f_x = c5_theta;
  c5_sinTheta = c5_f_x;
  c5_g_x = c5_sinTheta;
  c5_sinTheta = c5_g_x;
  c5_sinTheta = muDoubleScalarSin(c5_sinTheta);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 22);
  c5_h_x = c5_theta;
  c5_cosTheta = c5_h_x;
  c5_i_x = c5_cosTheta;
  c5_cosTheta = c5_i_x;
  c5_cosTheta = muDoubleScalarCos(c5_cosTheta);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 24);
  c5_j_x = c5_psi;
  c5_sinPsi = c5_j_x;
  c5_k_x = c5_sinPsi;
  c5_sinPsi = c5_k_x;
  c5_sinPsi = muDoubleScalarSin(c5_sinPsi);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 25);
  c5_l_x = c5_psi;
  c5_cosPsi = c5_l_x;
  c5_m_x = c5_cosPsi;
  c5_cosPsi = c5_m_x;
  c5_cosPsi = muDoubleScalarCos(c5_cosPsi);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 27);
  for (c5_i10 = 0; c5_i10 < 8; c5_i10++) {
    c5_corners[c5_i10] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 28);
  c5_x = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 29);
  c5_y = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 30);
  c5_i = 1.0;
  c5_b_i = 0;
  while (c5_b_i < 4) {
    c5_i = 1.0 + (real_T)c5_b_i;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 31);
    switch ((int32_T)_SFD_INTEGER_CHECK("i", c5_i)) {
     case 1:
      CV_EML_SWITCH(0, 1, 0, 1);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 33);
      c5_x = 0.0;
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 34);
      c5_y = 0.0;
      break;

     case 2:
      CV_EML_SWITCH(0, 1, 0, 2);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 36);
      c5_x = 640.0;
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 37);
      c5_y = 0.0;
      break;

     case 3:
      CV_EML_SWITCH(0, 1, 0, 3);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 39);
      c5_x = 640.0;
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 40);
      c5_y = 480.0;
      break;

     case 4:
      CV_EML_SWITCH(0, 1, 0, 4);
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 42);
      c5_x = 0.0;
      _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 43);
      c5_y = 480.0;
      break;

     default:
      CV_EML_SWITCH(0, 1, 0, 0);
      break;
    }

    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 46);
    c5_deltaX = c5_x - c5_x_c;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 47);
    c5_deltaY = c5_y - c5_y_c;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 49);
    c5_alpha = -c5_atan2(chartInstance, c5_deltaY, -c5_deltaX);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 50);
    c5_a = c5_mpower(chartInstance, c5_deltaX) + c5_mpower(chartInstance,
      c5_deltaY);
    c5_b_sqrt(chartInstance, &c5_a);
    c5_d = c5_a * 0.001328497723204;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 51);
    c5_n_x = c5_d;
    c5_o_x = c5_n_x;
    c5_o_x = muDoubleScalarAtan(c5_o_x);
    c5_beta = -c5_o_x;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 52);
    c5_p_x = c5_alpha;
    c5_sinAlpha = c5_p_x;
    c5_q_x = c5_sinAlpha;
    c5_sinAlpha = c5_q_x;
    c5_sinAlpha = muDoubleScalarSin(c5_sinAlpha);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 53);
    c5_r_x = c5_alpha;
    c5_cosAlpha = c5_r_x;
    c5_s_x = c5_cosAlpha;
    c5_cosAlpha = c5_s_x;
    c5_cosAlpha = muDoubleScalarCos(c5_cosAlpha);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 54);
    c5_t_x = c5_beta;
    c5_sinBeta = c5_t_x;
    c5_u_x = c5_sinBeta;
    c5_sinBeta = c5_u_x;
    c5_sinBeta = muDoubleScalarSin(c5_sinBeta);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 55);
    c5_v_x = c5_beta;
    c5_cosBeta = c5_v_x;
    c5_w_x = c5_cosBeta;
    c5_cosBeta = c5_w_x;
    c5_cosBeta = muDoubleScalarCos(c5_cosBeta);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 57);
    c5_b_a = c5_sinAlpha;
    c5_b = c5_sinBeta;
    c5_b_y = c5_b_a * c5_b;
    c5_c_a = c5_b_y;
    c5_b_b = c5_cosPsi;
    c5_c_y = c5_c_a * c5_b_b;
    c5_d_a = c5_c_y;
    c5_c_b = c5_cosTheta;
    c5_d_y = c5_d_a * c5_c_b;
    c5_e_a = c5_cosAlpha;
    c5_d_b = c5_sinBeta;
    c5_e_y = c5_e_a * c5_d_b;
    c5_f_a = -c5_sinPsi;
    c5_e_b = c5_cosPhi;
    c5_f_y = c5_f_a * c5_e_b;
    c5_g_a = c5_cosPsi;
    c5_f_b = c5_sinTheta;
    c5_g_y = c5_g_a * c5_f_b;
    c5_h_a = c5_g_y;
    c5_g_b = c5_sinPhi;
    c5_h_y = c5_h_a * c5_g_b;
    c5_i_a = c5_e_y;
    c5_h_b = c5_f_y + c5_h_y;
    c5_i_y = c5_i_a * c5_h_b;
    c5_j_a = c5_sinPsi;
    c5_i_b = c5_sinPhi;
    c5_j_y = c5_j_a * c5_i_b;
    c5_k_a = c5_cosPsi;
    c5_j_b = c5_sinTheta;
    c5_k_y = c5_k_a * c5_j_b;
    c5_l_a = c5_k_y;
    c5_k_b = c5_cosPhi;
    c5_l_y = c5_l_a * c5_k_b;
    c5_m_a = c5_cosBeta;
    c5_l_b = c5_j_y + c5_l_y;
    c5_m_y = c5_m_a * c5_l_b;
    c5_z1 = (c5_d_y - c5_i_y) + c5_m_y;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 58);
    c5_n_a = c5_sinAlpha;
    c5_m_b = c5_sinBeta;
    c5_n_y = c5_n_a * c5_m_b;
    c5_o_a = c5_n_y;
    c5_n_b = c5_sinPsi;
    c5_o_y = c5_o_a * c5_n_b;
    c5_p_a = c5_o_y;
    c5_o_b = c5_cosTheta;
    c5_p_y = c5_p_a * c5_o_b;
    c5_q_a = c5_cosAlpha;
    c5_p_b = c5_sinBeta;
    c5_q_y = c5_q_a * c5_p_b;
    c5_r_a = c5_cosPsi;
    c5_q_b = c5_cosPhi;
    c5_r_y = c5_r_a * c5_q_b;
    c5_s_a = c5_sinPsi;
    c5_r_b = c5_sinTheta;
    c5_s_y = c5_s_a * c5_r_b;
    c5_t_a = c5_s_y;
    c5_s_b = c5_sinPhi;
    c5_t_y = c5_t_a * c5_s_b;
    c5_u_a = c5_q_y;
    c5_t_b = c5_r_y + c5_t_y;
    c5_u_y = c5_u_a * c5_t_b;
    c5_v_a = -c5_cosPsi;
    c5_u_b = c5_sinPhi;
    c5_v_y = c5_v_a * c5_u_b;
    c5_w_a = c5_sinPsi;
    c5_v_b = c5_sinTheta;
    c5_w_y = c5_w_a * c5_v_b;
    c5_x_a = c5_w_y;
    c5_w_b = c5_sinPhi;
    c5_x_y = c5_x_a * c5_w_b;
    c5_y_a = c5_cosBeta;
    c5_x_b = c5_v_y + c5_x_y;
    c5_y_y = c5_y_a * c5_x_b;
    c5_z2 = (c5_p_y - c5_u_y) + c5_y_y;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 59);
    c5_ab_a = -c5_sinAlpha;
    c5_y_b = c5_sinBeta;
    c5_ab_y = c5_ab_a * c5_y_b;
    c5_bb_a = c5_ab_y;
    c5_ab_b = c5_sinTheta;
    c5_bb_y = c5_bb_a * c5_ab_b;
    c5_cb_a = c5_cosAlpha;
    c5_bb_b = c5_sinBeta;
    c5_cb_y = c5_cb_a * c5_bb_b;
    c5_db_a = c5_cb_y;
    c5_cb_b = c5_cosTheta;
    c5_db_y = c5_db_a * c5_cb_b;
    c5_eb_a = c5_db_y;
    c5_db_b = c5_sinPhi;
    c5_eb_y = c5_eb_a * c5_db_b;
    c5_fb_a = c5_cosBeta;
    c5_eb_b = c5_cosTheta;
    c5_fb_y = c5_fb_a * c5_eb_b;
    c5_gb_a = c5_fb_y;
    c5_fb_b = c5_cosPhi;
    c5_gb_y = c5_gb_a * c5_fb_b;
    c5_z3 = (c5_bb_y - c5_eb_y) + c5_gb_y;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 61);
    c5_x_x = c5_height;
    c5_y_x = c5_x_x;
    c5_hb_y = muDoubleScalarAbs(c5_y_x);
    c5_gb_b = c5_cosTheta;
    c5_ib_y = 0.0 * c5_gb_b;
    c5_hb_a = c5_ib_y;
    c5_hb_b = c5_cosPhi;
    c5_jb_y = c5_hb_a * c5_hb_b;
    c5_A = c5_hb_y - c5_jb_y;
    c5_B = c5_z3;
    c5_ab_x = c5_A;
    c5_kb_y = c5_B;
    c5_bb_x = c5_ab_x;
    c5_lb_y = c5_kb_y;
    c5_d2 = c5_bb_x / c5_lb_y;
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 63);
    c5_ib_a = c5_d2;
    c5_ib_b = c5_z1;
    c5_mb_y = c5_ib_a * c5_ib_b;
    c5_jb_a = c5_sinPsi;
    c5_jb_b = c5_sinPhi;
    c5_nb_y = c5_jb_a * c5_jb_b;
    c5_kb_a = c5_cosPsi;
    c5_kb_b = c5_sinTheta;
    c5_ob_y = c5_kb_a * c5_kb_b;
    c5_lb_a = c5_ob_y;
    c5_lb_b = c5_cosPhi;
    c5_pb_y = c5_lb_a * c5_lb_b;
    c5_mb_b = c5_nb_y + c5_pb_y;
    c5_qb_y = 0.0 * c5_mb_b;
    c5_corners[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("corners", (int32_T)
      _SFD_INTEGER_CHECK("i", c5_i), 1, 4, 1, 0) - 1] = (c5_mb_y + c5_qb_y) +
      c5_eta[0];
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 64);
    c5_mb_a = c5_d2;
    c5_nb_b = c5_z2;
    c5_rb_y = c5_mb_a * c5_nb_b;
    c5_nb_a = -c5_cosPsi;
    c5_ob_b = c5_sinPhi;
    c5_sb_y = c5_nb_a * c5_ob_b;
    c5_ob_a = c5_sinPsi;
    c5_pb_b = c5_sinTheta;
    c5_tb_y = c5_ob_a * c5_pb_b;
    c5_pb_a = c5_tb_y;
    c5_qb_b = c5_cosPhi;
    c5_ub_y = c5_pb_a * c5_qb_b;
    c5_rb_b = c5_sb_y + c5_ub_y;
    c5_vb_y = 0.0 * c5_rb_b;
    c5_corners[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("corners", (int32_T)
      _SFD_INTEGER_CHECK("i", c5_i), 1, 4, 1, 0) + 3] = (c5_rb_y + c5_vb_y) +
      c5_eta[1];
    c5_b_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 66);
  c5_result = FALSE;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 67);
  c5_j = 4.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 68);
  c5_i = 1.0;
  c5_c_i = 0;
  while (c5_c_i < 4) {
    c5_i = 1.0 + (real_T)c5_c_i;
    CV_EML_FOR(0, 1, 1, 1);
    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 69);
    guard1 = FALSE;
    if (CV_EML_COND(0, 1, 0, c5_corners[(int32_T)(real_T)
                    _SFD_EML_ARRAY_BOUNDS_CHECK("corners", (int32_T)
          _SFD_INTEGER_CHECK("i", c5_i), 1, 4, 1, 0) + 3] > c5_nodePos[1] !=
                    c5_corners[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
          "corners", (int32_T)_SFD_INTEGER_CHECK("j", c5_j), 1, 4, 1, 0) + 3] >
                    c5_nodePos[1])) {
      c5_qb_a = c5_corners[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK(
        "corners", (int32_T)_SFD_INTEGER_CHECK("j", c5_j), 1, 4, 1, 0) - 1] -
        c5_corners[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("corners",
        (int32_T)_SFD_INTEGER_CHECK("i", c5_i), 1, 4, 1, 0) - 1];
      c5_sb_b = c5_nodePos[1] - c5_corners[(int32_T)(real_T)
        _SFD_EML_ARRAY_BOUNDS_CHECK("corners", (int32_T)_SFD_INTEGER_CHECK("i",
        c5_i), 1, 4, 1, 0) + 3];
      c5_wb_y = c5_qb_a * c5_sb_b;
      c5_b_A = c5_wb_y;
      c5_b_B = c5_corners[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("corners",
        (int32_T)_SFD_INTEGER_CHECK("j", c5_j), 1, 4, 1, 0) + 3] - c5_corners
        [(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("corners", (int32_T)
        _SFD_INTEGER_CHECK("i", c5_i), 1, 4, 1, 0) + 3];
      c5_cb_x = c5_b_A;
      c5_xb_y = c5_b_B;
      c5_db_x = c5_cb_x;
      c5_yb_y = c5_xb_y;
      c5_ac_y = c5_db_x / c5_yb_y;
      if (CV_EML_COND(0, 1, 1, c5_nodePos[0] < c5_ac_y + c5_corners[(int32_T)
                      (real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("corners", (int32_T)
            _SFD_INTEGER_CHECK("i", c5_i), 1, 4, 1, 0) - 1])) {
        CV_EML_MCDC(0, 1, 0, TRUE);
        CV_EML_IF(0, 1, 0, TRUE);
        _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 70);
        c5_result = !c5_result;
      } else {
        guard1 = TRUE;
      }
    } else {
      guard1 = TRUE;
    }

    if (guard1 == TRUE) {
      CV_EML_MCDC(0, 1, 0, FALSE);
      CV_EML_IF(0, 1, 0, FALSE);
    }

    _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 72);
    c5_j = c5_i;
    c5_c_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 1, 0);
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 74);
  c5_measuredNodePos = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 75);
  c5_measuredNodeHeading = 1.0;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 76);
  c5_inFrame = c5_result;
  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 77);
  for (c5_i11 = 0; c5_i11 < 4; c5_i11++) {
    c5_cornersX[c5_i11] = c5_corners[c5_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, 78);
  for (c5_i12 = 0; c5_i12 < 4; c5_i12++) {
    c5_cornersY[c5_i12] = c5_corners[c5_i12 + 4];
  }

  _SFD_EML_CALL(0U, chartInstance->c5_sfEvent, -78);
  _SFD_SYMBOL_SCOPE_POP();
  *c5_b_measuredNodePos = c5_measuredNodePos;
  *c5_b_measuredNodeHeading = c5_measuredNodeHeading;
  *c5_b_inFrame = c5_inFrame;
  for (c5_i13 = 0; c5_i13 < 4; c5_i13++) {
    (*c5_b_cornersX)[c5_i13] = c5_cornersX[c5_i13];
  }

  for (c5_i14 = 0; c5_i14 < 4; c5_i14++) {
    (*c5_b_cornersY)[c5_i14] = c5_cornersY[c5_i14];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 4U, chartInstance->c5_sfEvent);
}

static void initSimStructsc5_test(SFc5_testInstanceStruct *chartInstance)
{
}

static void registerMessagesc5_test(SFc5_testInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c5_machineNumber, uint32_T
  c5_chartNumber)
{
}

static const mxArray *c5_sf_marshallOut(void *chartInstanceVoid, void *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i15;
  real_T c5_b_inData[4];
  int32_T c5_i16;
  real_T c5_u[4];
  const mxArray *c5_y = NULL;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  for (c5_i15 = 0; c5_i15 < 4; c5_i15++) {
    c5_b_inData[c5_i15] = (*(real_T (*)[4])c5_inData)[c5_i15];
  }

  for (c5_i16 = 0; c5_i16 < 4; c5_i16++) {
    c5_u[c5_i16] = c5_b_inData[c5_i16];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, FALSE);
  return c5_mxArrayOutData;
}

static void c5_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_cornersY, const char_T *c5_identifier, real_T c5_y[4])
{
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_cornersY), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_cornersY);
}

static void c5_b_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[4])
{
  real_T c5_dv2[4];
  int32_T c5_i17;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv2, 1, 0, 0U, 1, 0U, 1, 4);
  for (c5_i17 = 0; c5_i17 < 4; c5_i17++) {
    c5_y[c5_i17] = c5_dv2[c5_i17];
  }

  sf_mex_destroy(&c5_u);
}

static void c5_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_cornersY;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y[4];
  int32_T c5_i18;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_cornersY = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_cornersY), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_cornersY);
  for (c5_i18 = 0; c5_i18 < 4; c5_i18++) {
    (*(real_T (*)[4])c5_outData)[c5_i18] = c5_y[c5_i18];
  }

  sf_mex_destroy(&c5_mxArrayInData);
}

static const mxArray *c5_b_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  boolean_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(boolean_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 11, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, FALSE);
  return c5_mxArrayOutData;
}

static boolean_T c5_c_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_inFrame, const char_T *c5_identifier)
{
  boolean_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_inFrame), &c5_thisId);
  sf_mex_destroy(&c5_inFrame);
  return c5_y;
}

static boolean_T c5_d_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  boolean_T c5_y;
  boolean_T c5_b0;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_b0, 1, 11, 0U, 0, 0U, 0);
  c5_y = c5_b0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_inFrame;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  boolean_T c5_y;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_inFrame = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_inFrame), &c5_thisId);
  sf_mex_destroy(&c5_inFrame);
  *(boolean_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static const mxArray *c5_c_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  real_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(real_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, FALSE);
  return c5_mxArrayOutData;
}

static real_T c5_e_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_measuredNodeHeading, const char_T *c5_identifier)
{
  real_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_measuredNodeHeading),
    &c5_thisId);
  sf_mex_destroy(&c5_measuredNodeHeading);
  return c5_y;
}

static real_T c5_f_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  real_T c5_y;
  real_T c5_d0;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_d0, 1, 0, 0U, 0, 0U, 0);
  c5_y = c5_d0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_measuredNodeHeading;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_measuredNodeHeading = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_measuredNodeHeading),
    &c5_thisId);
  sf_mex_destroy(&c5_measuredNodeHeading);
  *(real_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static const mxArray *c5_d_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i19;
  real_T c5_b_inData[6];
  int32_T c5_i20;
  real_T c5_u[6];
  const mxArray *c5_y = NULL;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  for (c5_i19 = 0; c5_i19 < 6; c5_i19++) {
    c5_b_inData[c5_i19] = (*(real_T (*)[6])c5_inData)[c5_i19];
  }

  for (c5_i20 = 0; c5_i20 < 6; c5_i20++) {
    c5_u[c5_i20] = c5_b_inData[c5_i20];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 1, 6), FALSE);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, FALSE);
  return c5_mxArrayOutData;
}

static const mxArray *c5_e_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i21;
  real_T c5_b_inData[3];
  int32_T c5_i22;
  real_T c5_u[3];
  const mxArray *c5_y = NULL;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  for (c5_i21 = 0; c5_i21 < 3; c5_i21++) {
    c5_b_inData[c5_i21] = (*(real_T (*)[3])c5_inData)[c5_i21];
  }

  for (c5_i22 = 0; c5_i22 < 3; c5_i22++) {
    c5_u[c5_i22] = c5_b_inData[c5_i22];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, FALSE);
  return c5_mxArrayOutData;
}

static const mxArray *c5_f_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_i23;
  int32_T c5_i24;
  int32_T c5_i25;
  real_T c5_b_inData[8];
  int32_T c5_i26;
  int32_T c5_i27;
  int32_T c5_i28;
  real_T c5_u[8];
  const mxArray *c5_y = NULL;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_i23 = 0;
  for (c5_i24 = 0; c5_i24 < 2; c5_i24++) {
    for (c5_i25 = 0; c5_i25 < 4; c5_i25++) {
      c5_b_inData[c5_i25 + c5_i23] = (*(real_T (*)[8])c5_inData)[c5_i25 + c5_i23];
    }

    c5_i23 += 4;
  }

  c5_i26 = 0;
  for (c5_i27 = 0; c5_i27 < 2; c5_i27++) {
    for (c5_i28 = 0; c5_i28 < 4; c5_i28++) {
      c5_u[c5_i28 + c5_i26] = c5_b_inData[c5_i28 + c5_i26];
    }

    c5_i26 += 4;
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 0, 0U, 1U, 0U, 2, 4, 2), FALSE);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, FALSE);
  return c5_mxArrayOutData;
}

static void c5_g_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance, const
  mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId, real_T c5_y[8])
{
  real_T c5_dv3[8];
  int32_T c5_i29;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), c5_dv3, 1, 0, 0U, 1, 0U, 2, 4, 2);
  for (c5_i29 = 0; c5_i29 < 8; c5_i29++) {
    c5_y[c5_i29] = c5_dv3[c5_i29];
  }

  sf_mex_destroy(&c5_u);
}

static void c5_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_corners;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  real_T c5_y[8];
  int32_T c5_i30;
  int32_T c5_i31;
  int32_T c5_i32;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_corners = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_corners), &c5_thisId, c5_y);
  sf_mex_destroy(&c5_corners);
  c5_i30 = 0;
  for (c5_i31 = 0; c5_i31 < 2; c5_i31++) {
    for (c5_i32 = 0; c5_i32 < 4; c5_i32++) {
      (*(real_T (*)[8])c5_outData)[c5_i32 + c5_i30] = c5_y[c5_i32 + c5_i30];
    }

    c5_i30 += 4;
  }

  sf_mex_destroy(&c5_mxArrayInData);
}

const mxArray *sf_c5_test_get_eml_resolved_functions_info(void)
{
  const mxArray *c5_nameCaptureInfo;
  c5_ResolvedFunctionInfo c5_info[28];
  const mxArray *c5_m0 = NULL;
  int32_T c5_i33;
  c5_ResolvedFunctionInfo *c5_r0;
  c5_nameCaptureInfo = NULL;
  c5_nameCaptureInfo = NULL;
  c5_info_helper(c5_info);
  sf_mex_assign(&c5_m0, sf_mex_createstruct("nameCaptureInfo", 1, 28), FALSE);
  for (c5_i33 = 0; c5_i33 < 28; c5_i33++) {
    c5_r0 = &c5_info[c5_i33];
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c5_r0->context)), "context", "nameCaptureInfo",
                    c5_i33);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c5_r0->name)), "name", "nameCaptureInfo", c5_i33);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c5_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c5_i33);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", c5_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c5_r0->resolved)), "resolved", "nameCaptureInfo",
                    c5_i33);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c5_i33);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c5_i33);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c5_i33);
    sf_mex_addfield(c5_m0, sf_mex_create("nameCaptureInfo", &c5_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c5_i33);
  }

  sf_mex_assign(&c5_nameCaptureInfo, c5_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c5_nameCaptureInfo);
  return c5_nameCaptureInfo;
}

static void c5_info_helper(c5_ResolvedFunctionInfo c5_info[28])
{
  c5_info[0].context = "";
  c5_info[0].name = "sin";
  c5_info[0].dominantType = "double";
  c5_info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c5_info[0].fileTimeLo = 1343830386U;
  c5_info[0].fileTimeHi = 0U;
  c5_info[0].mFileTimeLo = 0U;
  c5_info[0].mFileTimeHi = 0U;
  c5_info[1].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c5_info[1].name = "eml_scalar_sin";
  c5_info[1].dominantType = "double";
  c5_info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  c5_info[1].fileTimeLo = 1286818736U;
  c5_info[1].fileTimeHi = 0U;
  c5_info[1].mFileTimeLo = 0U;
  c5_info[1].mFileTimeHi = 0U;
  c5_info[2].context = "";
  c5_info[2].name = "cos";
  c5_info[2].dominantType = "double";
  c5_info[2].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c5_info[2].fileTimeLo = 1343830372U;
  c5_info[2].fileTimeHi = 0U;
  c5_info[2].mFileTimeLo = 0U;
  c5_info[2].mFileTimeHi = 0U;
  c5_info[3].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c5_info[3].name = "eml_scalar_cos";
  c5_info[3].dominantType = "double";
  c5_info[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  c5_info[3].fileTimeLo = 1286818722U;
  c5_info[3].fileTimeHi = 0U;
  c5_info[3].mFileTimeLo = 0U;
  c5_info[3].mFileTimeHi = 0U;
  c5_info[4].context = "";
  c5_info[4].name = "atan2";
  c5_info[4].dominantType = "double";
  c5_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m";
  c5_info[4].fileTimeLo = 1343830372U;
  c5_info[4].fileTimeHi = 0U;
  c5_info[4].mFileTimeLo = 0U;
  c5_info[4].mFileTimeHi = 0U;
  c5_info[5].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m";
  c5_info[5].name = "eml_scalar_eg";
  c5_info[5].dominantType = "double";
  c5_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c5_info[5].fileTimeLo = 1286818796U;
  c5_info[5].fileTimeHi = 0U;
  c5_info[5].mFileTimeLo = 0U;
  c5_info[5].mFileTimeHi = 0U;
  c5_info[6].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m";
  c5_info[6].name = "eml_scalexp_alloc";
  c5_info[6].dominantType = "double";
  c5_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c5_info[6].fileTimeLo = 1352424860U;
  c5_info[6].fileTimeHi = 0U;
  c5_info[6].mFileTimeLo = 0U;
  c5_info[6].mFileTimeHi = 0U;
  c5_info[7].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan2.m";
  c5_info[7].name = "eml_scalar_atan2";
  c5_info[7].dominantType = "double";
  c5_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan2.m";
  c5_info[7].fileTimeLo = 1286818720U;
  c5_info[7].fileTimeHi = 0U;
  c5_info[7].mFileTimeLo = 0U;
  c5_info[7].mFileTimeHi = 0U;
  c5_info[8].context = "";
  c5_info[8].name = "mpower";
  c5_info[8].dominantType = "double";
  c5_info[8].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c5_info[8].fileTimeLo = 1286818842U;
  c5_info[8].fileTimeHi = 0U;
  c5_info[8].mFileTimeLo = 0U;
  c5_info[8].mFileTimeHi = 0U;
  c5_info[9].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  c5_info[9].name = "power";
  c5_info[9].dominantType = "double";
  c5_info[9].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  c5_info[9].fileTimeLo = 1348191930U;
  c5_info[9].fileTimeHi = 0U;
  c5_info[9].mFileTimeLo = 0U;
  c5_info[9].mFileTimeHi = 0U;
  c5_info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  c5_info[10].name = "eml_scalar_eg";
  c5_info[10].dominantType = "double";
  c5_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c5_info[10].fileTimeLo = 1286818796U;
  c5_info[10].fileTimeHi = 0U;
  c5_info[10].mFileTimeLo = 0U;
  c5_info[10].mFileTimeHi = 0U;
  c5_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  c5_info[11].name = "eml_scalexp_alloc";
  c5_info[11].dominantType = "double";
  c5_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  c5_info[11].fileTimeLo = 1352424860U;
  c5_info[11].fileTimeHi = 0U;
  c5_info[11].mFileTimeLo = 0U;
  c5_info[11].mFileTimeHi = 0U;
  c5_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  c5_info[12].name = "floor";
  c5_info[12].dominantType = "double";
  c5_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c5_info[12].fileTimeLo = 1343830380U;
  c5_info[12].fileTimeHi = 0U;
  c5_info[12].mFileTimeLo = 0U;
  c5_info[12].mFileTimeHi = 0U;
  c5_info[13].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  c5_info[13].name = "eml_scalar_floor";
  c5_info[13].dominantType = "double";
  c5_info[13].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  c5_info[13].fileTimeLo = 1286818726U;
  c5_info[13].fileTimeHi = 0U;
  c5_info[13].mFileTimeLo = 0U;
  c5_info[13].mFileTimeHi = 0U;
  c5_info[14].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  c5_info[14].name = "eml_scalar_eg";
  c5_info[14].dominantType = "double";
  c5_info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c5_info[14].fileTimeLo = 1286818796U;
  c5_info[14].fileTimeHi = 0U;
  c5_info[14].mFileTimeLo = 0U;
  c5_info[14].mFileTimeHi = 0U;
  c5_info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  c5_info[15].name = "mtimes";
  c5_info[15].dominantType = "double";
  c5_info[15].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c5_info[15].fileTimeLo = 1289519692U;
  c5_info[15].fileTimeHi = 0U;
  c5_info[15].mFileTimeLo = 0U;
  c5_info[15].mFileTimeHi = 0U;
  c5_info[16].context = "";
  c5_info[16].name = "sqrt";
  c5_info[16].dominantType = "double";
  c5_info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c5_info[16].fileTimeLo = 1343830386U;
  c5_info[16].fileTimeHi = 0U;
  c5_info[16].mFileTimeLo = 0U;
  c5_info[16].mFileTimeHi = 0U;
  c5_info[17].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c5_info[17].name = "eml_error";
  c5_info[17].dominantType = "char";
  c5_info[17].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  c5_info[17].fileTimeLo = 1343830358U;
  c5_info[17].fileTimeHi = 0U;
  c5_info[17].mFileTimeLo = 0U;
  c5_info[17].mFileTimeHi = 0U;
  c5_info[18].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m";
  c5_info[18].name = "eml_scalar_sqrt";
  c5_info[18].dominantType = "double";
  c5_info[18].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m";
  c5_info[18].fileTimeLo = 1286818738U;
  c5_info[18].fileTimeHi = 0U;
  c5_info[18].mFileTimeLo = 0U;
  c5_info[18].mFileTimeHi = 0U;
  c5_info[19].context = "";
  c5_info[19].name = "mtimes";
  c5_info[19].dominantType = "double";
  c5_info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c5_info[19].fileTimeLo = 1289519692U;
  c5_info[19].fileTimeHi = 0U;
  c5_info[19].mFileTimeLo = 0U;
  c5_info[19].mFileTimeHi = 0U;
  c5_info[20].context = "";
  c5_info[20].name = "atan";
  c5_info[20].dominantType = "double";
  c5_info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan.m";
  c5_info[20].fileTimeLo = 1343830372U;
  c5_info[20].fileTimeHi = 0U;
  c5_info[20].mFileTimeLo = 0U;
  c5_info[20].mFileTimeHi = 0U;
  c5_info[21].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/atan.m";
  c5_info[21].name = "eml_scalar_atan";
  c5_info[21].dominantType = "double";
  c5_info[21].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_atan.m";
  c5_info[21].fileTimeLo = 1286818718U;
  c5_info[21].fileTimeHi = 0U;
  c5_info[21].mFileTimeLo = 0U;
  c5_info[21].mFileTimeHi = 0U;
  c5_info[22].context = "";
  c5_info[22].name = "abs";
  c5_info[22].dominantType = "double";
  c5_info[22].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c5_info[22].fileTimeLo = 1343830366U;
  c5_info[22].fileTimeHi = 0U;
  c5_info[22].mFileTimeLo = 0U;
  c5_info[22].mFileTimeHi = 0U;
  c5_info[23].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  c5_info[23].name = "eml_scalar_abs";
  c5_info[23].dominantType = "double";
  c5_info[23].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  c5_info[23].fileTimeLo = 1286818712U;
  c5_info[23].fileTimeHi = 0U;
  c5_info[23].mFileTimeLo = 0U;
  c5_info[23].mFileTimeHi = 0U;
  c5_info[24].context = "";
  c5_info[24].name = "mrdivide";
  c5_info[24].dominantType = "double";
  c5_info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c5_info[24].fileTimeLo = 1357951548U;
  c5_info[24].fileTimeHi = 0U;
  c5_info[24].mFileTimeLo = 1319729966U;
  c5_info[24].mFileTimeHi = 0U;
  c5_info[25].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  c5_info[25].name = "rdivide";
  c5_info[25].dominantType = "double";
  c5_info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c5_info[25].fileTimeLo = 1346510388U;
  c5_info[25].fileTimeHi = 0U;
  c5_info[25].mFileTimeLo = 0U;
  c5_info[25].mFileTimeHi = 0U;
  c5_info[26].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c5_info[26].name = "eml_scalexp_compatible";
  c5_info[26].dominantType = "double";
  c5_info[26].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  c5_info[26].fileTimeLo = 1286818796U;
  c5_info[26].fileTimeHi = 0U;
  c5_info[26].mFileTimeLo = 0U;
  c5_info[26].mFileTimeHi = 0U;
  c5_info[27].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  c5_info[27].name = "eml_div";
  c5_info[27].dominantType = "double";
  c5_info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  c5_info[27].fileTimeLo = 1313347810U;
  c5_info[27].fileTimeHi = 0U;
  c5_info[27].mFileTimeLo = 0U;
  c5_info[27].mFileTimeHi = 0U;
}

static real_T c5_atan2(SFc5_testInstanceStruct *chartInstance, real_T c5_y,
  real_T c5_x)
{
  real_T c5_b_y;
  real_T c5_b_x;
  c5_eml_scalar_eg(chartInstance);
  c5_b_y = c5_y;
  c5_b_x = c5_x;
  return muDoubleScalarAtan2(c5_b_y, c5_b_x);
}

static void c5_eml_scalar_eg(SFc5_testInstanceStruct *chartInstance)
{
}

static real_T c5_mpower(SFc5_testInstanceStruct *chartInstance, real_T c5_a)
{
  real_T c5_b_a;
  real_T c5_c_a;
  real_T c5_ak;
  real_T c5_d_a;
  real_T c5_e_a;
  real_T c5_b;
  c5_b_a = c5_a;
  c5_c_a = c5_b_a;
  c5_eml_scalar_eg(chartInstance);
  c5_ak = c5_c_a;
  c5_d_a = c5_ak;
  c5_eml_scalar_eg(chartInstance);
  c5_e_a = c5_d_a;
  c5_b = c5_d_a;
  return c5_e_a * c5_b;
}

static real_T c5_sqrt(SFc5_testInstanceStruct *chartInstance, real_T c5_x)
{
  real_T c5_b_x;
  c5_b_x = c5_x;
  c5_b_sqrt(chartInstance, &c5_b_x);
  return c5_b_x;
}

static void c5_eml_error(SFc5_testInstanceStruct *chartInstance)
{
  int32_T c5_i34;
  static char_T c5_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c5_u[30];
  const mxArray *c5_y = NULL;
  int32_T c5_i35;
  static char_T c5_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c5_b_u[4];
  const mxArray *c5_b_y = NULL;
  for (c5_i34 = 0; c5_i34 < 30; c5_i34++) {
    c5_u[c5_i34] = c5_cv0[c5_i34];
  }

  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", c5_u, 10, 0U, 1U, 0U, 2, 1, 30), FALSE);
  for (c5_i35 = 0; c5_i35 < 4; c5_i35++) {
    c5_b_u[c5_i35] = c5_cv1[c5_i35];
  }

  c5_b_y = NULL;
  sf_mex_assign(&c5_b_y, sf_mex_create("y", c5_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c5_y, 14, c5_b_y));
}

static const mxArray *c5_g_sf_marshallOut(void *chartInstanceVoid, void
  *c5_inData)
{
  const mxArray *c5_mxArrayOutData = NULL;
  int32_T c5_u;
  const mxArray *c5_y = NULL;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_mxArrayOutData = NULL;
  c5_u = *(int32_T *)c5_inData;
  c5_y = NULL;
  sf_mex_assign(&c5_y, sf_mex_create("y", &c5_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c5_mxArrayOutData, c5_y, FALSE);
  return c5_mxArrayOutData;
}

static int32_T c5_h_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  int32_T c5_y;
  int32_T c5_i36;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_i36, 1, 6, 0U, 0, 0U, 0);
  c5_y = c5_i36;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c5_mxArrayInData, const char_T *c5_varName, void *c5_outData)
{
  const mxArray *c5_b_sfEvent;
  const char_T *c5_identifier;
  emlrtMsgIdentifier c5_thisId;
  int32_T c5_y;
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)chartInstanceVoid;
  c5_b_sfEvent = sf_mex_dup(c5_mxArrayInData);
  c5_identifier = c5_varName;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_sfEvent),
    &c5_thisId);
  sf_mex_destroy(&c5_b_sfEvent);
  *(int32_T *)c5_outData = c5_y;
  sf_mex_destroy(&c5_mxArrayInData);
}

static uint8_T c5_i_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_b_is_active_c5_test, const char_T *c5_identifier)
{
  uint8_T c5_y;
  emlrtMsgIdentifier c5_thisId;
  c5_thisId.fIdentifier = c5_identifier;
  c5_thisId.fParent = NULL;
  c5_y = c5_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c5_b_is_active_c5_test),
    &c5_thisId);
  sf_mex_destroy(&c5_b_is_active_c5_test);
  return c5_y;
}

static uint8_T c5_j_emlrt_marshallIn(SFc5_testInstanceStruct *chartInstance,
  const mxArray *c5_u, const emlrtMsgIdentifier *c5_parentId)
{
  uint8_T c5_y;
  uint8_T c5_u0;
  sf_mex_import(c5_parentId, sf_mex_dup(c5_u), &c5_u0, 1, 3, 0U, 0, 0U, 0);
  c5_y = c5_u0;
  sf_mex_destroy(&c5_u);
  return c5_y;
}

static void c5_b_sqrt(SFc5_testInstanceStruct *chartInstance, real_T *c5_x)
{
  if (*c5_x < 0.0) {
    c5_eml_error(chartInstance);
  }

  *c5_x = muDoubleScalarSqrt(*c5_x);
}

static void init_dsm_address_info(SFc5_testInstanceStruct *chartInstance)
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

void sf_c5_test_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3954162564U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1495588255U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1226055220U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2522345534U);
}

mxArray *sf_c5_test_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("qP4YQTnYzy2Tqd6p9df6aC");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

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
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(1));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(4);
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
      pr[0] = (double)(4);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c5_test_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c5_test(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x6'type','srcId','name','auxInfo'{{M[1],M[10],T\"cornersX\",},{M[1],M[11],T\"cornersY\",},{M[1],M[8],T\"inFrame\",},{M[1],M[9],T\"measuredNodeHeading\",},{M[1],M[5],T\"measuredNodePos\",},{M[8],M[0],T\"is_active_c5_test\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 6, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c5_test_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc5_testInstanceStruct *chartInstance;
    chartInstance = (SFc5_testInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _testMachineNumber_,
           5,
           1,
           1,
           7,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"nodePos");
          _SFD_SET_DATA_PROPS(1,2,0,1,"measuredNodePos");
          _SFD_SET_DATA_PROPS(2,2,0,1,"measuredNodeHeading");
          _SFD_SET_DATA_PROPS(3,1,1,0,"eta");
          _SFD_SET_DATA_PROPS(4,2,0,1,"inFrame");
          _SFD_SET_DATA_PROPS(5,2,0,1,"cornersX");
          _SFD_SET_DATA_PROPS(6,2,0,1,"cornersY");
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
        _SFD_CV_INIT_EML(0,1,1,1,0,0,1,2,0,2,1);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,2136);
        _SFD_CV_INIT_EML_IF(0,1,0,1790,1970,-1,2004);
        _SFD_CV_INIT_EML_FOR(0,1,0,528,540,1741);
        _SFD_CV_INIT_EML_FOR(0,1,1,1774,1786,2019);

        {
          static int caseStart[] = { -1, 570, 635, 702, 771 };

          static int caseExprEnd[] = { 8, 576, 641, 708, 777 };

          _SFD_CV_INIT_EML_SWITCH(0,1,0,548,557,837,5,&(caseStart[0]),
            &(caseExprEnd[0]));
        }

        {
          static int condStart[] = { 1794, 1857 };

          static int condEnd[] = { 1852, 1968 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,0,1794,1969,2,0,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

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
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_e_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_c_sf_marshallOut,(MexInFcnForType)c5_c_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_c_sf_marshallOut,(MexInFcnForType)c5_c_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 6;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_UINT8,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c5_b_sf_marshallOut,(MexInFcnForType)c5_b_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)
            c5_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c5_sf_marshallOut,(MexInFcnForType)
            c5_sf_marshallIn);
        }

        {
          real_T *c5_measuredNodePos;
          real_T *c5_measuredNodeHeading;
          boolean_T *c5_inFrame;
          real_T (*c5_nodePos)[3];
          real_T (*c5_eta)[6];
          real_T (*c5_cornersX)[4];
          real_T (*c5_cornersY)[4];
          c5_cornersY = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 5);
          c5_cornersX = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
          c5_inFrame = (boolean_T *)ssGetOutputPortSignal(chartInstance->S, 3);
          c5_eta = (real_T (*)[6])ssGetInputPortSignal(chartInstance->S, 1);
          c5_measuredNodeHeading = (real_T *)ssGetOutputPortSignal
            (chartInstance->S, 2);
          c5_measuredNodePos = (real_T *)ssGetOutputPortSignal(chartInstance->S,
            1);
          c5_nodePos = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c5_nodePos);
          _SFD_SET_DATA_VALUE_PTR(1U, c5_measuredNodePos);
          _SFD_SET_DATA_VALUE_PTR(2U, c5_measuredNodeHeading);
          _SFD_SET_DATA_VALUE_PTR(3U, *c5_eta);
          _SFD_SET_DATA_VALUE_PTR(4U, c5_inFrame);
          _SFD_SET_DATA_VALUE_PTR(5U, *c5_cornersX);
          _SFD_SET_DATA_VALUE_PTR(6U, *c5_cornersY);
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
  return "8KKzpvk7szkAFPHojwK82C";
}

static void sf_opaque_initialize_c5_test(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc5_testInstanceStruct*) chartInstanceVar)->S,0);
  initialize_params_c5_test((SFc5_testInstanceStruct*) chartInstanceVar);
  initialize_c5_test((SFc5_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c5_test(void *chartInstanceVar)
{
  enable_c5_test((SFc5_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c5_test(void *chartInstanceVar)
{
  disable_c5_test((SFc5_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c5_test(void *chartInstanceVar)
{
  sf_c5_test((SFc5_testInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c5_test(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c5_test((SFc5_testInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c5_test();/* state var info */
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

extern void sf_internal_set_sim_state_c5_test(SimStruct* S, const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c5_test();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c5_test((SFc5_testInstanceStruct*)chartInfo->chartInstance,
                        mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c5_test(SimStruct* S)
{
  return sf_internal_get_sim_state_c5_test(S);
}

static void sf_opaque_set_sim_state_c5_test(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c5_test(S, st);
}

static void sf_opaque_terminate_c5_test(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc5_testInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_optimization_info();
    }

    finalize_c5_test((SFc5_testInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc5_test((SFc5_testInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c5_test(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c5_test((SFc5_testInstanceStruct*)(((ChartInfoStruct *)
      ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c5_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      5);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,5,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,5,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,5);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,5,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,5,5);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=5; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,5);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3384835291U));
  ssSetChecksum1(S,(1691968672U));
  ssSetChecksum2(S,(359669321U));
  ssSetChecksum3(S,(2985012226U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c5_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c5_test(SimStruct *S)
{
  SFc5_testInstanceStruct *chartInstance;
  chartInstance = (SFc5_testInstanceStruct *)utMalloc(sizeof
    (SFc5_testInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc5_testInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c5_test;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c5_test;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c5_test;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c5_test;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c5_test;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c5_test;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c5_test;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c5_test;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c5_test;
  chartInstance->chartInfo.mdlStart = mdlStart_c5_test;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c5_test;
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

void c5_test_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c5_test(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c5_test(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c5_test(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c5_test_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
