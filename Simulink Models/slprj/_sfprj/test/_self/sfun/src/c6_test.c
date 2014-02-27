/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_sfun.h"
#include "c6_test.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "test_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c6_debug_family_names[9] = { "nargin", "nargout",
  "measuredNodePos", "measuredNodeHeading", "inFrame", "pos_ref_in",
  "pos_ref_out", "heading_ref", "heading" };

/* Function Declarations */
static void initialize_c6_test(SFc6_testInstanceStruct *chartInstance);
static void initialize_params_c6_test(SFc6_testInstanceStruct *chartInstance);
static void enable_c6_test(SFc6_testInstanceStruct *chartInstance);
static void disable_c6_test(SFc6_testInstanceStruct *chartInstance);
static void c6_update_debugger_state_c6_test(SFc6_testInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c6_test(SFc6_testInstanceStruct
  *chartInstance);
static void set_sim_state_c6_test(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_st);
static void finalize_c6_test(SFc6_testInstanceStruct *chartInstance);
static void sf_c6_test(SFc6_testInstanceStruct *chartInstance);
static void initSimStructsc6_test(SFc6_testInstanceStruct *chartInstance);
static void registerMessagesc6_test(SFc6_testInstanceStruct *chartInstance);
static void init_script_number_translation(uint32_T c6_machineNumber, uint32_T
  c6_chartNumber);
static const mxArray *c6_sf_marshallOut(void *chartInstanceVoid, void *c6_inData);
static real_T c6_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_b_heading, const char_T *c6_identifier);
static real_T c6_b_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId);
static void c6_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static const mxArray *c6_b_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static real_T c6_c_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_heading_ref, const char_T *c6_identifier);
static real_T c6_d_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId);
static void c6_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static const mxArray *c6_c_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static void c6_e_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_pos_ref_out, const char_T *c6_identifier, real_T c6_y[3]);
static void c6_f_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[3]);
static void c6_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static const mxArray *c6_d_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static const mxArray *c6_e_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData);
static int32_T c6_g_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId);
static void c6_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData);
static uint8_T c6_h_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_b_is_active_c6_test, const char_T *c6_identifier);
static uint8_T c6_i_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId);
static void init_dsm_address_info(SFc6_testInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c6_test(SFc6_testInstanceStruct *chartInstance)
{
  chartInstance->c6_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c6_heading_not_empty = FALSE;
  chartInstance->c6_is_active_c6_test = 0U;
}

static void initialize_params_c6_test(SFc6_testInstanceStruct *chartInstance)
{
}

static void enable_c6_test(SFc6_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c6_test(SFc6_testInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c6_update_debugger_state_c6_test(SFc6_testInstanceStruct
  *chartInstance)
{
}

static const mxArray *get_sim_state_c6_test(SFc6_testInstanceStruct
  *chartInstance)
{
  const mxArray *c6_st;
  const mxArray *c6_y = NULL;
  real_T c6_hoistedGlobal;
  real_T c6_u;
  const mxArray *c6_b_y = NULL;
  int32_T c6_i0;
  real_T c6_b_u[3];
  const mxArray *c6_c_y = NULL;
  real_T c6_b_hoistedGlobal;
  real_T c6_c_u;
  const mxArray *c6_d_y = NULL;
  uint8_T c6_c_hoistedGlobal;
  uint8_T c6_d_u;
  const mxArray *c6_e_y = NULL;
  real_T *c6_heading_ref;
  real_T (*c6_pos_ref_out)[3];
  c6_heading_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c6_pos_ref_out = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c6_st = NULL;
  c6_st = NULL;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_createcellarray(4), FALSE);
  c6_hoistedGlobal = *c6_heading_ref;
  c6_u = c6_hoistedGlobal;
  c6_b_y = NULL;
  sf_mex_assign(&c6_b_y, sf_mex_create("y", &c6_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c6_y, 0, c6_b_y);
  for (c6_i0 = 0; c6_i0 < 3; c6_i0++) {
    c6_b_u[c6_i0] = (*c6_pos_ref_out)[c6_i0];
  }

  c6_c_y = NULL;
  sf_mex_assign(&c6_c_y, sf_mex_create("y", c6_b_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_setcell(c6_y, 1, c6_c_y);
  c6_b_hoistedGlobal = chartInstance->c6_heading;
  c6_c_u = c6_b_hoistedGlobal;
  c6_d_y = NULL;
  if (!chartInstance->c6_heading_not_empty) {
    sf_mex_assign(&c6_d_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c6_d_y, sf_mex_create("y", &c6_c_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_setcell(c6_y, 2, c6_d_y);
  c6_c_hoistedGlobal = chartInstance->c6_is_active_c6_test;
  c6_d_u = c6_c_hoistedGlobal;
  c6_e_y = NULL;
  sf_mex_assign(&c6_e_y, sf_mex_create("y", &c6_d_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c6_y, 3, c6_e_y);
  sf_mex_assign(&c6_st, c6_y, FALSE);
  return c6_st;
}

static void set_sim_state_c6_test(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_st)
{
  const mxArray *c6_u;
  real_T c6_dv0[3];
  int32_T c6_i1;
  real_T *c6_heading_ref;
  real_T (*c6_pos_ref_out)[3];
  c6_heading_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c6_pos_ref_out = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c6_doneDoubleBufferReInit = TRUE;
  c6_u = sf_mex_dup(c6_st);
  *c6_heading_ref = c6_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c6_u, 0)), "heading_ref");
  c6_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c6_u, 1)),
                        "pos_ref_out", c6_dv0);
  for (c6_i1 = 0; c6_i1 < 3; c6_i1++) {
    (*c6_pos_ref_out)[c6_i1] = c6_dv0[c6_i1];
  }

  chartInstance->c6_heading = c6_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c6_u, 2)), "heading");
  chartInstance->c6_is_active_c6_test = c6_h_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c6_u, 3)), "is_active_c6_test");
  sf_mex_destroy(&c6_u);
  c6_update_debugger_state_c6_test(chartInstance);
  sf_mex_destroy(&c6_st);
}

static void finalize_c6_test(SFc6_testInstanceStruct *chartInstance)
{
}

static void sf_c6_test(SFc6_testInstanceStruct *chartInstance)
{
  int32_T c6_i2;
  int32_T c6_i3;
  int32_T c6_i4;
  real_T c6_hoistedGlobal;
  boolean_T c6_b_hoistedGlobal;
  int32_T c6_i5;
  real_T c6_measuredNodePos[3];
  real_T c6_measuredNodeHeading;
  boolean_T c6_inFrame;
  int32_T c6_i6;
  real_T c6_pos_ref_in[3];
  uint32_T c6_debug_family_var_map[9];
  real_T c6_nargin = 4.0;
  real_T c6_nargout = 2.0;
  real_T c6_pos_ref_out[3];
  real_T c6_heading_ref;
  int32_T c6_i7;
  real_T *c6_b_measuredNodeHeading;
  boolean_T *c6_b_inFrame;
  real_T *c6_b_heading_ref;
  real_T (*c6_b_pos_ref_out)[3];
  real_T (*c6_b_pos_ref_in)[3];
  real_T (*c6_b_measuredNodePos)[3];
  c6_b_pos_ref_in = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 3);
  c6_b_heading_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c6_b_inFrame = (boolean_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c6_b_measuredNodeHeading = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c6_b_pos_ref_out = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c6_b_measuredNodePos = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 5U, chartInstance->c6_sfEvent);
  for (c6_i2 = 0; c6_i2 < 3; c6_i2++) {
    _SFD_DATA_RANGE_CHECK((*c6_b_measuredNodePos)[c6_i2], 0U);
  }

  for (c6_i3 = 0; c6_i3 < 3; c6_i3++) {
    _SFD_DATA_RANGE_CHECK((*c6_b_pos_ref_out)[c6_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c6_b_measuredNodeHeading, 2U);
  _SFD_DATA_RANGE_CHECK((real_T)*c6_b_inFrame, 3U);
  _SFD_DATA_RANGE_CHECK(*c6_b_heading_ref, 4U);
  for (c6_i4 = 0; c6_i4 < 3; c6_i4++) {
    _SFD_DATA_RANGE_CHECK((*c6_b_pos_ref_in)[c6_i4], 5U);
  }

  chartInstance->c6_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 5U, chartInstance->c6_sfEvent);
  c6_hoistedGlobal = *c6_b_measuredNodeHeading;
  c6_b_hoistedGlobal = *c6_b_inFrame;
  for (c6_i5 = 0; c6_i5 < 3; c6_i5++) {
    c6_measuredNodePos[c6_i5] = (*c6_b_measuredNodePos)[c6_i5];
  }

  c6_measuredNodeHeading = c6_hoistedGlobal;
  c6_inFrame = c6_b_hoistedGlobal;
  for (c6_i6 = 0; c6_i6 < 3; c6_i6++) {
    c6_pos_ref_in[c6_i6] = (*c6_b_pos_ref_in)[c6_i6];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 9U, 9U, c6_debug_family_names,
    c6_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_nargin, 0U, c6_b_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_nargout, 1U, c6_b_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c6_measuredNodePos, 2U, c6_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c6_measuredNodeHeading, 3U, c6_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c6_inFrame, 4U, c6_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c6_pos_ref_in, 5U, c6_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c6_pos_ref_out, 6U, c6_c_sf_marshallOut,
    c6_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c6_heading_ref, 7U, c6_b_sf_marshallOut,
    c6_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&chartInstance->c6_heading, 8U,
    c6_sf_marshallOut, c6_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 4);
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 6);
  if (CV_EML_IF(0, 1, 0, !chartInstance->c6_heading_not_empty)) {
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 7);
    chartInstance->c6_heading = 0.0;
    chartInstance->c6_heading_not_empty = TRUE;
  }

  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 9);
  c6_pos_ref_out[0] = 0.0;
  c6_pos_ref_out[1] = 0.0;
  c6_pos_ref_out[2] = c6_pos_ref_in[2];
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 10);
  if (CV_EML_IF(0, 1, 1, (real_T)c6_inFrame == 1.0)) {
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 11);
    c6_pos_ref_out[0] = c6_measuredNodePos[0];
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 12);
    c6_pos_ref_out[1] = c6_measuredNodePos[1];
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 13);
    c6_pos_ref_out[2] = c6_measuredNodePos[2];
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 14);
    chartInstance->c6_heading = c6_measuredNodeHeading;
  } else {
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 16);
    c6_pos_ref_out[0] = c6_pos_ref_in[0];
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 17);
    c6_pos_ref_out[1] = c6_pos_ref_in[1];
    _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 18);
    c6_pos_ref_out[2] = c6_pos_ref_in[2];
  }

  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, 21);
  c6_heading_ref = chartInstance->c6_heading;
  _SFD_EML_CALL(0U, chartInstance->c6_sfEvent, -21);
  _SFD_SYMBOL_SCOPE_POP();
  for (c6_i7 = 0; c6_i7 < 3; c6_i7++) {
    (*c6_b_pos_ref_out)[c6_i7] = c6_pos_ref_out[c6_i7];
  }

  *c6_b_heading_ref = c6_heading_ref;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 5U, chartInstance->c6_sfEvent);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_testMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc6_test(SFc6_testInstanceStruct *chartInstance)
{
}

static void registerMessagesc6_test(SFc6_testInstanceStruct *chartInstance)
{
}

static void init_script_number_translation(uint32_T c6_machineNumber, uint32_T
  c6_chartNumber)
{
}

static const mxArray *c6_sf_marshallOut(void *chartInstanceVoid, void *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  real_T c6_u;
  const mxArray *c6_y = NULL;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_u = *(real_T *)c6_inData;
  c6_y = NULL;
  if (!chartInstance->c6_heading_not_empty) {
    sf_mex_assign(&c6_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0), FALSE);
  } else {
    sf_mex_assign(&c6_y, sf_mex_create("y", &c6_u, 0, 0U, 0U, 0U, 0), FALSE);
  }

  sf_mex_assign(&c6_mxArrayOutData, c6_y, FALSE);
  return c6_mxArrayOutData;
}

static real_T c6_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_b_heading, const char_T *c6_identifier)
{
  real_T c6_y;
  emlrtMsgIdentifier c6_thisId;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_b_heading),
    &c6_thisId);
  sf_mex_destroy(&c6_b_heading);
  return c6_y;
}

static real_T c6_b_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId)
{
  real_T c6_y;
  real_T c6_d0;
  if (mxIsEmpty(c6_u)) {
    chartInstance->c6_heading_not_empty = FALSE;
  } else {
    chartInstance->c6_heading_not_empty = TRUE;
    sf_mex_import(c6_parentId, sf_mex_dup(c6_u), &c6_d0, 1, 0, 0U, 0, 0U, 0);
    c6_y = c6_d0;
  }

  sf_mex_destroy(&c6_u);
  return c6_y;
}

static void c6_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_b_heading;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_b_heading = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_b_heading),
    &c6_thisId);
  sf_mex_destroy(&c6_b_heading);
  *(real_T *)c6_outData = c6_y;
  sf_mex_destroy(&c6_mxArrayInData);
}

static const mxArray *c6_b_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  real_T c6_u;
  const mxArray *c6_y = NULL;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_u = *(real_T *)c6_inData;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", &c6_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, FALSE);
  return c6_mxArrayOutData;
}

static real_T c6_c_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_heading_ref, const char_T *c6_identifier)
{
  real_T c6_y;
  emlrtMsgIdentifier c6_thisId;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_heading_ref),
    &c6_thisId);
  sf_mex_destroy(&c6_heading_ref);
  return c6_y;
}

static real_T c6_d_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId)
{
  real_T c6_y;
  real_T c6_d1;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), &c6_d1, 1, 0, 0U, 0, 0U, 0);
  c6_y = c6_d1;
  sf_mex_destroy(&c6_u);
  return c6_y;
}

static void c6_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_heading_ref;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_heading_ref = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_heading_ref),
    &c6_thisId);
  sf_mex_destroy(&c6_heading_ref);
  *(real_T *)c6_outData = c6_y;
  sf_mex_destroy(&c6_mxArrayInData);
}

static const mxArray *c6_c_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_i8;
  real_T c6_b_inData[3];
  int32_T c6_i9;
  real_T c6_u[3];
  const mxArray *c6_y = NULL;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  for (c6_i8 = 0; c6_i8 < 3; c6_i8++) {
    c6_b_inData[c6_i8] = (*(real_T (*)[3])c6_inData)[c6_i8];
  }

  for (c6_i9 = 0; c6_i9 < 3; c6_i9++) {
    c6_u[c6_i9] = c6_b_inData[c6_i9];
  }

  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", c6_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, FALSE);
  return c6_mxArrayOutData;
}

static void c6_e_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_pos_ref_out, const char_T *c6_identifier, real_T c6_y[3])
{
  emlrtMsgIdentifier c6_thisId;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_pos_ref_out), &c6_thisId,
                        c6_y);
  sf_mex_destroy(&c6_pos_ref_out);
}

static void c6_f_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance, const
  mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId, real_T c6_y[3])
{
  real_T c6_dv1[3];
  int32_T c6_i10;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), c6_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c6_i10 = 0; c6_i10 < 3; c6_i10++) {
    c6_y[c6_i10] = c6_dv1[c6_i10];
  }

  sf_mex_destroy(&c6_u);
}

static void c6_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_pos_ref_out;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  real_T c6_y[3];
  int32_T c6_i11;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_pos_ref_out = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_pos_ref_out), &c6_thisId,
                        c6_y);
  sf_mex_destroy(&c6_pos_ref_out);
  for (c6_i11 = 0; c6_i11 < 3; c6_i11++) {
    (*(real_T (*)[3])c6_outData)[c6_i11] = c6_y[c6_i11];
  }

  sf_mex_destroy(&c6_mxArrayInData);
}

static const mxArray *c6_d_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  boolean_T c6_u;
  const mxArray *c6_y = NULL;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_u = *(boolean_T *)c6_inData;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", &c6_u, 11, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, FALSE);
  return c6_mxArrayOutData;
}

const mxArray *sf_c6_test_get_eml_resolved_functions_info(void)
{
  const mxArray *c6_nameCaptureInfo = NULL;
  c6_nameCaptureInfo = NULL;
  sf_mex_assign(&c6_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), FALSE);
  return c6_nameCaptureInfo;
}

static const mxArray *c6_e_sf_marshallOut(void *chartInstanceVoid, void
  *c6_inData)
{
  const mxArray *c6_mxArrayOutData = NULL;
  int32_T c6_u;
  const mxArray *c6_y = NULL;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_mxArrayOutData = NULL;
  c6_u = *(int32_T *)c6_inData;
  c6_y = NULL;
  sf_mex_assign(&c6_y, sf_mex_create("y", &c6_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c6_mxArrayOutData, c6_y, FALSE);
  return c6_mxArrayOutData;
}

static int32_T c6_g_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId)
{
  int32_T c6_y;
  int32_T c6_i12;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), &c6_i12, 1, 6, 0U, 0, 0U, 0);
  c6_y = c6_i12;
  sf_mex_destroy(&c6_u);
  return c6_y;
}

static void c6_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c6_mxArrayInData, const char_T *c6_varName, void *c6_outData)
{
  const mxArray *c6_b_sfEvent;
  const char_T *c6_identifier;
  emlrtMsgIdentifier c6_thisId;
  int32_T c6_y;
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)chartInstanceVoid;
  c6_b_sfEvent = sf_mex_dup(c6_mxArrayInData);
  c6_identifier = c6_varName;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_g_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_b_sfEvent),
    &c6_thisId);
  sf_mex_destroy(&c6_b_sfEvent);
  *(int32_T *)c6_outData = c6_y;
  sf_mex_destroy(&c6_mxArrayInData);
}

static uint8_T c6_h_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_b_is_active_c6_test, const char_T *c6_identifier)
{
  uint8_T c6_y;
  emlrtMsgIdentifier c6_thisId;
  c6_thisId.fIdentifier = c6_identifier;
  c6_thisId.fParent = NULL;
  c6_y = c6_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c6_b_is_active_c6_test),
    &c6_thisId);
  sf_mex_destroy(&c6_b_is_active_c6_test);
  return c6_y;
}

static uint8_T c6_i_emlrt_marshallIn(SFc6_testInstanceStruct *chartInstance,
  const mxArray *c6_u, const emlrtMsgIdentifier *c6_parentId)
{
  uint8_T c6_y;
  uint8_T c6_u0;
  sf_mex_import(c6_parentId, sf_mex_dup(c6_u), &c6_u0, 1, 3, 0U, 0, 0U, 0);
  c6_y = c6_u0;
  sf_mex_destroy(&c6_u);
  return c6_y;
}

static void init_dsm_address_info(SFc6_testInstanceStruct *chartInstance)
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

void sf_c6_test_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2468929871U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3540964573U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2027743532U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1140265104U);
}

mxArray *sf_c6_test_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("pE5acqqrmbe6DnVIyFmfyC");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

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
      pr[0] = (double)(3);
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

mxArray *sf_c6_test_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c6_test(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x4'type','srcId','name','auxInfo'{{M[1],M[8],T\"heading_ref\",},{M[1],M[5],T\"pos_ref_out\",},{M[4],M[0],T\"heading\",S'l','i','p'{{M1x2[122 129],M[0],}}},{M[8],M[0],T\"is_active_c6_test\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 4, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c6_test_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc6_testInstanceStruct *chartInstance;
    chartInstance = (SFc6_testInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _testMachineNumber_,
           6,
           1,
           1,
           6,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"measuredNodePos");
          _SFD_SET_DATA_PROPS(1,2,0,1,"pos_ref_out");
          _SFD_SET_DATA_PROPS(2,1,1,0,"measuredNodeHeading");
          _SFD_SET_DATA_PROPS(3,1,1,0,"inFrame");
          _SFD_SET_DATA_PROPS(4,2,0,1,"heading_ref");
          _SFD_SET_DATA_PROPS(5,1,1,0,"pos_ref_in");
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
        _SFD_CV_INIT_EML(0,1,1,2,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,528);
        _SFD_CV_INIT_EML_IF(0,1,0,132,151,-1,173);
        _SFD_CV_INIT_EML_IF(0,1,1,209,224,383,500);
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
            1.0,0,0,(MexFcnForType)c6_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c6_c_sf_marshallOut,(MexInFcnForType)
            c6_c_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c6_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_UINT8,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c6_d_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c6_b_sf_marshallOut,(MexInFcnForType)c6_b_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c6_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c6_measuredNodeHeading;
          boolean_T *c6_inFrame;
          real_T *c6_heading_ref;
          real_T (*c6_measuredNodePos)[3];
          real_T (*c6_pos_ref_out)[3];
          real_T (*c6_pos_ref_in)[3];
          c6_pos_ref_in = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S,
            3);
          c6_heading_ref = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c6_inFrame = (boolean_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c6_measuredNodeHeading = (real_T *)ssGetInputPortSignal
            (chartInstance->S, 1);
          c6_pos_ref_out = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S,
            1);
          c6_measuredNodePos = (real_T (*)[3])ssGetInputPortSignal
            (chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c6_measuredNodePos);
          _SFD_SET_DATA_VALUE_PTR(1U, *c6_pos_ref_out);
          _SFD_SET_DATA_VALUE_PTR(2U, c6_measuredNodeHeading);
          _SFD_SET_DATA_VALUE_PTR(3U, c6_inFrame);
          _SFD_SET_DATA_VALUE_PTR(4U, c6_heading_ref);
          _SFD_SET_DATA_VALUE_PTR(5U, *c6_pos_ref_in);
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
  return "WmWjeQtDkWPTWsxyaUEPvH";
}

static void sf_opaque_initialize_c6_test(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc6_testInstanceStruct*) chartInstanceVar)->S,0);
  initialize_params_c6_test((SFc6_testInstanceStruct*) chartInstanceVar);
  initialize_c6_test((SFc6_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c6_test(void *chartInstanceVar)
{
  enable_c6_test((SFc6_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c6_test(void *chartInstanceVar)
{
  disable_c6_test((SFc6_testInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c6_test(void *chartInstanceVar)
{
  sf_c6_test((SFc6_testInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c6_test(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c6_test((SFc6_testInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c6_test();/* state var info */
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

extern void sf_internal_set_sim_state_c6_test(SimStruct* S, const mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c6_test();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c6_test((SFc6_testInstanceStruct*)chartInfo->chartInstance,
                        mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c6_test(SimStruct* S)
{
  return sf_internal_get_sim_state_c6_test(S);
}

static void sf_opaque_set_sim_state_c6_test(SimStruct* S, const mxArray *st)
{
  sf_internal_set_sim_state_c6_test(S, st);
}

static void sf_opaque_terminate_c6_test(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc6_testInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_optimization_info();
    }

    finalize_c6_test((SFc6_testInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc6_test((SFc6_testInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c6_test(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c6_test((SFc6_testInstanceStruct*)(((ChartInfoStruct *)
      ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c6_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      6);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,6,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,6,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,6);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,6,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,6,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 4; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,6);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1820330611U));
  ssSetChecksum1(S,(3876188706U));
  ssSetChecksum2(S,(4026847322U));
  ssSetChecksum3(S,(593117327U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c6_test(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c6_test(SimStruct *S)
{
  SFc6_testInstanceStruct *chartInstance;
  chartInstance = (SFc6_testInstanceStruct *)utMalloc(sizeof
    (SFc6_testInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc6_testInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c6_test;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c6_test;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c6_test;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c6_test;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c6_test;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c6_test;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c6_test;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c6_test;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c6_test;
  chartInstance->chartInfo.mdlStart = mdlStart_c6_test;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c6_test;
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

void c6_test_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c6_test(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c6_test(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c6_test(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c6_test_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
