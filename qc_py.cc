// -*- coding:utf-8-unix;mode:c++ -*-
// Copyright [2019] tadashi9e
#include <iostream>
#include <Python.h>
#include "qc.h"

static qc::qbit_id_t n_qbits = 0;

static void
update_n_qbits(qc::qbit_id_t id) {
  n_qbits = std::max(n_qbits, id);
}

static PyObject*
qc_reset(PyObject* self, PyObject* args) {
  qc::reset();
  n_qbits = 0;
  return Py_None;
}
static PyObject*
qc_gate_x(PyObject* self, PyObject* args) {
  int target;
  if (!PyArg_ParseTuple(args, "i", &target)) {
    return NULL;
  }
  qc::pauli_x(target);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_y(PyObject* self, PyObject* args) {
  int target;
  if (!PyArg_ParseTuple(args, "i", &target)) {
    return NULL;
  }
  qc::pauli_y(target);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_z(PyObject* self, PyObject* args) {
  int target;
  if (!PyArg_ParseTuple(args, "i", &target)) {
    return NULL;
  }
  qc::pauli_z(target);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_cx(PyObject* self, PyObject* args) {
  int control;
  int target;
  if (!PyArg_ParseTuple(args, "ii", &control, &target)) {
    return NULL;
  }
  qc::cx(control, target);
  update_n_qbits(control);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_cy(PyObject* self, PyObject* args) {
  int control;
  int target;
  if (!PyArg_ParseTuple(args, "ii", &control, &target)) {
    return NULL;
  }
  qc::cy(control, target);
  update_n_qbits(control);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_cz(PyObject* self, PyObject* args) {
  int control;
  int target;
  if (!PyArg_ParseTuple(args, "ii", &control, &target)) {
    return NULL;
  }
  qc::cz(control, target);
  update_n_qbits(control);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_ccx(PyObject* self, PyObject* args) {
  int control1;
  int control2;
  int target;
  if (!PyArg_ParseTuple(args, "iii", &control1, &control2, &target)) {
    return NULL;
  }
  qc::ccx(control1, control2, target);
  update_n_qbits(control1);
  update_n_qbits(control2);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_ccy(PyObject* self, PyObject* args) {
  int control1;
  int control2;
  int target;
  if (!PyArg_ParseTuple(args, "iii", &control1, &control2, &target)) {
    return NULL;
  }
  qc::ccy(control1, control2, target);
  update_n_qbits(control1);
  update_n_qbits(control2);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_ccz(PyObject* self, PyObject* args) {
  int control1;
  int control2;
  int target;
  if (!PyArg_ParseTuple(args, "iii", &control1, &control2, &target)) {
    return NULL;
  }
  qc::ccz(control1, control2, target);
  update_n_qbits(control1);
  update_n_qbits(control2);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_rz(PyObject* self, PyObject* args) {
  int target;
  double real;
  double imag;
  if (!PyArg_ParseTuple(args, "idd", &target, &real, &imag)) {
    return NULL;
  }
  std::complex<double> const rot(real, imag);
  qc::cphase(target, rot);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_h(PyObject* self, PyObject* args) {
  int target;
  if (!PyArg_ParseTuple(args, "i", &target)) {
    return NULL;
  }
  qc::hadamard(target);
  update_n_qbits(target);
  return Py_None;
}
static PyObject*
qc_gate_h_up_to(PyObject* self, PyObject* args) {
  int max_qbit_id;
  if (!PyArg_ParseTuple(args, "i", &max_qbit_id)) {
    return NULL;
  }
  qc::hadamard_up_to(max_qbit_id + 1);
  update_n_qbits(max_qbit_id);
  return Py_None;
}
static PyObject*
qc_gate_measure(PyObject* self, PyObject* args) {
  int target;
  if (!PyArg_ParseTuple(args, "i", &target)) {
    return NULL;
  }
  bool const b(qc::measure(target));
  update_n_qbits(target);
  return b ? Py_True : Py_False;
}
static PyObject*
qc_dump(PyObject* self, PyObject* args) {
  qc::dump(n_qbits+1);
  return Py_None;
}

static PyMethodDef QcMethods[] = {
  {"qc_reset", qc_reset, METH_VARARGS, "initialize"},
  {"qc_gate_x", qc_gate_x, METH_VARARGS, "pauli_x"},
  {"qc_gate_y", qc_gate_y, METH_VARARGS, "gate_y"},
  {"qc_gate_z", qc_gate_z, METH_VARARGS, "gate_z"},
  {"qc_gate_cx", qc_gate_cx, METH_VARARGS, "gate_cx"},
  {"qc_gate_cy", qc_gate_cy, METH_VARARGS, "gate_cy"},
  {"qc_gate_cz", qc_gate_cz, METH_VARARGS, "gate_cz"},
  {"qc_gate_ccx", qc_gate_ccx, METH_VARARGS, "gate_ccx"},
  {"qc_gate_ccy", qc_gate_ccy, METH_VARARGS, "gate_ccy"},
  {"qc_gate_ccz", qc_gate_ccz, METH_VARARGS, "gate_ccz"},
  {"qc_gate_rz", qc_gate_rz, METH_VARARGS, "gate_rz"},
  {"qc_gate_h", qc_gate_h, METH_VARARGS, "gate_h"},
  {"qc_gate_h_up_to", qc_gate_h_up_to, METH_VARARGS, "gate_h_up_to"},
  {"qc_gate_measure", qc_gate_measure, METH_VARARGS, "gate_measure"},
  {"qc_dump",   qc_dump, METH_VARARGS, "dump"},
  {NULL, NULL, 0, NULL}  // Sentinel
};

static struct PyModuleDef qcpymodule = {
  PyModuleDef_HEAD_INIT, "qc_py", NULL, -1, QcMethods
};

PyMODINIT_FUNC
PyInit_qc_py(void) {
  return PyModule_Create(&qcpymodule);
}
