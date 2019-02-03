# -*- coding:utf-8;mode:python -*-
# Copyright [2019] tadashi9e

u'''Interactive quantum computer simulator.
'''
import qc_py
import math

try:
    from typing import Optional
except:
    pass

is_interactive = True  # type: bool

def _interactive():
    # type: () -> None
    if is_interactive:
        dump()

def reset():
    # type: () -> None
    u'''Reset quantum computer status
    '''
    qc_py.qc_reset()
    _interactive()
def set_interactive(b):
    u'''Turn on/off interactive mode
    '''
    # type: (bool) -> None
    global is_interactive
    is_interactive = b
def x(target):
    # type: (int) -> None
    u'''Apply Pauli-X gate
    '''
    qc_py.qc_gate_x(target)
    _interactive()
def y(target):
    # type: (int) -> None
    u'''Apply Pauli-Y gate
    '''
    qc_py.qc_gate_y(target)
    _interactive()
def z(target):
    # type: (int) -> None
    u'''Apply Pauli-Z gate
    '''
    qc_py.qc_gate_z(target)
    _interactive()
def cx(control, target):
    # type: (int, int) -> None
    u'''Controlled Pauli-X gate
    '''
    qc_py.qc_gate_cx(control, target)
    _interactive()
def cy(control, target):
    # type: (int, int) -> None
    u'''Controlled Pauli-Y gate
    '''
    qc_py.qc_gate_ccy(control, target)
    _interactive()
def cz(control, target):
    # type: (int, int) -> None
    u'''Controlled Pauli-Z gate
    '''
    qc_py.qc_gate_cz(control, target)
    _interactive()
def ccx(control1, control2, target):
    # type: (int, int, int) -> None
    u'''Toffoli gate. Multi-controlled Pauli-X gate
    '''
    qc_py.qc_gate_ccx(control1, control2, target)
    _interactive()
def ccy(control1, control2, target):
    # type: (int, int, int) -> None
    u'''Multi-controlled Pauli-Y gate
    '''
    qc_py.qc_gate_ccy(control1, control2, target)
    _interactive()
def ccz(control1, control2, target):
    # type: (int, int, int) -> None
    u'''Multi-controlled Pauli-Z gate
    '''
    qc_py.qc_gate_ccz(control1, control2, target)
    _interactive()
def rz(target, phase):
    # type: (int, float) -> None
    u'''Phase shift gate
    '''
    real = math.cos(phase)
    imag = math.sin(phase)
    qc_py.qc_gate_rz(target, real, imag)
    _interactive()
def h(min_target, max_target=None):
    # type: (int, Optional[int]) -> None
    u'''Hadamard gate
    '''
    if not max_target:
        qc_py.qc_gate_h(min_target)
        _interactive()
        return
    if min_target == 0:
        qc_py.qc_gate_h_up_to(max_target)
    else:
        for target in range(min_target, max_target + 1):
            qc_py.qc_gate_h(target)
    _interactive()
def m(target):
    # type: (int) -> bool
    u'''Measure qubit
    '''
    b = qc_py.qc_gate_measure(target)
    _interactive()
    return b
def dump():
    # type: () -> None
    u'''Inspect quantum computer status
    '''
    qc_py.qc_dump()

reset()
