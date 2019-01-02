// -*- coding:utf-8-unix;mode:c++ -*-
// Copyright [2019] tadashi9e
#include <iostream>
#include "qc.h"

void epr(qc::qbit const& a,
         qc::qbit const& b) {
  hadamard(a);
  qc::cx(a, b);
}

void teleport(qc::qbit const& msg,
              qc::qbit const& here,
              qc::qbit const& there) {
  epr(here, there);
  qc::dump("EPR");
  qc::cx(msg, here);
  qc::dump("CNOT(msg, here)");
  qc::hadamard(msg);
  qc::dump("H(msg)");
  bool m_here = qc::measure(here);
  qc::dump(std::string("M(here)=") + (m_here ? "1" : "0"));
  if (m_here) {
    qc::pauli_x(there);
    qc::dump("X(there)");
  }
  bool m_msg = qc::measure(msg);
  qc::dump(std::string("M(msg)=") + (m_here ? "1" : "0"));
  if (m_msg) {
    qc::pauli_z(there);
    qc::dump("Z(there)");
  }
}

int main(int argc, char* argv[]) {
  // ---- 量子変数
  qc::qbit msg("msg");
  qc::qbit here("here");
  qc::qbit there("there");
  // ---- 初期化
  qc::reset();
  // qc::dump("reset");
  // qc::hadamard(msg);  qc::dump("hadamard(msg)");
  // qc::pauli_x(msg);  qc::dump("pauli_x(msg)");
  qc::dump("start teleportation");
  // ---- 量子テレポート実行
  teleport(msg, here, there);
  // ---- 量子テレポート結果の表示
  qc::dump("finish");
}
