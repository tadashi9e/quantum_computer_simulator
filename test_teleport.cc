#include "qc.h"
#include <iostream>

void epr(qc::qbit& a,
         qc::qbit& b) {
  hadamard(a);
  qc::cnot(a, b);
}

void teleport(qc::qbit& msg,
              qc::qbit& here,
              qc::qbit& there) {
  epr(here, there);
  qc::dump("EPR");
  qc::cnot(msg, here);
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
  //qc::dump("reset");
  //qc::hadamard(msg);
  //qc::pauli_x(msg);
  qc::dump("start teleportation");
  // ---- 量子テレポート実行
  teleport(msg, here, there);
  // ---- 量子テレポート結果の表示
  qc::dump("finish");
}
