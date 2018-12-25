// -*- coding:utf-8;mode:c++ -*-
#ifndef QC_H
#define QC_H

#include <complex>
#include <boost/shared_ptr.hpp>

namespace qc {

/**
 * qbit。
 * 
 * 出現順に ID 番号が振られる。
 * 宣言は qc::reset 呼び出し前でなければならない。
 * デバッグ表示の便宜のために名称を設定することができる。
 */
class qbit {
public:
  /**
   * qbit のコンストラクタ
   */
  qbit();
  /**
   * qbit のコンストラクタ（名称を与える）
   */
  qbit(std::string const& name);
  void set_name(std::string const& name);
  std::string const& get_name();
  int get_id() const;
  std::string str() const;
  struct q_impl;
private:
  boost::shared_ptr<q_impl> impl;
};

template<typename T>
T amplitude(std::complex<T> const& c) {
  return c.real() * c.real() + c.imag() * c.imag();
}

extern void dump(std::string const& title);
extern void reset();
extern double measure(qbit const& q, bool is_up);
extern bool measure(qbit const& q);
extern void hadamard(qbit const& q);
extern void dump_hadamard(qbit const& q);
extern void pauli_x(qbit const& q);
extern void dump_pauli_x(qbit const& q);
extern void pauli_y(qbit const& q);
extern void dump_pauli_y(qbit const& q);
extern void pauli_z(qbit const& q);
extern void dump_pauli_z(qbit const& q);
extern void cnot(qbit const& control_q, qbit const& target_q);
extern void dump_cnot(qbit const& control_q, qbit const& target_q);
void cphase(qbit const& q, std::complex<double> const& phase);


}

#endif
