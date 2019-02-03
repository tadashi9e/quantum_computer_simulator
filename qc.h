// -*- coding:utf-8-unix;mode:c++ -*-
// Copyright [2019] tadashi9e
#ifndef QC_H_
#define QC_H_

#include <complex>
#include <memory>
#include <string>

namespace qc {

typedef std::int16_t qbit_id_t;

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
  explicit qbit(std::string const& name);
  void set_name(std::string const& name);
  std::string const& get_name();
  qbit_id_t get_id() const;
  std::string str() const;
  struct q_impl;
 private:
  std::shared_ptr<q_impl> impl;
};

extern void dump(std::string const& title);
extern void dump(qbit_id_t n_bits);
extern void reset();

class frozen;
typedef std::shared_ptr<frozen> frozen_ptr;

extern frozen_ptr backup();
extern void restore(frozen_ptr const& frozenptr);

extern double measure(qbit_id_t id, bool is_up);
extern bool measure(qbit_id_t id);
extern void hadamard(qbit_id_t target_id);
extern void hadamard_up_to(qbit_id_t n_qbits);
extern void cphase(qbit_id_t target_id,
                   std::complex<double> const& phase);
extern void pauli_x(qbit_id_t target_id);
extern void cx(qbit_id_t control_id, qbit_id_t target_id);
extern void ccx(qbit_id_t control1_id, qbit_id_t control2_id,
                qbit_id_t target_id);
extern void pauli_y(qbit_id_t target_id);
extern void cy(qbit_id_t control_id,
               qbit_id_t target_id);
extern void ccy(qbit_id_t control_id1, qbit_id_t control_id2,
                qbit_id_t target_id);
extern void pauli_z(qbit_id_t target_id);
extern void cz(qbit_id_t control_id,
               qbit_id_t target_id);
extern void ccz(qbit_id_t control1_id, qbit_id_t control2_id,
                qbit_id_t target_id);
extern void cphase(qbit_id_t id, std::complex<double> const& phase);

extern double measure(qbit const& q, bool is_up);
extern bool measure(qbit const& q);
extern void hadamard(qbit const& q);
extern void hadamard_for_all();
extern void pauli_x(qbit const& q);
extern void cx(qbit const& control_q, qbit const& target_q);
extern void ccx(qbit const& control1_q, qbit const& control2_q,
                qbit const& target_q);
extern void pauli_y(qbit const& q);
extern void cy(qbit const& control_q, qbit const& target_q);
extern void ccy(qbit const& control1_q, qbit const& control2_q,
                qbit const& target_q);
extern void pauli_z(qbit const& q);
extern void cz(qbit const& control_q, qbit const& target_q);
extern void ccz(qbit const& control1_q, qbit const& control2_q,
                qbit const& target_q);
extern void cphase(qbit const& q, std::complex<double> const& phase);


}  // namespace qc

#endif  // QC_H_
