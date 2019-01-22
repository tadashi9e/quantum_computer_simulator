// -*- coding:utf-8-unix;mode:c++ -*-
// Copyright [2019] tadashi9e
#include "qc.h"
#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>
#include <random>

namespace qc {

class Random01 {
  std::random_device rnd_dev;
  std::mt19937 rnd_gen;
 public:
  Random01() : rnd_gen(rnd_dev()) {
  }
  double get_random() {
    std::uniform_real_distribution<> rnd01(0.0, 1.0);
    double rnd = rnd01(rnd_gen);
    return rnd;
  }
};

Random01 random01;
typedef std::uint64_t basis_t;
typedef std::unordered_map<basis_t, std::complex<double> > amplitudes_t;

/**
 * 右から順に出現順に量子変数を並べたときの
 * |000...000> から |111...111> までの確率振幅。
 */
static amplitudes_t q_amplitudes;
/**
 * 計算対象の量子変数
 */
static std::vector<qbit> qbits;

/**
 * 新たな量子変数に与える番号
 */
static qbit_id_t
new_qbit_id() {
  return static_cast<qbit_id_t>(qbits.size());
}

struct qbit::q_impl {
  qbit_id_t id;
  std::string name;
  q_impl() : id(-1), name() {
  }
  void setup() {
    id = new_qbit_id();
  }
  void set_name(std::string const& name) {
    this->name = name;
  }
  std::string const& get_name() {
    return name;
  }
};

qbit::qbit()
  : impl(std::make_shared<q_impl>()) {
  impl->setup();
  qbits.push_back(*this);
}

qbit::qbit(std::string const& name)
  : impl(std::make_shared<q_impl>()) {
  impl->setup();
  impl->set_name(name);
  qbits.push_back(*this);
}

void qbit::set_name(std::string const& name) {
  impl->set_name(name);
}

std::string const& qbit::get_name() {
  return impl->get_name();
}

/**
 * 量子変数の番号を返す。
 */
qbit_id_t qbit::get_id() const {
  return impl->id;
}

std::string qbit::str() const {
  return std::string("[") + std::to_string(impl->id)
    + "]" + impl->name;
}

static std::string bitset_of(basis_t n, basis_t max_n) {
  std::string result;
  basis_t pow_i = 1;
  for (qbit_id_t i = 0; ; ++i) {
    result += (n & (static_cast<basis_t>(1) << i)) ? '1' : '0';
    pow_i *= 2;
    if (pow_i >= max_n) {
      break;
    }
  }
  return result;
}

static double amplitude(std::complex<double> const& q_amp) {
  return q_amp.real() * q_amp.real() + q_amp.imag() * q_amp.imag();
}

/**
 * 量子変数の状態を表示する。
 */
  void dump(std::string const& title) {
  std::cout << "===== " << title << " =====" << std::endl;
  std::cout << "-- qbits --" << std::endl;
  for(qbit& qbit : qbits) {
    std::cout << qbit.str() << std::endl;
  }
  std::cout << "-- states --" << std::endl;
  std::vector<basis_t> bases;
  bases.reserve(q_amplitudes.size());
  for(amplitudes_t::value_type const& v : q_amplitudes) {
    bases.push_back(v.first);
  }
  std::sort(bases.begin(), bases.end());
  for(basis_t basis : bases) {
    std::complex<double> const& q_amp(q_amplitudes[basis]);
    std::cout << bitset_of(basis, static_cast<basis_t>(1) << qbits.size())
              << ": |"
              << q_amp
              << "| = " << amplitude(q_amp)
              << std::endl;
  }
}

/**
 * |000...000> の確率が 1, それ以外は 0 として
 * 量子変数を初期化する。
 */
void
reset() {
  q_amplitudes.clear();
  q_amplitudes[0] = std::complex<double>(1, 0);
}

class frozen {
 public:
  explicit frozen(amplitudes_t const& q_amp) : q_amp(q_amp) {
  }
  amplitudes_t const& get() const {
    return q_amp;
  }
 private:
  amplitudes_t q_amp;
};

frozen_ptr backup() {
  return std::make_shared<frozen>(q_amplitudes);
}
void restore(frozen_ptr const& frozenptr) {
  q_amplitudes = frozenptr->get();
}

// ----------------------------------------------------------------------

/**
 * /    \
 * | 1 0|
 * | 0 1|
 * \    /
 */
static std::complex<double> op_unit(qbit_id_t j, qbit_id_t i) {
  return (i == j)
    ? std::complex<double>(1, 0) : std::complex<double>();
}

/**
 * /     \
 * | 1  1| / sqrt(2)
 * | 1 -1|
 * \     /
 */
static std::complex<double> op_hadamard(qbit_id_t j, qbit_id_t i) {
  return ((i == 0 || j == 0)
          ? std::complex<double>(1, 0)
          : std::complex<double>(-1, 0)) / sqrt(2.0);
}

/**
 * /             \
 * | 1 0         |
 * | 0 exp(i phi)|
 * \             /
 */
static std::complex<double>
op_cphase(std::complex<double> const& phase, qbit_id_t j, qbit_id_t i) {
  if (j == 0 && i == 0) {
    return std::complex<double>(1, 0);
  }
  if (j == 1 && i == 1) {
    return phase;
  }
  return std::complex<double>();
}

/**
 * /    \
 * | 0 1|
 * | 1 0|
 * \    /
 */
static std::complex<double> op_pauli_x(qbit_id_t j, qbit_id_t i) {
  return (i != j) ? std::complex<double>(1, 0) : std::complex<double>();
}

/**
 * /     \
 * | 0 -i|
 * | i  0|
 * \     /
 */
static std::complex<double> op_pauli_y(qbit_id_t j, qbit_id_t i) {
  if (i == j) {
    return std::complex<double>();
  }
  if (i == 0) {
    return std::complex<double>(0, 1);
  } else {
    return std::complex<double>(0, -1);
  }
}

/**
 * /     \
 * | 1  0|
 * | 0 -1|
 * \     /
 */
static std::complex<double> op_pauli_z(qbit_id_t j, qbit_id_t i) {
  if (i != j) {
    return std::complex<double>();
  }
  if (i == 0) {
    return std::complex<double>(1, 0);
  } else {
    return std::complex<double>(-1, 0);
  }
}

static std::complex<double> op_downside(qbit_id_t j, qbit_id_t i) {
  return (i == 0 && j == 0)
    ? std::complex<double>(1, 0)
    : std::complex<double>();
}

static std::complex<double> op_upside(qbit_id_t j, qbit_id_t i) {
  return (i == 1 && j == 1)
    ? std::complex<double>(1, 0)
    : std::complex<double>();
}

/**
 * 二演算子のテンソル積
 */
static std::complex<double>
tensor_product(std::function<std::complex<double>(qbit_id_t, qbit_id_t)> op1,
               qbit_id_t j1, qbit_id_t i1,
               std::function<std::complex<double>(qbit_id_t, qbit_id_t)> op2,
               qbit_id_t j2, qbit_id_t i2) {
  return op1(j1, i1) * op2(j2, i2);
}

/**
 * 三演算子のテンソル積
 */
static std::complex<double>
tensor_product(std::function<std::complex<double>(qbit_id_t, qbit_id_t)> op1,
               qbit_id_t j1, qbit_id_t i1,
               std::function<std::complex<double>(qbit_id_t, qbit_id_t)> op2,
               qbit_id_t j2, qbit_id_t i2,
               std::function<std::complex<double>(qbit_id_t, qbit_id_t)> op3,
               qbit_id_t j3, qbit_id_t i3) {
  return (op1(j1, i1) *
          op2(j2, i2) *
          op3(j3, i3));
}

/**
 * CNOT ゲート
 */
static std::complex<double> op_cx(qbit_id_t control_j, qbit_id_t control_i,
                                  qbit_id_t target_j, qbit_id_t target_i) {
  return (tensor_product(&op_downside, control_j, control_i,
                         &op_unit,     target_j,  target_i) +
          tensor_product(&op_upside,   control_j, control_i,
                         &op_pauli_x,  target_j,  target_i));
}
static std::complex<double> op_cz(qbit_id_t control_j, qbit_id_t control_i,
                                  qbit_id_t target_j, qbit_id_t target_i) {
  return (tensor_product(&op_downside, control_j, control_i,
                         &op_unit,     target_j,   target_i) +
          tensor_product(&op_upside,   control_j,  control_i,
                         &op_pauli_z,  target_j,   target_i));
}

/**
 * Toffoli ゲート
 */
static std::complex<double> op_ccx(qbit_id_t control1_j, qbit_id_t control1_i,
                                   qbit_id_t control2_j, qbit_id_t control2_i,
                                   qbit_id_t target_j, qbit_id_t target_i) {
  return (tensor_product(&op_downside, control1_j, control1_i,
                         &op_downside, control2_j, control2_i,
                         &op_unit,     target_j,   target_i) +
          tensor_product(&op_upside,   control1_j, control1_i,
                         &op_downside, control2_j, control2_i,
                         &op_unit,     target_j,   target_i) +
          tensor_product(&op_downside, control1_j, control1_i,
                         &op_upside,   control2_j, control2_i,
                         &op_unit,     target_j,   target_i) +
          tensor_product(&op_upside,   control1_j, control1_i,
                         &op_upside,   control2_j, control2_i,
                         &op_pauli_x,  target_j,   target_i));
}
static std::complex<double> op_ccz(qbit_id_t control1_j, qbit_id_t control1_i,
                                   qbit_id_t control2_j, qbit_id_t control2_i,
                                   qbit_id_t target_j, qbit_id_t target_i) {
  return (tensor_product(&op_downside, control1_j, control1_i,
                         &op_downside, control2_j, control2_i,
                         &op_unit,     target_j,   target_i) +
          tensor_product(&op_upside,   control1_j, control1_i,
                         &op_downside, control2_j, control2_i,
                         &op_unit,     target_j,   target_i) +
          tensor_product(&op_downside, control1_j, control1_i,
                         &op_upside,   control2_j, control2_i,
                         &op_unit,     target_j,   target_i) +
          tensor_product(&op_upside,   control1_j, control1_i,
                         &op_upside,   control2_j, control2_i,
                         &op_pauli_z,  target_j,   target_i));
}

static amplitudes_t
apply_tensor_product(std::function<std::complex<double>(qbit_id_t, qbit_id_t)> op,
                     qbit_id_t id,
                     amplitudes_t const& amplitudes) {
  amplitudes_t amplitudes2;
  amplitudes2.reserve(amplitudes.size());
  basis_t imask(~(static_cast<basis_t>(0x01) << id));
  for(amplitudes_t::value_type const& v : amplitudes) {
    basis_t const basis(v.first);
    std::complex<double> const& amplitude(v.second);
    qbit_id_t i((basis >> id) & 0x01);
    for (qbit_id_t j = 0; j < 2; ++j) {
      basis_t const target_id((basis & imask) | (static_cast<basis_t>(j) << id));
      std::complex<double> w(op(j, i));
      if (w == std::complex<double>()) {
        continue;
      }
      amplitudes2[target_id] += (w * amplitude);
    }
  }
  return amplitudes2;
}

static amplitudes_t
apply_tensor_product(std::function<std::complex<double>(qbit_id_t, qbit_id_t,
                                                        qbit_id_t, qbit_id_t)> op2,
                     qbit_id_t id1, qbit_id_t id2,
                     amplitudes_t const& amplitudes) {
  amplitudes_t amplitudes2;
  amplitudes2.reserve(amplitudes.size());
  basis_t imask(~((static_cast<basis_t>(0x01) << id1) |
                  (static_cast<basis_t>(0x01) << id2)));
  for(amplitudes_t::value_type const& v : amplitudes) {
    basis_t const basis(v.first);
    std::complex<double> const& amplitude(v.second);
    qbit_id_t const i1((basis >> id1) & 0x01);
    qbit_id_t const i2((basis >> id2) & 0x01);
    for (qbit_id_t j2 = 0; j2 < 2; ++j2) {
      for (qbit_id_t j1 = 0; j1 < 2; ++j1) {
        std::complex<double> w(op2(j1, i1, j2, i2));
        if (w == std::complex<double>()) {
          continue;
        }
        basis_t const target_id((basis & imask) |
                                (static_cast<basis_t>(j1) << id1) |
                                (static_cast<basis_t>(j2) << id2));
        amplitudes2[target_id] += w * amplitude;
      }
    }
  }
  return amplitudes2;
}

static amplitudes_t
apply_tensor_product(std::function<std::complex<double>(qbit_id_t, qbit_id_t,
                                                        qbit_id_t, qbit_id_t,
                                                        qbit_id_t, qbit_id_t)> op3,
                     qbit_id_t id1, qbit_id_t id2, qbit_id_t id3,
                     amplitudes_t const& amplitudes) {
  amplitudes_t amplitudes2;
  amplitudes2.reserve(amplitudes.size());
  basis_t imask(~((static_cast<basis_t>(0x01) << id1) |
                  (static_cast<basis_t>(0x01) << id2) |
                  (static_cast<basis_t>(0x01) << id3)));
  for(amplitudes_t::value_type const& v : amplitudes) {
    basis_t basis(v.first);
    std::complex<double> const& amplitude(v.second);
    qbit_id_t const i1((basis >> id1) & 0x01);
    qbit_id_t const i2((basis >> id2) & 0x01);
    qbit_id_t const i3((basis >> id3) & 0x01);
    for (qbit_id_t j3 = 0; j3 < 2; ++j3) {
      for (qbit_id_t j2 = 0; j2 < 2; ++j2) {
        for (qbit_id_t j1 = 0; j1 < 2; ++j1) {
          std::complex<double> w(op3(j1, i1, j2, i2, j3, i3));
          if (w == std::complex<double>()) {
            continue;
          }
          basis_t const target_id((basis & imask) |
                                  (static_cast<basis_t>(j1) << id1) |
                                  (static_cast<basis_t>(j2) << id2) |
                                  (static_cast<basis_t>(j3) << id3));
          amplitudes2[target_id] += w * amplitude;
        }
      }
    }
  }
  return amplitudes2;
}

static double
op_measure(qbit_id_t id, bool is_up) {
  double p = 0;
  basis_t mask = static_cast<basis_t>(0x01) << id;
  for(amplitudes_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    std::complex<double> const& q_amp(v.second);
    bool b = ((basis & mask) == 0) ? false : true;
    if (b == is_up) {
      p += amplitude(q_amp);
    }
  }
  if (p != 0) {
    double w = 1.0 / sqrt(p);
    for (amplitudes_t::iterator iter = q_amplitudes.begin();
        iter != q_amplitudes.end();) {
      basis_t const basis(iter->first);
      std::complex<double>& q_amp(iter->second);
      bool b = ((basis & mask) == 0) ? false : true;
      if (b != is_up) {
        iter = q_amplitudes.erase(iter);
      } else {
        q_amp *= w;
        ++iter;
      }
    }
  }
  return p;
}

static bool
op_measure(qbit_id_t id) {
  double p0 = 0;
  double p1 = 0;
  basis_t mask = static_cast<basis_t>(0x01) << id;
  for(amplitudes_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    std::complex<double> const& q_amp(v.second);
    bool b = ((basis & mask) == 0) ? false : true;
    if (b) {
      p1 += amplitude(q_amp);
    } else {
      p0 += amplitude(q_amp);
    }
  }
  bool is_up = (p1 > random01.get_random()) ? true : false;
  double p = is_up ? p1 : p0;
  double w = 1.0 / sqrt(p);
  for (amplitudes_t::iterator iter = q_amplitudes.begin();
      iter != q_amplitudes.end();) {
    basis_t const basis(iter->first);
    std::complex<double>& q_amp(iter->second);
    bool b = ((basis & mask) == 0) ? false : true;
    if (b != is_up) {
      iter = q_amplitudes.erase(iter);
    } else {
      q_amp *= w;
      ++iter;
    }
  }
  return is_up;
}

double
measure(qbit const& q, bool is_up) {
  qbit_id_t const id(q.get_id());
  return op_measure(id, is_up);
}

bool
measure(qbit const& q) {
  qbit_id_t const id(q.get_id());
  return op_measure(id);
}

void
hadamard(qbit const& q) {
  qbit_id_t const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_hadamard, id, q_amplitudes);
}
void
hadamard_for_all() {
  if (q_amplitudes.size() == 1 &&
      q_amplitudes[0] == std::complex<double>(1, 0)) {
    basis_t const total_basis(static_cast<basis_t>(0x01) << qbits.size());
    double const average = 1.0 / sqrt(static_cast<double>(total_basis));
    q_amplitudes.reserve(total_basis);
    for (basis_t basis = 0; basis < total_basis; ++basis) {
      q_amplitudes[basis] = average;
    }
  } else {
    for (qbit_id_t id = 0; id < static_cast<qbit_id_t>(qbits.size()); ++id) {
      q_amplitudes = apply_tensor_product(op_hadamard, id, q_amplitudes);
    }
  }
}
void
cphase(qbit const& q, std::complex<double> const& phase) {
  qbit_id_t const id(q.get_id());
  std::function<std::complex<double>(qbit_id_t, qbit_id_t)>
    op(std::bind(&op_cphase, phase, std::placeholders::_1, std::placeholders::_2));
  q_amplitudes = apply_tensor_product(op, id, q_amplitudes);
}
void pauli_x(qbit const& q) {
  qbit_id_t const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_pauli_x, id, q_amplitudes);
}
void pauli_y(qbit const& q) {
  qbit_id_t const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_pauli_y, id, q_amplitudes);
}
void pauli_z(qbit const& q) {
  qbit_id_t const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_pauli_z, id, q_amplitudes);
}
void cx(qbit const& control_q, qbit const& target_q) {
  qbit_id_t const control_id(control_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  q_amplitudes = apply_tensor_product(op_cx, control_id, target_id,
                                      q_amplitudes);
}
void cz(qbit const& control_q, qbit const& target_q) {
  qbit_id_t const control_id(control_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  q_amplitudes = apply_tensor_product(op_cz, control_id, target_id,
                                      q_amplitudes);
}
void ccx(qbit const& control1_q,
         qbit const& control2_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control1_q.get_id());
  qbit_id_t const control2_id(control2_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  q_amplitudes = apply_tensor_product(op_ccx, control1_id, control2_id,
                                      target_id,
                                      q_amplitudes);
}
void ccz(qbit const& control1_q,
         qbit const& control2_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control1_q.get_id());
  qbit_id_t const control2_id(control2_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  q_amplitudes = apply_tensor_product(op_ccz, control1_id, control2_id,
                                      target_id,
                                      q_amplitudes);
}

}  // namespace qc
