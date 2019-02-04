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
typedef std::int16_t qbit_id_t;
typedef std::uint64_t basis_t;
typedef std::complex<double> amplitude_t;
typedef std::unordered_map<basis_t, amplitude_t> amplitude_map_t;

static bool is_1(basis_t basis, qbit_id_t n) {
  return (basis & (static_cast<basis_t>(1) << n)) ? true : false;
}
static basis_t set_1(basis_t basis, qbit_id_t n) {
  return basis | (static_cast<basis_t>(1) << n);
}
static basis_t set_0(basis_t basis, qbit_id_t n) {
  return basis & ~(static_cast<basis_t>(1) << n);
}
static int count_1_bits(basis_t n) {
  // Humming weight
  const uint64_t m1  = 0x5555555555555555;  // binary: 0101...
  const uint64_t m2  = 0x3333333333333333;  // binary: 00110011..
  const uint64_t m4  = 0x0f0f0f0f0f0f0f0f;  // binary:  4 zeros,  4 ones ...
  const uint64_t m8  = 0x00ff00ff00ff00ff;  // binary:  8 zeros,  8 ones ...
  const uint64_t m16 = 0x0000ffff0000ffff;  // binary: 16 zeros, 16 ones ...
  const uint64_t m32 = 0x00000000ffffffff;  // binary: 32 zeros, 32 ones
  n = (n & m1) + ((n >>  1) & m1);
  n = (n & m2) + ((n >>  2) & m2);
  n = (n & m4) + ((n >>  4) & m4);
  n = (n & m8) + ((n >>  8) & m8);
  n = (n & m16) + ((n >> 16) & m16);
  n = (n & m32) + ((n >> 32) & m32);
  return n;
}
static void add_amplitude(amplitude_map_t& amplitudes,
                          basis_t basis,
                          amplitude_t amplitude) {
  amplitude_map_t::iterator
    iter = amplitudes.find(basis);
  if (iter == amplitudes.end()) {
    amplitudes[basis] = amplitude;
    return;
  }
  amplitude_t& a(iter->second);
  a += amplitude;
}
static std::string bitset_of(basis_t basis, int n_qbits) {
  std::string s;
  for (qbit_id_t i = 0; i < n_qbits; ++i) {
    s += is_1(basis, i) ? "1" : "0";
  }
  return s;
}

static double
probability(amplitude_t const& a) {
  return a.real() * a.real() + a.imag() * a.imag();
}

class QC {
 public:
  QC() {
    reset();
  }
  void reset() {
    amplitudes.clear();
    amplitudes[0] = amplitude_t(1.0, 0.0);
  }
  void gate_x(qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      basis_t const basis2(is_1(basis, n) ? set_0(basis, n) : set_1(basis, n));
      add_amplitude(amplitudes2, basis2, amplitude);
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_cx(qbit_id_t c, qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (!is_1(basis, c)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      basis_t const basis2(is_1(basis, n) ? set_0(basis, n) : set_1(basis, n));
      add_amplitude(amplitudes2, basis2, amplitude);
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_ccx(qbit_id_t c1, qbit_id_t c2, qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (!is_1(basis, c1)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      if (!is_1(basis, c2)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      basis_t const basis2(is_1(basis, n) ? set_0(basis, n) : set_1(basis, n));
      add_amplitude(amplitudes2, basis2, amplitude);
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_y(qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      basis_t const basis2(is_1(basis, n) ? set_0(basis, n) : set_1(basis, n));
      if (is_1(basis, n)) {
        add_amplitude(amplitudes2, basis2, amplitude * amplitude_t(0.0, -1.0));
      } else {
        add_amplitude(amplitudes2, basis2, amplitude * amplitude_t(0.0, 1.0));
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_cy(qbit_id_t c, qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (!is_1(basis, c)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      basis_t const basis2(is_1(basis, n) ? set_0(basis, n) : set_1(basis, n));
      if (is_1(basis, n)) {
        add_amplitude(amplitudes2, basis2, amplitude * amplitude_t(0.0, -1.0));
      } else {
        add_amplitude(amplitudes2, basis2, amplitude * amplitude_t(0.0, 1.0));
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_ccy(qbit_id_t c1, qbit_id_t c2, qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (!is_1(basis, c1)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      if (!is_1(basis, c2)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      basis_t const basis2(is_1(basis, n) ? set_0(basis, n) : set_1(basis, n));
      if (is_1(basis, n)) {
        add_amplitude(amplitudes2, basis2, amplitude * amplitude_t(0.0, -1.0));
      } else {
        add_amplitude(amplitudes2, basis2, amplitude * amplitude_t(0.0, 1.0));
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_z(qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (is_1(basis, n)) {
        add_amplitude(amplitudes2, basis, amplitude * amplitude_t(-1.0, 0.0));
      } else {
        add_amplitude(amplitudes2, basis, amplitude * amplitude_t(1.0, 0.0));
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_cz(qbit_id_t c, qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (!is_1(basis, c)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      if (is_1(basis, n)) {
        add_amplitude(amplitudes2, basis, amplitude * amplitude_t(-1.0, 0.0));
      } else {
        add_amplitude(amplitudes2, basis, amplitude * amplitude_t(1.0, 0.0));
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_ccz(qbit_id_t c1, qbit_id_t c2, qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (!is_1(basis, c1)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      if (!is_1(basis, c2)) {
        add_amplitude(amplitudes2, basis, amplitude);
        continue;
      }
      if (is_1(basis, n)) {
        add_amplitude(amplitudes2, basis, amplitude * amplitude_t(-1.0, 0.0));
      } else {
        add_amplitude(amplitudes2, basis, amplitude * amplitude_t(1.0, 0.0));
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  /**
   * /                                          \
   * | 1  1| / sqrt(2)
   * | 1 -1|
   * \     /
   */
  void gate_h(qbit_id_t n) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (is_1(basis, n)) {
        basis_t const basis10(set_0(basis, n));
        add_amplitude(amplitudes2, basis10, amplitude / sqrt(2.0));
        basis_t const basis11(basis);
        add_amplitude(amplitudes2, basis11, -amplitude / sqrt(2.0));
      } else {
        basis_t const basis00(basis);
        add_amplitude(amplitudes2, basis00, amplitude / sqrt(2.0));
        basis_t const basis01(set_1(basis, n));
        add_amplitude(amplitudes2, basis01, amplitude / sqrt(2.0));
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_h_up_to(qbit_id_t n_qbits) {
    basis_t const total_basis(static_cast<basis_t>(0x01) << n_qbits);
    amplitude_map_t amplitudes2;
    double const p(pow(sqrt(2.0), n_qbits));
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& q_amp(v.second);
      amplitude_t const q_amp2(q_amp / p);
      for (basis_t basis2 = 0; basis2 < total_basis; ++basis2) {
        if (count_1_bits(basis & basis2) & 0x01) {
          // odd
          add_amplitude(amplitudes2, basis2, -q_amp2);
        } else {
          // even
          add_amplitude(amplitudes2, basis2, q_amp2);
        }
      }
    }
    amplitudes = std::move(amplitudes2);
  }
  /**
   * /                                          \
   * | 1 0         |
   * | 0 exp(i phi)|
   * \             /
   */
  void gate_rz(qbit_id_t n, std::complex<double> const& phase) {
    amplitude_map_t amplitudes2;
    amplitudes2.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      add_amplitude(amplitudes2, basis,
                    is_1(basis, n) ? amplitude * phase : amplitude);
    }
    amplitudes = std::move(amplitudes2);
  }
  void gate_rz(qbit_id_t n, double theta) {
    std::complex<double> phase(cos(theta), sin(theta));
    gate_rz(n, phase);
  }
  double gate_measure(qbit_id_t n, bool is_up) {
    double p = 0.0;
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (is_1(basis, n) == is_up) {
        p += probability(amplitude);
      }
    }
    if (p > 0.0) {
      amplitude_map_t amplitudes2;
      double const w(1.0 / sqrt(p));
      for (amplitude_map_t::value_type const& v : amplitudes) {
        basis_t const basis(v.first);
        amplitude_t const& amplitude(v.second);
        if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
          continue;
        }
        if (is_1(basis, n) == is_up) {
          amplitudes2[basis] = amplitude / w;
        } else {
          amplitudes2[basis] = 0.0;
        }
      }
      amplitudes = std::move(amplitudes2);
    }
    return p;
  }
  bool gate_measure(qbit_id_t n) {
    double p1 = 0.0;
    double p0 = 0.0;
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (is_1(basis, n)) {
        p1 += probability(amplitude);
      } else {
        p0 += probability(amplitude);
      }
    }
    bool is_up = (p1 > random01.get_random()) ? true : false;
    amplitude_map_t amplitudes2;
    double const w(1.0 / sqrt(is_up ? p1 : p0));
    for (amplitude_map_t::value_type const& v : amplitudes) {
      basis_t const basis(v.first);
      amplitude_t const& amplitude(v.second);
      if (amplitude.real() == 0.0 && amplitude.imag() == 0.0) {
        continue;
      }
      if (is_1(basis, n) == is_up) {
        amplitudes2[basis] = amplitude * w;
      } else {
        amplitudes2[basis] = 0.0;
      }
    }
    amplitudes = std::move(amplitudes2);
    return is_up;
  }
  std::string dump(int n_qbits,  double limit) {
    std::vector<basis_t> bases;
    bases.reserve(amplitudes.size());
    for (amplitude_map_t::value_type const& v : amplitudes) {
      bases.push_back(v.first);
    }
    std::sort(bases.begin(), bases.end());
    std::string s;
    for (basis_t basis : bases) {
      amplitude_t const& amplitude(amplitudes[basis]);
      double const p = probability(amplitude);
      if (p < limit) {
        continue;
      }
      s += "|" + bitset_of(basis, n_qbits) + ">";
      s += " * ( " + std::to_string(amplitude.real()) + " , " +
        std::to_string(amplitude.imag()) + " )";
      s += " | " + std::to_string(p);
      s += "\n";
    }
    return s;
  }

 private:
  amplitude_map_t amplitudes;
};

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

static QC qc;

/**
 * 量子変数の状態を表示する。
 */
void dump(qbit_id_t n_qbits, double limit) {
  std::cout << qc.dump(n_qbits, limit);
}
void dump(std::string const& title, double limit) {
  std::cout << "===== " << title << " =====" << std::endl;
  std::cout << "-- qbits --" << std::endl;
  for (qbit& qbit : qbits) {
    std::cout << qbit.str() << std::endl;
  }
  dump(qbits.size(), limit);
}

/**
 * |000...000> の確率が 1, それ以外は 0 として
 * 量子変数を初期化する。
 */
void
reset() {
  qc.reset();
}

class frozen {
 public:
  explicit frozen(QC const& qc) : qc_(qc) {
  }
  QC const& get() const {
    return qc_;
  }
 private:
  QC const qc_;
};

frozen_ptr backup() {
  return std::make_shared<frozen>(qc);
}
void restore(frozen_ptr const& frozenptr) {
  qc = frozenptr->get();
}

// ----------------------------------------------------------------------

double
measure(qbit_id_t id, bool is_up) {
  return qc.gate_measure(id, is_up);
}

bool
measure(qbit_id_t id) {
  return qc.gate_measure(id);
}

void
hadamard(qbit_id_t target_id) {
  qc.gate_h(target_id);
}
void
hadamard_up_to(qbit_id_t n_qbits) {
  qc.gate_h_up_to(n_qbits);
}
void
cphase(qbit_id_t target_id, std::complex<double> const& phase) {
  qc.gate_rz(target_id, phase);
}
/**
 * /    \
 * | 0 1|
 * | 1 0|
 * \    /
 */
void pauli_x(qbit_id_t target_id) {
  qc.gate_x(target_id);
}
void cx(qbit_id_t control_id, qbit_id_t target_id) {
  qc.gate_cx(control_id, target_id);
}
void ccx(qbit_id_t control1_id, qbit_id_t control2_id,
         qbit_id_t target_id) {
  qc.gate_ccx(control1_id, control2_id, target_id);
}
/**
 * /     \
 * | 0 -i|
 * | i  0|
 * \     /
 */
void pauli_y(qbit_id_t target_id) {
  qc.gate_y(target_id);
}
void cy(qbit_id_t control_id,
        qbit_id_t target_id) {
  qc.gate_cy(control_id, target_id);
}
void ccy(qbit_id_t control1_id, qbit_id_t control2_id,
         qbit_id_t target_id) {
  qc.gate_ccy(control1_id, control2_id, target_id);
}
/**
 * /     \
 * | 1  0|
 * | 0 -1|
 * \     /
 */
void pauli_z(qbit_id_t target_id) {
  qc.gate_z(target_id);
}
void cz(qbit_id_t control_id,
        qbit_id_t target_id) {
  qc.gate_cz(control_id, target_id);
}
void ccz(qbit_id_t control1_id, qbit_id_t control2_id,
         qbit_id_t target_id) {
  qc.gate_ccz(control1_id, control2_id, target_id);
}
// ----------------------------------------------------------------------
double
measure(qbit const& q, bool is_up) {
  qbit_id_t const id(q.get_id());
  return measure(id, is_up);
}
bool
measure(qbit const& q) {
  qbit_id_t const id(q.get_id());
  return measure(id);
}
void
hadamard(qbit const& q) {
  qbit_id_t const target_id(q.get_id());
  hadamard(target_id);
}
void
hadamard_for_all() {
  hadamard_up_to(static_cast<qbit_id_t>(qbits.size()));
}

void
cphase(qbit const& q, std::complex<double> const& phase) {
  qbit_id_t const target_id(q.get_id());
  cphase(target_id, phase);
}
void pauli_x(qbit const& q) {
  qbit_id_t const target_id(q.get_id());
  pauli_x(target_id);
}
void cx(qbit const& control_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  cx(control1_id, target_id);
}
void ccx(qbit const& control1_q,
         qbit const& control2_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control1_q.get_id());
  qbit_id_t const control2_id(control2_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  ccx(control1_id, control2_id, target_id);
}
void pauli_y(qbit const& q) {
  qbit_id_t const target_id(q.get_id());
  pauli_y(target_id);
}
void cy(qbit const& control_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  cy(control1_id, target_id);
}
void ccy(qbit const& control1_q,
         qbit const& control2_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control1_q.get_id());
  qbit_id_t const control2_id(control2_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  ccy(control1_id, control2_id, target_id);
}
void pauli_z(qbit const& q) {
  qbit_id_t const target_id(q.get_id());
  pauli_z(target_id);
}
void cz(qbit const& control_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  cz(control1_id, target_id);
}
void ccz(qbit const& control1_q,
         qbit const& control2_q,
         qbit const& target_q) {
  qbit_id_t const control1_id(control1_q.get_id());
  qbit_id_t const control2_id(control2_q.get_id());
  qbit_id_t const target_id(target_q.get_id());
  ccz(control1_id, control2_id, target_id);
}

}  // namespace qc
