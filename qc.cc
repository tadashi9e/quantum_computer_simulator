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
typedef std::complex<double> amplitude_t;
typedef std::unordered_map<basis_t, amplitude_t> amplitude_map_t;

/**
 * 右から順に出現順に量子変数を並べたときの
 * |000...000> から |111...111> までの確率振幅。
 */
static amplitude_map_t q_amplitudes;
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

static bool is_1(basis_t n, qbit_id_t i) {
  return (n & (static_cast<basis_t>(1) << i)) ? true : false;
}

static basis_t set_1(basis_t n, qbit_id_t i) {
  return n | (static_cast<basis_t>(1) << i);
}

static basis_t set_0(basis_t n, qbit_id_t i) {
  return n & ~(static_cast<basis_t>(1) << i);
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
                          amplitude_t const& amplitude) {
  amplitude_map_t::iterator iter = amplitudes.find(basis);
  if (iter == amplitudes.end()) {
    amplitudes[basis] = amplitude;
    return;
  }
  amplitude_t& a(iter->second);
  a += amplitude;
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

static double amplitude(amplitude_t const& q_amp) {
  return q_amp.real() * q_amp.real() + q_amp.imag() * q_amp.imag();
}

/**
 * 量子変数の状態を表示する。
 */
void dump(qbit_id_t n_bits) {
  std::vector<basis_t> bases;
  bases.reserve(q_amplitudes.size());
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    bases.push_back(v.first);
  }
  std::sort(bases.begin(), bases.end());
  for (basis_t basis : bases) {
    amplitude_t const& q_amp(q_amplitudes[basis]);
    std::cout << bitset_of(basis, static_cast<basis_t>(1) << n_bits)
              << ": |"
              << q_amp
              << "| = " << amplitude(q_amp)
              << std::endl;
  }
}
void dump(std::string const& title) {
  std::cout << "===== " << title << " =====" << std::endl;
  std::cout << "-- qbits --" << std::endl;
  for (qbit& qbit : qbits) {
    std::cout << qbit.str() << std::endl;
  }
  dump(qbits.size());
}

/**
 * |000...000> の確率が 1, それ以外は 0 として
 * 量子変数を初期化する。
 */
void
reset() {
  q_amplitudes.clear();
  q_amplitudes[0] = amplitude_t(1, 0);
}

class frozen {
 public:
  explicit frozen(amplitude_map_t const& q_amp) : q_amp(q_amp) {
  }
  amplitude_map_t const& get() const {
    return q_amp;
  }
 private:
  amplitude_map_t q_amp;
};

frozen_ptr backup() {
  return std::make_shared<frozen>(q_amplitudes);
}
void restore(frozen_ptr const& frozenptr) {
  q_amplitudes = frozenptr->get();
}

// ----------------------------------------------------------------------

double
measure(qbit_id_t id, bool is_up) {
  double p = 0;
  basis_t mask = static_cast<basis_t>(0x01) << id;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    bool const b(((basis & mask) == 0) ? false : true);
    if (b == is_up) {
      p += amplitude(q_amp);
    }
  }
  if (p != 0) {
    double const w(1.0 / sqrt(p));
    for (amplitude_map_t::iterator iter = q_amplitudes.begin();
        iter != q_amplitudes.end();) {
      basis_t const basis(iter->first);
      amplitude_t& q_amp(iter->second);
      bool const b(((basis & mask) == 0) ? false : true);
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

bool
measure(qbit_id_t id) {
  double p0 = 0;
  double p1 = 0;
  basis_t mask = static_cast<basis_t>(0x01) << id;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    bool const b(((basis & mask) == 0) ? false : true);
    if (b) {
      p1 += amplitude(q_amp);
    } else {
      p0 += amplitude(q_amp);
    }
  }
  bool is_up((p1 > random01.get_random()) ? true : false);
  double const p(is_up ? p1 : p0);
  double const w(1.0 / sqrt(p));
  for (amplitude_map_t::iterator iter = q_amplitudes.begin();
       iter != q_amplitudes.end();) {
    basis_t const basis(iter->first);
    amplitude_t& q_amp(iter->second);
    bool const b(((basis & mask) == 0) ? false : true);
    if (b != is_up) {
      iter = q_amplitudes.erase(iter);
    } else {
      q_amp *= w;
      ++iter;
    }
  }
  return is_up;
}

/**
 * /     \
 * | 1  1| / sqrt(2)
 * | 1 -1|
 * \     /
 */
void
hadamard(qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  double const p(sqrt(2.0));
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    amplitude_t const q_amp2(q_amp / p);
    if (is_1(basis, target_id)) {
      add_amplitude(q_amplitudes2, set_0(basis, target_id), q_amp2);
      add_amplitude(q_amplitudes2, set_1(basis, target_id), -q_amp2);
    } else {
      add_amplitude(q_amplitudes2, set_0(basis, target_id), q_amp2);
      add_amplitude(q_amplitudes2, set_1(basis, target_id), q_amp2);
    }
  }
  q_amplitudes = std::move(q_amplitudes2);
}
void
hadamard_up_to(qbit_id_t n_qbits) {
  basis_t const total_basis(static_cast<basis_t>(0x01) << n_qbits);
  amplitude_map_t q_amplitudes2;
  double const p(pow(sqrt(2.0), n_qbits));
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    amplitude_t const q_amp2(q_amp / p);
    for (basis_t basis2 = 0; basis2 < total_basis; ++basis2) {
      if (count_1_bits(basis & basis2) & 0x01) {
        // odd
        add_amplitude(q_amplitudes2, basis2, -q_amp2);
      } else {
        // even
        add_amplitude(q_amplitudes2, basis2, q_amp2);
      }
    }
  }
  q_amplitudes = std::move(q_amplitudes2);
}
/**
 * /             \
 * | 1 0         |
 * | 0 exp(i phi)|
 * \             /
 */
void
cphase(qbit_id_t target_id, std::complex<double> const& phase) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    amplitude_t const& q_amp2(is_1(basis, target_id) ?
                              (q_amp * phase) : q_amp);
    add_amplitude(q_amplitudes2, basis, q_amp2);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
/**
 * /    \
 * | 0 1|
 * | 1 0|
 * \    /
 */
void pauli_x(qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    basis_t const basis2(is_1(basis, target_id) ?
                         set_0(basis, target_id) : set_1(basis, target_id));
    add_amplitude(q_amplitudes2, basis2, q_amp);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
void cx(qbit_id_t control_id, qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    if (!is_1(basis, control_id)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    basis_t const basis2(is_1(basis, target_id) ?
                         set_0(basis, target_id) : set_1(basis, target_id));
    add_amplitude(q_amplitudes2, basis2, q_amp);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
void ccx(qbit_id_t control1_id, qbit_id_t control2_id,
         qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    if (!is_1(basis, control1_id)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    if (!is_1(basis, control2_id)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    basis_t const basis2(is_1(basis, target_id) ?
                         set_0(basis, target_id) : set_1(basis, target_id));
    add_amplitude(q_amplitudes2, basis2, q_amp);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
/**
 * /     \
 * | 0 -i|
 * | i  0|
 * \     /
 */
void pauli_y(qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    basis_t const basis2(is_1(basis, target_id) ?
                         set_0(basis, target_id) : set_1(basis, target_id));
    amplitude_t const q_amp2(q_amp *
                             (is_1(basis, target_id) ?
                              amplitude_t(0.0, -1.0) : amplitude_t(0.0, 1.0)));
    add_amplitude(q_amplitudes2, basis2, q_amp2);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
void cy(qbit_id_t control_id,
        qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    if (!is_1(basis, control_id)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    basis_t const basis2(is_1(basis, target_id) ?
                         set_0(basis, target_id) : set_1(basis, target_id));
    amplitude_t const q_amp2(q_amp *
                             (is_1(basis, target_id) ?
                              amplitude_t(0.0, -1.0) : amplitude_t(0.0, 1.0)));
    add_amplitude(q_amplitudes2, basis2, q_amp2);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
void ccy(qbit_id_t control_id1, qbit_id_t control_id2,
         qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    if (!is_1(basis, control_id1)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    if (!is_1(basis, control_id2)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    basis_t const basis2(is_1(basis, target_id) ?
                         set_0(basis, target_id) : set_1(basis, target_id));
    amplitude_t const q_amp2(q_amp *
                             (is_1(basis, target_id) ?
                              amplitude_t(0.0, -1.0) : amplitude_t(0.0, 1.0)));
    add_amplitude(q_amplitudes2, basis2, q_amp2);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
/**
 * /     \
 * | 1  0|
 * | 0 -1|
 * \     /
 */
void pauli_z(qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    amplitude_t const q_amp2(is_1(basis, target_id) ? -q_amp : q_amp);
    add_amplitude(q_amplitudes2, basis, q_amp2);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
void cz(qbit_id_t control_id,
        qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    if (!is_1(basis, control_id)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    amplitude_t const q_amp2(is_1(basis, target_id) ? -q_amp : q_amp);
    add_amplitude(q_amplitudes2, basis, q_amp2);
  }
  q_amplitudes = std::move(q_amplitudes2);
}
void ccz(qbit_id_t control1_id, qbit_id_t control2_id,
         qbit_id_t target_id) {
  amplitude_map_t q_amplitudes2;
  for (amplitude_map_t::value_type const& v : q_amplitudes) {
    basis_t const basis(v.first);
    amplitude_t const& q_amp(v.second);
    if (!is_1(basis, control1_id)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    if (!is_1(basis, control2_id)) {
      add_amplitude(q_amplitudes2, basis, q_amp);
      continue;
    }
    amplitude_t const q_amp2(is_1(basis, target_id) ? -q_amp : q_amp);
    add_amplitude(q_amplitudes2, basis, q_amp2);
  }
  q_amplitudes = std::move(q_amplitudes2);
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
