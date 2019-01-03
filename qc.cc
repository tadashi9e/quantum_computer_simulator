// -*- coding:utf-8-unix;mode:c++ -*-
// Copyright [2019] tadashi9e
#include "qc.h"
#include <sys/time.h>
#include <iostream>
#include <map>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

namespace qc {

class Random01 {
  boost::random::random_device rnd_dev;
  boost::random::mt19937 rnd_gen;
 public:
  Random01() : rnd_gen(rnd_dev()) {
  }
  double get_random() {
    boost::random::uniform_01<double> rnd01;
    double rnd = rnd01(rnd_gen);
    return rnd;
  }
};

Random01 random01;

typedef std::map<size_t, std::complex<double> > amplitudes_t;
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
static int
new_qbit_id() {
  int id = qbits.size();
  return id;
}

struct qbit::q_impl {
  int id;
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
  : impl(boost::make_shared<q_impl>()) {
  impl->setup();
  qbits.push_back(*this);
}

qbit::qbit(std::string const& name)
  : impl(boost::make_shared<q_impl>()) {
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
int qbit::get_id() const {
  return impl->id;
}

std::string qbit::str() const {
  return std::string("[") + boost::lexical_cast<std::string>(impl->id)
    + "]" + impl->name;
}

static std::string bitset_of(size_t n, size_t max_n) {
  std::string result;
  size_t pow_i = 1;
  for (size_t i = 0; ; ++i) {
    result += (n & ((size_t)1 << i)) ? '1' : '0';
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
  BOOST_FOREACH(qbit& qbit, qbits) {
    std::cout << qbit.str() << std::endl;
  }
  std::cout << "-- states --" << std::endl;
  BOOST_FOREACH(amplitudes_t::value_type const& v, q_amplitudes) {
    size_t basis(v.first);
    std::complex<double> const& q_amp(v.second);
    std::cout << bitset_of(basis, 1 << qbits.size())
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
  return boost::make_shared<frozen>(q_amplitudes);
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
static std::complex<double> op_unit(size_t j, size_t i) {
  return (i == j)
    ? std::complex<double>(1, 0) : std::complex<double>();
}

/**
 * /     \
 * | 1  1| / sqrt(2)
 * | 1 -1|
 * \     /
 */
static std::complex<double> op_hadamard(size_t j, size_t i) {
  assert(i < 2);
  assert(j < 2);
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
op_cphase(std::complex<double> const& phase, size_t j, size_t i) {
  assert(i < 2);
  assert(j < 2);
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
static std::complex<double> op_pauli_x(size_t j, size_t i) {
  assert(i < 2);
  assert(j < 2);
  return (i != j) ? std::complex<double>(1, 0) : std::complex<double>();
}

/**
 * /     \
 * | 0 -i|
 * | i  0|
 * \     /
 */
static std::complex<double> op_pauli_y(size_t j, size_t i) {
  assert(i < 2);
  assert(j < 2);
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
static std::complex<double> op_pauli_z(size_t j, size_t i) {
  assert(i < 2);
  assert(j < 2);
  if (i != j) {
    return std::complex<double>();
  }
  if (i == 0) {
    return std::complex<double>(1, 0);
  } else {
    return std::complex<double>(-1, 0);
  }
}

static std::complex<double> op_downside(size_t j, size_t i) {
  return (i == 0 && j == 0)
    ? std::complex<double>(1, 0)
    : std::complex<double>();
}

static std::complex<double> op_upside(size_t j, size_t i) {
  return (i == 1 && j == 1)
    ? std::complex<double>(1, 0)
    : std::complex<double>();
}

/**
 * 二演算子のテンソル積
 */
static std::complex<double>
tensor_product(boost::function<std::complex<double>(size_t, size_t)> op1,
               size_t j1, size_t i1,
               boost::function<std::complex<double>(size_t, size_t)> op2,
               size_t j2, size_t i2) {
  return op1(j1, i1) * op2(j2, i2);
}

/**
 * 三演算子のテンソル積
 */
static std::complex<double>
tensor_product(boost::function<std::complex<double>(size_t, size_t)> op1,
               size_t j1, size_t i1,
               boost::function<std::complex<double>(size_t, size_t)> op2,
               size_t j2, size_t i2,
               boost::function<std::complex<double>(size_t, size_t)> op3,
               size_t j3, size_t i3) {
  return (op1(j1, i1) *
          op2(j2, i2) *
          op3(j3, i3));
}

/**
 * CNOT ゲート
 */
static std::complex<double> op_cx(size_t j1, size_t i1,
                                  size_t j2, size_t i2) {
  return (tensor_product(&op_downside, j1, i1, &op_unit,    j2, i2) +
          tensor_product(&op_upside,   j1, i1, &op_pauli_x, j2, i2));
}

/**
 * Toffoli ゲート
 */
static std::complex<double> op_ccx(size_t j1, size_t i1,
                                   size_t j2, size_t i2,
                                   size_t j3, size_t i3) {
  return (tensor_product(&op_downside, j1, i1,
                         &op_downside, j2, i2,
                         &op_unit,     j3, i3) +
          tensor_product(&op_upside,   j1, i1,
                         &op_downside, j2, i2,
                         &op_unit,     j3, i3) +
          tensor_product(&op_downside, j1, i1,
                         &op_upside,   j2, i2,
                         &op_unit,     j3, i3) +
          tensor_product(&op_upside,   j1, i1,
                         &op_upside,   j2, i2,
                         &op_pauli_x,  j3, i3));
}

static amplitudes_t
apply_tensor_product(boost::function<std::complex<double>(size_t, size_t)> op,
                     int id,
                     amplitudes_t amplitudes) {
  amplitudes_t amplitudes2;
  size_t imask(~(static_cast<size_t>(0x01) << id));
  BOOST_FOREACH(amplitudes_t::value_type const& v, amplitudes) {
    size_t const basis(v.first);
    std::complex<double> const& amplitude(v.second);
    size_t i((basis >> id) & 0x01);
    for (size_t j = 0; j < 2; ++j) {
      size_t const target_id((basis & imask) | (j << id));
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
apply_tensor_product(boost::function<std::complex<double>(size_t, size_t,
                                                          size_t, size_t)> op2,
                     int id1, int id2,
                     amplitudes_t amplitudes) {
  amplitudes_t amplitudes2;
  size_t imask(~((static_cast<size_t>(0x01) << id1) |
                 (static_cast<size_t>(0x01) << id2)));
  BOOST_FOREACH(amplitudes_t::value_type const& v, amplitudes) {
    size_t const basis(v.first);
    std::complex<double> const& amplitude(v.second);
    size_t const i1((basis >> id1) & 0x01);
    size_t const i2((basis >> id2) & 0x01);
    for (size_t j2 = 0; j2 < 2; ++j2) {
      for (size_t j1 = 0; j1 < 2; ++j1) {
        std::complex<double> w(op2(j1, i1, j2, i2));
        if (w == std::complex<double>()) {
          continue;
        }
        size_t const target_id((basis & imask) | (j1 << id1) | (j2 << id2));
        amplitudes2[target_id] += w * amplitude;
      }
    }
  }
  return amplitudes2;
}

static amplitudes_t
apply_tensor_product(boost::function<std::complex<double>(size_t, size_t,
                                                          size_t, size_t,
                                                          size_t, size_t)> op3,
                     int id1, int id2, int id3,
                     amplitudes_t amplitudes) {
  amplitudes_t amplitudes2;
  size_t imask(~((static_cast<size_t>(0x01) << id1) |
                 (static_cast<size_t>(0x01) << id2) |
                 (static_cast<size_t>(0x01) << id3)));
  BOOST_FOREACH(amplitudes_t::value_type const& v, amplitudes) {
    size_t basis(v.first);
    std::complex<double> const& amplitude(v.second);
    size_t const i1((basis >> id1) & 0x01);
    size_t const i2((basis >> id2) & 0x01);
    size_t const i3((basis >> id3) & 0x01);
    for (size_t j3 = 0; j3 < 2; ++j3) {
      for (size_t j2 = 0; j2 < 2; ++j2) {
        for (size_t j1 = 0; j1 < 2; ++j1) {
          std::complex<double> w(op3(j1, i1, j2, i2, j3, i3));
          if (w == std::complex<double>()) {
            continue;
          }
          size_t const target_id((basis & imask) |
                                 (j1 << id1) | (j2 << id2) | (j3 << id3));
          amplitudes2[target_id] += w * amplitude;
        }
      }
    }
  }
  return amplitudes2;
}

static double
op_measure(int id, bool is_up) {
  double p = 0;
  size_t mask = static_cast<size_t>(0x01) << id;
  BOOST_FOREACH(amplitudes_t::value_type const& v, q_amplitudes) {
    size_t const basis(v.first);
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
      size_t const basis(iter->first);
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
op_measure(int id) {
  double p0 = 0;
  double p1 = 0;
  size_t mask = static_cast<size_t>(0x01) << id;
  BOOST_FOREACH(amplitudes_t::value_type const& v, q_amplitudes) {
    size_t const basis(v.first);
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
    size_t const basis(iter->first);
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
  int const id(q.get_id());
  return op_measure(id, is_up);
}

bool
measure(qbit const& q) {
  int const id(q.get_id());
  return op_measure(id);
}

void
hadamard(qbit const& q) {
  int const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_hadamard, id, q_amplitudes);
}
void
cphase(qbit const& q, std::complex<double> const& phase) {
  int const id(q.get_id());
  boost::function<std::complex<double>(size_t, size_t)>
    op(boost::bind(&op_cphase, phase, _1, _2));
  q_amplitudes = apply_tensor_product(op, id, q_amplitudes);
}
void pauli_x(qbit const& q) {
  int const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_pauli_x, id, q_amplitudes);
}
void pauli_y(qbit const& q) {
  int const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_pauli_y, id, q_amplitudes);
}
void pauli_z(qbit const& q) {
  int const id(q.get_id());
  q_amplitudes = apply_tensor_product(op_pauli_z, id, q_amplitudes);
}
void cx(qbit const& control_q, qbit const& target_q) {
  int const control_id(control_q.get_id());
  int const target_id(target_q.get_id());
  q_amplitudes = apply_tensor_product(op_cx, control_id, target_id,
                                      q_amplitudes);
}
void ccx(qbit const& control1_q,
         qbit const& control2_q,
         qbit const& target_q) {
  int const control1_id(control1_q.get_id());
  int const control2_id(control2_q.get_id());
  int const target_id(target_q.get_id());
  q_amplitudes = apply_tensor_product(op_ccx, control1_id, control2_id,
                                      target_id,
                                      q_amplitudes);
}

}  // namespace qc
