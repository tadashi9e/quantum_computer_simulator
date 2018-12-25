#include "qc.h"
#include <sys/time.h>
#include <iostream>
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

/**
 * 右から順に出現順に量子変数を並べたときの
 * |000...000> から |111...111> までの確率振幅。
 */
static std::vector<std::complex<double> > q_amplitudes;
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
  while(q_amplitudes.size() < ((size_t)1 << (id + 1))) {
    q_amplitudes.push_back(std::complex<double>());
  }
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
  for(size_t i = 0;;++i) {
    result += (n & ((size_t)1 << i)) ? '1' : '0';
    pow_i *= 2;
    if (pow_i >= max_n) {
      break;
    }
  }
  return result;
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
  for(size_t i = 0;i < q_amplitudes.size();++i) {
    std::cout << bitset_of(i, q_amplitudes.size())
              << ": |"
              << q_amplitudes[i]
              << "| = " << amplitude(q_amplitudes[i])
              << std::endl;
  }
}

/**
 * |000...000> の確率が 1, それ以外は 0 として
 * 量子変数を初期化する。
 */
void
reset() {
  if (q_amplitudes.size() == 0) {
    return;
  }
  q_amplitudes[0] = std::complex<double>(1, 0);
  for(size_t i = 1;i < q_amplitudes.size();++i) {
    q_amplitudes[i] = std::complex<double>();
  }
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
  return (i == 0 && j == 0) ? std::complex<double>(1, 0) : std::complex<double>();
}

static std::complex<double> op_upside(size_t j, size_t i) {
  return (i == 1 && j == 1) ? std::complex<double>(1, 0) : std::complex<double>();
}

/**
 * 二演算子のテンソル積
 */
static std::complex<double>
tensor_product(boost::function<std::complex<double>(size_t, size_t)> op1,
               boost::function<std::complex<double>(size_t, size_t)> op2,
               size_t j, size_t i) {
  assert(j < 4);
  assert(i < 4);
  return op1(j >> 1, i >> 1) * op2(j & 0x01, i & 0x01);
}

/**
 * CNOT ゲート
 */
static std::complex<double> op_cnot(size_t j, size_t i) {
  return tensor_product(&op_downside, &op_unit, i, j)
    + tensor_product(&op_upside, &op_pauli_x, i, j);
}

/**
 * 与えられた演算子を適用する
 */
static std::vector<std::complex<double> >
mat_multiply(boost::function<std::complex<double>(size_t, size_t)> op,
             std::vector<std::complex<double> > const& values) {
  std::vector<std::complex<double> > result;
  for(size_t i = 0;i < values.size();++i) {
    result.push_back(std::complex<double>());
  }
  for(size_t j = 0;j < values.size();++j) {
    std::complex<double> c;
    for(size_t i = 0;i < values.size();++i) {
      c += op(j, i) * values[i];
    }
    result[j] = c;
  }
  return result;
}

/**
 * id 癌目の量子変数に対して op を適用し、それ以外は I を適用する
 * 演算子のテンソル積
 *
 * I (X) I (X) .....I (X) OP (X) I.....(X) I (X) I
 */
static std::complex<double>
op_i_op1_i(boost::function<std::complex<double>(size_t, size_t)> op1,
                   int id, size_t j, size_t i) {
  size_t mask = ~(static_cast<size_t>(0x01) << id);
  if ((mask & i) != (mask & j)) {
    return std::complex<double>();
  }
  size_t j1 = (j >> id) & 0x01;
  size_t i1 = (i >> id) & 0x01;
  return op1(j1, i1);
}

/**
 * id1 癌目の量子変数に対して op1 を
 * id2 癌目の量子変数に対して op2 を適用し、それ以外は I を適用する
 * 演算子のテンソル積
 *
 * I (X) ... I (X) OP1 (X) I ...(X) I (X) OP2 (X) I ... (X) I
 */
static std::complex<double>
op_i_op2_i(boost::function<std::complex<double>(size_t,size_t)> op2,
                   int id2, int id1, size_t j, size_t i) {
  size_t mask2 = static_cast<size_t>(0x01) << id2;
  size_t mask1 = static_cast<size_t>(0x01) << id1;
  size_t mask = ~(mask1 | mask2);
  if ((mask & i) != (mask & j)) {
    return std::complex<double>();
  }
  size_t j2 = ((j & mask2) == 0) ? 0 : 1;
  size_t i2 = ((i & mask2) == 0) ? 0 : 1;
  size_t j1 = ((j & mask1) == 0) ? 0 : 1;
  size_t i1 = ((i & mask1) == 0) ? 0 : 1;
  return op2((i2 << 1) | i1, (j2 << 1) | j1);
}

void
dump_operator(boost::function<std::complex<double>(size_t, size_t)> op) {
  for(size_t j = 0;j < q_amplitudes.size();++j) {
    for(size_t i = 0;i < q_amplitudes.size();++i) {
      std::cout << op(j, i) << " , ";
    }
    std::cout << std::endl;
  }
}

static double
op_measure(int id, bool is_up) {
  double p = 0;
  size_t mask = static_cast<size_t>(0x01) << id;
  for(size_t n = 0;n < q_amplitudes.size();++n) {
    bool b = ((n & mask) == 0) ? false : true;
    if (b == is_up) {
      p += amplitude(q_amplitudes[n]);
    }
  }
  if (p != 0) {
    double w = 1.0 / sqrt(p);
    for(size_t n = 0;n < q_amplitudes.size();++n) {
      bool b = ((n & mask) == 0) ? false : true;
      if (b != is_up) {
        q_amplitudes[n] = std::complex<double>();
      } else {
        q_amplitudes[n] *= w;
      }
    }
  } else {
    double w = 1.0 / sqrt(q_amplitudes.size());
    for(size_t n = 0;n < q_amplitudes.size();++n) {
      bool b = ((n & mask) == 0) ? false : true;
      if (b != is_up) {
        q_amplitudes[n] = std::complex<double>();
      } else {
        q_amplitudes[n] = std::complex<double>(w, 0);
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
  for(size_t n = 0;n < q_amplitudes.size();++n) {
    bool b = ((n & mask) == 0) ? false : true;
    double a = amplitude(q_amplitudes[n]);
    if (b) {
      p1 += a;
    } else {
      p0 += a;
    }
  }
  bool is_up = (p1 > random01.get_random()) ? true : false;
  double p = is_up ? p1 : p0;
  double w = 1.0 / sqrt(p);
  for(size_t n = 0;n < q_amplitudes.size();++n) {
    bool b = ((n & mask) == 0) ? false : true;
    if (b != is_up) {
      q_amplitudes[n] = std::complex<double>();
    } else {
      q_amplitudes[n] *= w;
    }
  }
  return is_up;
}

double
measure(qbit const& q, bool is_up) {
  int id = q.get_id();
  return op_measure(id, is_up);
}

bool
measure(qbit const& q) {
  int id = q.get_id();
  return op_measure(id);
}

void
hadamard(qbit const& q) {
  int id = q.get_id();
  q_amplitudes = mat_multiply(boost::bind(&op_i_op1_i, &op_hadamard, id, _1, _2),
                              q_amplitudes);
}
void
dump_hadamard(qbit const& q) {
  int id = q.get_id();
  dump_operator(boost::bind(&op_i_op1_i, op_hadamard, id, _1, _2));
}

void
cphase(qbit const& q, std::complex<double> const& phase) {
  int id = q.get_id();
  boost::function<std::complex<double>(size_t, size_t)> op = boost::bind(&op_cphase, phase, _1, _2);
  q_amplitudes = mat_multiply(boost::bind(&op_i_op1_i,
                                      op, id, _1, _2),
                              q_amplitudes);
}

void
dump_cphase(qbit const& q, std::complex<double> const& phase) {
  int id = q.get_id();
  boost::function<std::complex<double>(size_t, size_t)>
    op = boost::bind(&op_cphase, phase, _1, _2);
  dump_operator(boost::bind(&op_i_op1_i,
                            op, id, _1, _2));
}

void pauli_x(qbit const& q) {
  int id = q.get_id();
  q_amplitudes = mat_multiply(boost::bind(&op_i_op1_i, op_pauli_x, id, _1, _2),
                              q_amplitudes);
}
void dump_pauli_x(qbit const& q) {
  int id = q.get_id();
  dump_operator(boost::bind(&op_i_op1_i, op_pauli_x, id, _1, _2));
}

void pauli_y(qbit const& q) {
  int id = q.get_id();
  q_amplitudes = mat_multiply(boost::bind(&op_i_op1_i, op_pauli_y, id, _1, _2),
                          q_amplitudes);
}
void dump_pauli_y(qbit const& q) {
  int id = q.get_id();
  dump_operator(boost::bind(&op_i_op1_i, op_pauli_y, id, _1, _2));
}

void pauli_z(qbit const& q) {
  int id = q.get_id();
  q_amplitudes = mat_multiply(boost::bind(&op_i_op1_i, op_pauli_z, id, _1, _2),
                              q_amplitudes);
}
void dump_pauli_z(qbit const& q) {
  int id = q.get_id();
  dump_operator(boost::bind(&op_i_op1_i, op_pauli_z, id, _1, _2));
}

void cnot(qbit const& control_q, qbit const& target_q) {
  int control_id = control_q.get_id();
  int target_id = target_q.get_id();
  q_amplitudes = mat_multiply(boost::bind(&op_i_op2_i, op_cnot,
                                          control_id, target_id,
                                          _1, _2),
                              q_amplitudes);
}
void dump_cnot(qbit const& control_q, qbit const& target_q) {
  int control_id = control_q.get_id();
  int target_id = target_q.get_id();
  dump_operator(boost::bind(&op_i_op2_i, op_cnot, control_id, target_id,
                            _1, _2));
}

}
