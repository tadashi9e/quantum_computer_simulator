// -*- coding:utf-8-unix;mode:c++ -*-
// Copyright [2019] tadashi9e
#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include "qc.h"

/**
 * 3bits の足し算
 */
static void
q_add3(qc::qbit const& a,
       qc::qbit const& b,
       qc::qbit const& c,
       qc::qbit const& s0,
       qc::qbit const& s1) {
  qc::cx(a, s0);
  qc::cx(b, s0);
  qc::cx(c, s0);
  qc::ccx(a, b, s1);
  qc::ccx(b, c, s1);
  qc::ccx(c, a, s1);
}

/**
 * 加算に利用した量子変数を 0 に戻す
 */
void
q_add3_cleanup(qc::qbit const& a,
               qc::qbit const& b,
               qc::qbit const& c,
               qc::qbit const& t) {
  qc::ccx(a, b, t);
  qc::ccx(b, c, t);
  qc::ccx(c, a, t);
}

/**
 * 与えられた量子ビットについて加算を実行する。
 *
 * t はすべて 0 に初期化されていなければならない。
 */
void
q_add(std::vector<qc::qbit> const& a,
      std::vector<qc::qbit> const& b,
      std::vector<qc::qbit> const& s,
      std::vector<qc::qbit> const& t) {
  int const n(std::min(a.size(), b.size()));
  for (int i = 0; i < n; ++i) {
    if (i < n - 1) {
      q_add3(a[i], b[i], t[i],
             s[i], t[i + 1]);
    } else {
      q_add3(a[i], b[i], t[i],
             s[i], s[i + 1]);
    }
  }
  // ゴミを消す
  for (int i = n-2; i >= 0; --i) {
    q_add3_cleanup(a[i], b[i], t[i], t[i + 1]);
  }
}

/**
 * 与えられた量子変数全てについて測定を行い、数値として返す。
 */
int
as_c_number(std::vector<qc::qbit> const& qs) {
  int n = 0;
  for (int i = qs.size()-1; i >= 0; --i) {
    n *= 2;
    if (qc::measure(qs[i])) {
      ++n;
    }
  }
  return n;
}

/**
 * 0+0=0
 * 0+1=1
 * 1+0=1
 * 1+1=2
 * を同時に実行する。
 *
 * 実行例:
 * $ ./test_add|sort |uniq -c
 *    26 0+0=0
 *    27 0+1=1
 *    25 1+0=1
 *    22 1+1=2
 */
int
main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "usage: " << argv[0] << " <nbits> <nloop>" << std::endl;
    exit(1);
  }
  int const nbits(atoi(argv[1]));
  int const nloop(atoi(argv[2]));
  std::vector<qc::qbit> a;
  std::vector<qc::qbit> b;
  qc::reset();
  a.resize(nbits);
  b.resize(nbits);
  qc::hadamard_for_all();
  std::vector<qc::qbit> t;
  std::vector<qc::qbit> s;
  t.resize(nbits);
  s.resize(nbits + 1);
  // a + b = c を計算する
  // qc::dump("add start");
  q_add(a, b, s, t);
  // qc::dump("after add");

  qc::frozen_ptr backup = qc::backup();
  for (int i = 0; i < nloop; ++i) {
    std::cout << as_c_number(a) << "+" << as_c_number(b) <<
      "=" << as_c_number(s) << std::endl;
    qc::restore(backup);
  }
}
