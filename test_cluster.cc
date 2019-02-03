// -*- coding:utf-8-unix;mode:c++ -*-
// Copyright [2019] tadashi9e
#include <iostream>
#include <vector>
#include "qc.h"

int const SIZE_X(4);
int const SIZE_Y(4);

template<int SZ_Y, int SZ_X>
void make_cluster_2d(qc::qbit (&v)[SZ_Y][SZ_X]) {
  for (size_t y = 0; y < SZ_Y; ++y) {
    for (size_t x = 0; x < SZ_X; ++x) {
      if (x + 1 < SZ_X) {
        qc::cz(v[y][x], v[y    ][x + 1]);
      }
      if (y + 1 < SZ_Y) {
        qc::cz(v[y][x], v[y + 1][x    ]);
      }
    }
  }
}

int
main(int argc, char* argv[]) {
  qc::qbit v[SIZE_Y][SIZE_X];
  qc::reset();
  qc::hadamard_for_all();
  make_cluster_2d(v);
  qc::dump("2D");
}
