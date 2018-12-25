
量子計算機シミュレータ
================

- qc.h  ヘッダファイル
- qc.cc ソースファイル
- test_teleport.cc テストプログラム

前提
====

C++ コンパイラおよび boost ライブラリ

使用方法
======

qc.h をヘッダファイルとして以下の内容のシミュレーションプログラムを記述する。

1. qc::qbit 型変数として、使用する量子変数を宣言する。
2. 量子変数を初期化する。すべての量子変数が |0> 状態となる。
3. 量子変数を操作する。
5. 観測操作を行い、結果を得る。

記述したプログラムを qc.cc とともにコンパイルする。

量子ビットの宣言
============

コンストラクタ qc::qbit() あるいは 
qc::qbit(std::string const& name) で量子ゲートを宣言する。
引数に名前を与えられた場合には、デバッグ時にその名前付きで表示する。

初期化
======

```c++
void qc::reset()
```

|000...000> の確率を 1 にする。それ以外は 0。

デバッグ表示
=========

```c++
void dump(std::string const& title);
```

量子ビットの一覧、および、その時点での各状態の振幅を表示する。

演算操作
======

アダマールゲート
------------

```c++
void hadamard(qbit const& q);
```

量子変数 q にアダマール変換を行う。

パウリ X ゲート
-----------

```c++
void pauli_x(qbit const& q);
```

パウリ Y ゲート
-----------

```c++
void pauli_y(qbit const& q);
```

パウリ Z ゲート
-----------

```c++
void pauli_z(qbit const& q);
```

位相ゲート
-------

```c++
void cphase(qbit const& q, std::complex<double> const& phase);
```

ここでは phase には位相 $\theta$ ではなく ${e}^{i \theta}$ を与える。

制御ノットゲート
-----------

```c++
void cnot(qbit const& control_q, qbit const& target_q);
```

観測
===

```c++
bool measure(qbit const& q);
```

観測を行い、その結果をブール値で返す。
