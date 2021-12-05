/*
    This file is part of TON Blockchain Library.

    TON Blockchain Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    TON Blockchain Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with TON Blockchain Library.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <assert.h>
#include <string.h>
#include <array>
#include <string>
#include <iostream>
#include <sstream>
#include "common/refcnt.hpp"
#include "common/bigint.hpp"
#include "common/refint.h"
#include "modbigint.cpp"

#include "td/utils/tests.h"

using BInt = modint::ModArray<18>;  // integers up to 2^537
using MRInt = modint::MixedRadix<18>;  // auxiliary integer representation for printing, comparing etc

MRInt p2_256, np2_256, p2_63, np2_63;
constexpr long long ll_min = -2 * (1LL << 62), ll_max = ~ll_min;

int mkint_chk_mode = -1, res_chk_mode = 0;

bool mr_in_range(const MRInt& x) {
  return x < p2_256 && x >= np2_256;
}

bool mr_is_small(const MRInt& x) {
  return x < p2_63 && x >= np2_63;
}

bool mr_fits_bits(const MRInt& x, int bits) {
  if (bits > 0) {
    return x < MRInt::pow2(bits - 1) && x >= MRInt::negpow2(bits - 1);
  } else {
    return !bits && !x.sgn();
  }
}

bool mr_ufits_bits(const MRInt& x, int bits) {
  return bits >= 0 && x.sgn() >= 0 && x < MRInt::pow2(bits);
}

bool extract_value_any_bool(BInt& val, const td::AnyIntView<td::BigIntInfo>& x, bool chk_norm = true) {
  int n = x.size();
  if (n <= 0 || n > x.max_size() || (!x.digits[n - 1] && n > 1)) {
    return false;
  }
  assert (n == 1 || x.digits[n - 1] != 0);
  val.set_zero();
  for (int i = n - 1; i >= 0; --i) {
    val.lshift_add(td::BigIntInfo::word_shift, x.digits[i]);
    if (chk_norm && (x.digits[i] < - td::BigIntInfo::Half || x.digits[i] >= td::BigIntInfo::Half)) {
      return false;   // unnormalized
    }
  }
  return true;
}

template <typename T>
bool extract_value_bool(BInt& val, const T& x, bool chk_norm = true) {
  return extract_value_any_bool(val, x.as_any_int(), chk_norm);
}

BInt extract_value_any(const td::AnyIntView<td::BigIntInfo>& x, bool chk_norm = true) {
  BInt res;
  CHECK(extract_value_any_bool(res, x, chk_norm));
  return res;
}

template <typename T>
BInt extract_value(const T& x, bool chk_norm = true) {
  return extract_value_any(x.as_any_int(), chk_norm);
}

template <typename T>
BInt extract_value_alt(const T& x) {
  BInt res;
  const int *md = res.mod_array();
  for (int i = 0; i < res.n / 2; i++) {
    T copy{x};
    int m1 = md[2 * i], m2 = md[2 * i + 1];
    long long rem = copy.divmod_short((long long)m1 * m2);
    res.a[2 * i] = (int)(rem % m1); 
    res.a[2 * i + 1] = (int)(rem % m2);
  }
  if (res.n & 1) {
    T copy{x};
    res.a[res.n - 1] = (int)copy.divmod_short(md[res.n - 1]);
  }
  return res;
}

constexpr int min_spec_int = -0xfd08, max_spec_int = 0xfd07;
// x = sgn*(ord*256+a*16+b) => sgn*((32+a)*2^(ord-2) + b - 8)
// x = -0xfd08 => -2^256 ... x = 0xfd07 => 2^256 - 1
td::RefInt256 make_special_int(int x, BInt *ptr = nullptr, unsigned char bin[64] = nullptr) {
  bool sgn = (x < 0);
  if (sgn) { x = -x; }
  int ord = (x >> 8) - 2, a = 32 + ((x >> 4) & 15), b = (x & 15) - 8;
  if (ord < 0) {
    a >>= -ord;
    ord = 0;
  }
  if (sgn) { a = -a;  b = -b; }
  if (ptr) {
    ptr->set_int(a);
    *ptr <<= ord;
    *ptr += b;
  }
  if (bin) {
    int acc = b, r = ord;
    for (int i = 63; i >= 0; --i) {
      if (r < 8) {
	acc += (a << r);
	r = 1024;
      }
      r -= 8;
      bin[i] = (unsigned char)(acc & 0xff);
      acc >>= 8;
    }
  }
  return (td::make_refint(a) << ord) + b;
}

void check_one_int_repr(td::RefInt256 x, int mode, int in_range, const BInt *valptr = nullptr, const unsigned char bin[64] = nullptr) {
  CHECK(x.not_null() && (in_range <= -2 || x->is_valid()));
  if (!x->is_valid()) {
    // not much to check when x is a NaN
    unsigned char bytes[64];
    if (valptr) {
      // check that the true answer at `valptr` is out of range
      CHECK(!mr_in_range(valptr->to_mixed_radix()));
      if (mode & 0x200) {
	// check BInt binary export
	valptr->to_binary(bytes, 64);
	if (bin) {
	  // check that the two true answers match
	  CHECK(!memcmp(bin, bytes, 64));
	} else {
	  bin = bytes;
	}
      }
    }
    if (bin) {
      // check that the true answer in `bin` is out of range
      int i = 0, sgn = (bin[0] >= 0x80 ? -1 : 0);
      while (i < 32 && bin[i] == (unsigned char)sgn);
      CHECK(i < 32);
      if (valptr && (mode & 0x100)) {
	// check BInt binary export
        BInt val2;
        val2.from_binary(bin, 64);
        CHECK(*valptr == val2);
      }
    }
    return;
  }
  unsigned char bytes[64];
  CHECK(x->export_bytes(bytes, 64));
  if (bin) {
    CHECK(!memcmp(bytes, bin, 64));
  }
  BInt val = extract_value(*x);
  if (valptr) {
    CHECK(val == *valptr);
  }
  if (mode & 1) {
    BInt val2 = extract_value_alt(*x);
    CHECK(val == val2);
  }
  MRInt mval(val);
  bool val_in_range = mr_in_range(mval);
  CHECK(x->fits_bits(257) == val_in_range);
  if (in_range >= 0) {
    CHECK((int)val_in_range == in_range);
  }
  if (mode & 2) {
    // check binary import
    td::BigInt256 y;
    y.import_bytes(bytes, 64);
    CHECK(y == *x);
  }
  if (mode & 0x100) {
    // check binary import for BInt
    BInt val2;
    val2.from_binary(bytes, 64);
    CHECK(val == val2);
  }
  if (mode & 0x200) {
    // check binary export for BInt
    unsigned char bytes2[64];
    mval.to_binary(bytes2, 64);
    CHECK(!memcmp(bytes, bytes2, 64));
  }
  // check sign
  int sgn = mval.sgn();
  CHECK(x->sgn() == sgn);
  // check if small (fits into 64 bits)
  bool is_small = mr_is_small(mval);
  CHECK(is_small == x->fits_bits(64));
  if (is_small) {
    // special check for small (64-bit) values
    long long xval = (long long)mval;
    CHECK(x->to_long() == xval);
    CHECK((long long)__builtin_bswap64(*(long long *)(bytes + 64 - 8)) == xval);
    // check comparison with long long
    CHECK(x == xval);
    CHECK(!cmp(x, xval));
    // check constructor from long long
    CHECK(!cmp(x, td::make_refint(xval)));
    if (xval != ll_min) {
      CHECK(x > xval - 1);
      CHECK(x > td::make_refint(xval - 1));
    }
    if (xval != ll_max) {
      CHECK(x < xval + 1);
      CHECK(x < td::make_refint(xval + 1));
    }
  }
  if (mode & 0x10) {
    // check decimal export
    std::string dec = mval.to_dec_string();
    CHECK(x->to_dec_string() == dec);
    // check decimal import
    td::BigInt256 y;
    int l = y.parse_dec(dec);
    CHECK((std::size_t)l == dec.size() && y == *x);
    if (mode & 0x1000) {
      // check decimal import for BInt
      BInt val2;
      CHECK(val2.from_dec_string(dec) && val2 == val);
    }
  }
  if (mode & 0x20) {
    // check binary bit size
    int sz = x->bit_size();
    CHECK(sz >= 0 && sz <= 300);
    CHECK(x->fits_bits(sz) && (!sz || !x->fits_bits(sz - 1)));
    CHECK(mr_fits_bits(mval, sz) && !mr_fits_bits(mval, sz - 1));
    int usz = x->bit_size(false);
    CHECK(sgn >= 0 || usz == 0x7fffffff);
    if (sgn >= 0) {
      CHECK(x->unsigned_fits_bits(usz) && (!usz || !x->unsigned_fits_bits(usz - 1)));
      CHECK(mr_ufits_bits(mval, usz) && !mr_ufits_bits(mval, usz - 1));
    } else {
      CHECK(!x->unsigned_fits_bits(256) && !x->unsigned_fits_bits(300));
    }
  }
}

void check_special_int(int idx) {
  BInt b;
  unsigned char binary[64];
  td::RefInt256 x = make_special_int(idx, &b, binary);
  check_one_int_repr(x, mkint_chk_mode, idx >= min_spec_int && idx <= max_spec_int, &b, binary);
}

void init_aux() {
  np2_256 = p2_256 = MRInt::pow2(256);
  np2_256.negate();
  CHECK(np2_256 == MRInt::negpow2(256));
  p2_63 = np2_63 = MRInt::pow2(63);
  np2_63.negate();
  CHECK(np2_63 == MRInt::negpow2(63));
}

std::vector<td::RefInt256> SpecInt;
BInt SpecIntB[max_spec_int - min_spec_int + 1];

void init_check_special_ints() {
  std::cerr << "check_special_ints" << std::endl;
  BInt b;
  unsigned char binary[64];
  for (int idx = min_spec_int - 512; idx <= max_spec_int + 512; idx++) {
    td::RefInt256 x = make_special_int(idx, &b, binary);
    check_one_int_repr(x, mkint_chk_mode, idx >= min_spec_int && idx <= max_spec_int, &b, binary);
    if (idx >= min_spec_int && idx <= max_spec_int) {
      SpecIntB[idx - min_spec_int] = b;
      SpecInt.push_back(std::move(x));
    }
  }
}

void check_res(td::RefInt256 y, const BInt& yv) {
  check_one_int_repr(std::move(y), res_chk_mode, -2, &yv);
}

void check_unary_ops_on(td::RefInt256 x, const BInt& xv) {
  // NEGATE
  BInt yv = -xv;
  check_res(-x, yv);
  // NOT 
  check_res(~x, yv -= 1);
}

void check_unary_ops() {
  std::cerr << "check unary ops" << std::endl;
  for (int idx = min_spec_int; idx <= max_spec_int; idx++) {
    check_unary_ops_on(SpecInt[idx - min_spec_int], SpecIntB[idx - min_spec_int]);
  }
}

void check_pow2_ops(int shift) {
  // POW2
  td::RefInt256 r{true};
  r.unique_write().set_pow2(shift);
  check_res(r, BInt::pow2(shift));
  // POW2DEC
  r.unique_write().set_pow2(shift).add_tiny(-1).normalize();
  check_res(r, BInt::pow2(shift) - 1);
  // NEGPOW2
  r.unique_write().set_pow2(shift).negate().normalize();
  check_res(r, -BInt::pow2(shift));
}

void check_pow2_ops() {
  std::cerr << "check power-2 ops" << std::endl;
  for (int i = 0; i <= 256; i++) {
    check_pow2_ops(i);
  }
}

void check_shift_ops_on(int shift, td::RefInt256 x, const BInt& xv) {
  // LSHIFT
  check_res(x << shift, xv << shift);
  // FITS
  MRInt mval(xv);
  CHECK(x->fits_bits(shift) == mr_fits_bits(mval, shift));
  // UFITS
  CHECK(x->unsigned_fits_bits(shift) == mr_ufits_bits(mval, shift));
  // ADDPOW2 / SUBPOW2
  auto y = x;
  y.write().add_pow2(shift).normalize();
  check_res(std::move(y), xv + BInt::pow2(shift));
  y = x;
  y.write().sub_pow2(shift).normalize();
  check_res(std::move(y), xv - BInt::pow2(shift));
  // RSHIFT, MODPOW2
  for (int round_mode = -1; round_mode <= 1; round_mode++) {
    auto r = x, q = td::rshift(x, shift, round_mode);   // RSHIFT
    CHECK(q.not_null() && q->is_valid());
    r.write().mod_pow2(shift, round_mode).normalize();  // MODPOW2
    CHECK(r.not_null() && r->is_valid());
    if (round_mode < 0) {
      CHECK(!cmp(x >> shift, q));  // operator>> should be equivalent to td::rshift
    }
    BInt qv = extract_value(*q), rv = extract_value(*r);
    // check main division equality (q << shift) + r == x
    CHECK((qv << shift) + rv == xv);
    MRInt rval(rv);
    // check remainder range
    switch (round_mode) {
    case 1:
      rval.negate();  // fallthrough
    case -1:
      CHECK(mr_ufits_bits(rval, shift));
      break;
    case 0:
      CHECK(mr_fits_bits(rval, shift));
    }
  }
}

void check_shift_ops() {
  std::cerr << "check left/right shift ops" << std::endl;
  for (int idx = min_spec_int; idx <= max_spec_int; idx++) {
  //for (int idx = -52239; idx <= max_spec_int; idx++) {
    std::cerr << "idx=" << idx << " : " << SpecIntB[idx - min_spec_int] << std::endl;
    for (int i = 0; i <= 256; i++) {
      check_shift_ops_on(i, SpecInt[idx - min_spec_int], SpecIntB[idx - min_spec_int]);
    }
  }
}

int main() {
  modint::init();
  init_aux();
  init_check_special_ints();
  check_pow2_ops();
  check_unary_ops();
  check_shift_ops();
  return 0;
}
