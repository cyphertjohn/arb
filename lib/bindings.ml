(**This module contains a functor which when provided with a stub.ml file will contatin all of the bindings to the arb libraries. 
   It is not meant to be used by users. *)

open Ctypes

module Bindings (F : Cstubs.FOREIGN) = struct
  open F
  
  (*Types and functions from gmp.h*)
  let mp_limb_t = ptr uint
  type mpz_struct
  let mpz_struct : mpz_struct structure typ = structure "__mpz_struct"
  let mp_alloc = field mpz_struct "_mp_alloc" int
  let mp_size = field mpz_struct "_mp_size" int
  let mp_d = field mpz_struct "_mp_d" (ptr mp_limb_t)
  let () = seal mpz_struct
  type mpz_t = mpz_struct structure ptr
  let mpz_t : mpz_t typ = ptr mpz_struct


  (*Types and functions from fmpz.h*)
  let fmpz = typedef long "fmpz"
  type fmpz_t
  let fmpz_t = ptr fmpz
   (*typedef (array 1 fmpz) "fmpz_t"*)
  (* For whatever reason Ctypes doesn't like the fmpz_t type. So we translate to a void ptr, 
     and let c handle it.*)
  let fmpz_init = foreign "fmpz_init" (fmpz_t @-> returning void)
  let fmpz_clear = foreign "fmpz_clear" (fmpz_t @-> returning void)
  let fmpz_get_mpz = foreign "fmpz_get_mpz" (mpz_t @-> fmpz_t @-> returning void)
  let fmpz_set_mpz = foreign "fmpz_set_mpz" (fmpz_t @-> mpz_t @-> returning void)
  let fmpz_get_str = foreign "fmpz_get_str" (ptr char @-> int @-> fmpz_t @-> returning string)
  let fmpz_set_str = foreign "fmpz_set_str" (fmpz_t @-> string @-> int @-> returning int)
  let fmpz_set = foreign "fmpz_set" (fmpz_t @-> fmpz_t @-> returning void)
  let fmpz_mul_2exp = foreign "fmpz_mul_2exp" (fmpz_t @-> fmpz_t @-> ulong @-> returning void)

  (*Types and functions from fmpq.h*)
  type fmpq
  let fmpq : fmpq structure typ = structure "fmpq"
  let num = field fmpq "num" fmpz
  let den = field fmpq "den" fmpz
  let () = seal fmpq
  type fmpq_t = fmpq structure ptr
  let fmpq_t : fmpq_t typ = ptr fmpq
  let fmpq_init = foreign "fmpq_init" (fmpq_t @-> returning void)
  let fmpq_clear = foreign "fmpq_clear" (fmpq_t @-> returning void)
  let fmpq_set_si = foreign "fmpq_set_si" (fmpq_t @-> long @-> ulong @-> returning void)



  (*Types and functions from arf.h*)
  let mp_size_t = typedef long "mp_size_t"
  let mp_limb_t = typedef ulong "mp_limp_t"
  let mp_ptr = typedef (ptr mp_limb_t) "mp_ptr"
  type mantissa_noptr_struct
  let mantissa_noptr_struct : mantissa_noptr_struct structure typ = structure "mantissa_noptr_struct"
  let d_mns = field mantissa_noptr_struct "d" (array 2 mp_limb_t)
  let () = seal mantissa_noptr_struct
  type mantissa_ptr_struct
  let mantissa_ptr_struct : mantissa_ptr_struct structure typ = structure "mantissa_ptr_struct"
  let alloc = field mantissa_ptr_struct "alloc" mp_size_t
  let d_mps = field mantissa_ptr_struct "d" mp_ptr
  let () = seal mantissa_ptr_struct
  type mantissa_struct
  let mantissa_struct :mantissa_struct union typ = union "mantissa_struct"
  let noptr = field mantissa_struct "noptr" mantissa_noptr_struct
  let ptr_ms = field mantissa_struct "ptr" mantissa_ptr_struct
  let () = seal mantissa_struct
  type arf_struct
  let arf_struct : arf_struct structure typ = structure "arf_struct"
  let exp_as = field arf_struct "exp" fmpz
  let size_as = field arf_struct "size" mp_size_t
  let d_as = field arf_struct "d" mantissa_struct
  let () = seal arf_struct
  let arf_t = typedef (array 1 arf_struct) "arf_t"
  type arf_ptr = arf_struct structure ptr
  let arf_ptr : arf_ptr typ = ptr arf_struct
  let arf_srcptr = ptr arf_struct
  let arf_init = foreign "arf_init" (arf_ptr @-> returning void)
  let arf_clear = foreign "arf_clear" (arf_ptr @-> returning void)
  let arf_get_fmpz_2exp = foreign "arf_get_fmpz_2exp" (fmpz_t @-> fmpz_t @-> arf_srcptr @-> returning void)


  (*Types and functions from mag.h*)
  type mag_struct
  let mag_struct : mag_struct structure typ = structure "mag_struct"
  type mag_t = mag_struct structure ptr
  let mag_t : mag_t typ = ptr mag_struct
  let exp_mags = field mag_struct "exp" fmpz
  let man_mags = field mag_struct "man" mp_limb_t
  let () = seal mag_struct
  let mag_init = foreign "mag_init" (mag_t @-> returning void)
  let mag_clear = foreign "mag_clear" (mag_t @-> returning void)
  let mag_get_fmpz = foreign "mag_get_fmpz" (fmpz_t @-> mag_t @-> returning void)

  (*Types and functions from arb.h*)
  type arb_struct
  let arb_struct : arb_struct structure typ = structure "arb_struct"
  let mid_arbs = field arb_struct "mid" arf_struct
  let rad_arbs = field arb_struct "rad" mag_struct
  let () = seal arb_struct
  let arb_t = typedef (array 1 arb_struct) "arb_t"
  type arb_ptr = arb_struct structure ptr
  let arb_ptr : arb_ptr typ = ptr arb_struct
  let arb_srcptr = ptr arb_struct
  let arb_init = foreign "arb_init" (arb_ptr @-> returning void)
  let arb_clear = foreign "arb_clear" (arb_ptr @-> returning void)
  let arb_get_interval_fmpz_2exp = foreign "arb_get_interval_fmpz_2exp" (fmpz_t @-> fmpz_t @-> fmpz_t @-> arb_ptr @-> returning void)
  (*let arb_get_mid_arb = foreign "arb_get_mid_arb" (arb_ptr @-> arb_srcptr @-> returning void)*)
  let arb_mid_ptr = foreign "arb_mid_ptr" (arb_ptr @-> returning arf_ptr)
  let arb_get_mag = foreign "arb_get_mag" (mag_t @-> arb_srcptr @-> returning void)


  (*Types and functions from acb.h*)
  type acb_struct
  let acb_struct : acb_struct structure typ = structure "acb_struct"
  let real = field acb_struct "real" arb_struct
  let imag = field acb_struct "imag" arb_struct
  let () = seal acb_struct
  let acb_t = typedef (array 1 acb_struct) "acb_t"
  type acb_ptr = acb_struct structure ptr
  let acb_ptr : acb_ptr typ = ptr acb_struct
  let acb_srcptr = ptr acb_struct
  let acb_init = foreign "acb_init" (acb_ptr @-> returning void)
  let acb_clear = foreign "acb_clear" (acb_ptr @-> returning void)
  let acb_get_real = foreign "acb_get_real" (arb_ptr @-> acb_srcptr @-> returning void)
  let acb_get_imag = foreign "acb_get_imag" (arb_ptr @-> acb_srcptr @-> returning void)
  let acb_set_fmpz_fmpz = foreign "acb_set_fmpz_fmpz" (acb_ptr @-> fmpz_t @-> fmpz_t @-> returning void)
  let acb_set_fmpq = foreign "acb_set_fmpq" (acb_ptr @-> fmpq_t @-> long @-> returning void)
  let acb_log = foreign "acb_log" (acb_ptr @-> acb_srcptr @-> long @-> returning void)
  let acb_const_pi = foreign "acb_const_pi" (acb_ptr @-> long @-> returning void)
  let acb_add = foreign "acb_add" (acb_ptr @-> acb_srcptr @-> acb_srcptr @-> long @-> returning void)
  let acb_sub = foreign "acb_sub" (acb_ptr @-> acb_srcptr @-> acb_srcptr @-> long @-> returning void)
  let acb_mul = foreign "acb_mul" (acb_ptr @-> acb_srcptr @-> acb_srcptr @-> long @-> returning void)
  let acb_div = foreign "acb_div" (acb_ptr @-> acb_srcptr @-> acb_srcptr @-> long @-> returning void)
  let acb_mul_si = foreign "acb_mul_si" (acb_ptr @-> acb_srcptr @-> long @-> long @-> returning void)
  let acb_div_si = foreign "acb_div_si" (acb_ptr @-> acb_srcptr @-> long @-> long @-> returning void)
  let _acb_vec_init = foreign "_acb_vec_init" (long @-> returning acb_ptr)
  let _acb_vec_clear = foreign "_acb_vec_clear" (acb_ptr @-> long @-> returning void)
  let acb_bits = foreign "acb_bits" (acb_ptr @-> returning long)
  let acb_trim = foreign "acb_trim" (acb_ptr @-> acb_srcptr @-> returning void)
  let acb_rel_error_bits = foreign "acb_rel_error_bits" (acb_srcptr @-> returning long)
  let acb_rel_accuracy_bits = foreign "acb_rel_accuracy_bits" (acb_srcptr @-> returning long)
  let acb_rel_one_accuracy_bits = foreign "acb_rel_one_accuracy_bits" (acb_srcptr @-> returning long)
  let acb_neg = foreign "acb_neg" (acb_ptr @-> acb_srcptr @-> returning void)
  let acb_pow_si = foreign "acb_pow_si" (acb_ptr @-> acb_srcptr @-> long @-> long @-> returning void)
  let acb_set_fmpz = foreign "acb_set_fmpz" (acb_ptr @-> fmpz_t @-> returning void)
  let acb_print = foreign "acb_print" (acb_ptr @-> returning void)


  (*Types and functions from flint/fmpz_poly.h*)
  type fmpz_poly_struct
  let fmpz_poly_struct : fmpz_poly_struct structure typ = structure "fmpz_poly_struct"
  let coefs = field fmpz_poly_struct "coeffs" fmpz_t
  let alloc_fmpz_poly = field fmpz_poly_struct "alloc" long
  let length = field fmpz_poly_struct "length" long
  let () = seal fmpz_poly_struct
  type fmpz_poly_t = fmpz_poly_struct structure ptr
  let fmpz_poly_t : fmpz_poly_t typ = ptr fmpz_poly_struct
  let fmpz_poly_init = foreign "fmpz_poly_init" (fmpz_poly_t @-> returning void)
  let fmpz_poly_clear = foreign "fmpz_poly_clear" (fmpz_poly_t @-> returning void)
  let fmpz_poly_set_coeff_fmpz = foreign "fmpz_poly_set_coeff_fmpz" (fmpz_poly_t @-> long @-> fmpz_t @-> returning void)
  let fmpz_poly_get_coeff_fmpz = foreign "fmpz_poly_get_coeff_fmpz" (fmpz_t @-> fmpz_poly_t @-> long @-> returning void)
  let fmpz_poly_degree = foreign "fmpz_poly_degree" (fmpz_poly_t @-> returning long)

  (*Function from arb_fmpz_poly.h*)
  let arb_fmpz_poly_complex_roots = foreign "arb_fmpz_poly_complex_roots" (acb_ptr @-> fmpz_poly_t @-> int @-> long @-> returning void)


  (*Types and functions from fmpz_mat.h*)
  type fmpz_mat_struct
  let fmpz_mat_struct : fmpz_mat_struct structure typ = structure "fmpz_mat_struct"
  let entries = field fmpz_mat_struct "entries" (ptr fmpz)
  let r = field fmpz_mat_struct "r" long
  let c = field fmpz_mat_struct "c" long
  let rows = field fmpz_mat_struct "rows" (ptr (ptr fmpz))
  let () = seal fmpz_mat_struct
  type fmpz_mat_t = fmpz_mat_struct structure ptr
  let fmpz_mat_t : fmpz_mat_t typ = ptr fmpz_mat_struct
  let fmpz_mat_entry = foreign "fmpz_mat_entry" (fmpz_mat_t @-> long @-> long @-> returning (ptr fmpz))
  let fmpz_mat_init = foreign "fmpz_mat_init" (fmpz_mat_t @-> long @-> long @-> returning void)
  let fmpz_mat_clear = foreign "fmpz_mat_clear" (fmpz_mat_t @-> returning void)
  let fmpz_mat_nrows = foreign "fmpz_mat_nrows" (fmpz_mat_t @-> returning long)
  let fmpz_mat_ncols = foreign "fmpz_mat_ncols" (fmpz_mat_t @-> returning long)
  let fmpz_mat_lll_original = foreign "fmpz_mat_lll_original" (fmpz_mat_t @-> fmpq_t @-> fmpq_t @-> returning void)
  let fmpz_mat_lll_storjohann = foreign "fmpz_mat_lll_storjohann" (fmpz_mat_t @-> fmpq_t @-> fmpq_t @-> returning void)
  let fmpz_mat_one = foreign "fmpz_mat_one" (fmpz_mat_t @-> returning void)


  
end
