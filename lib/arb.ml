module C = Bind.Bindings(Arb_stub)


module Fmpz = struct
  type t = Signed.long Ctypes.ptr

  let init () = 
    let a = Ctypes.CArray.start (Ctypes.CArray.make C.fmpz 1) in
    C.fmpz_init a;
    a

  let clear a = 
    C.fmpz_clear a

  let get_str a =
    let buff = Ctypes.from_voidp Ctypes.char Ctypes.null in
    C.fmpz_get_str buff 10 a

  let set_str a str radix = 
    let _ = C.fmpz_set_str a str radix in
    ()

  let init_set_str str radix = 
    let a = init () in
    let _ = set_str a str radix in
    a

  let fmpz_from_mpz mpz = 
    let r = Ctypes.CArray.start (Ctypes.CArray.make C.fmpz 1) in (* I don't think it needs to be init*)
    C.fmpz_set_mpz r mpz;
    r

  let mpz_from_fmpz fmpz = 
    let r = Ctypes.CArray.start (Ctypes.CArray.make C.mpz_struct 1) in (* I don't think it needs to be init*)
    C.fmpz_get_mpz r fmpz;
    r

  let set dest source = C.fmpz_set dest source

  let mul_2exp a e = 
    let r = init () in
    C.fmpz_mul_2exp r a (Unsigned.ULong.of_int e);
    r

  let init_set a = 
    let res = Ctypes.CArray.start (Ctypes.CArray.make C.fmpz 1) in
    C.fmpz_init_set res a;
    res

  let of_int i = 
    let r = Ctypes.CArray.start (Ctypes.CArray.make C.fmpz 1) in
    C.fmpz_init_set_si r (Signed.Long.of_int i);
    r

  let to_int z = 
    Signed.Long.to_int (C.fmpz_get_si z)

  let zero () =
    let r = init () in
    C.fmpz_zero r;
    r
  
  let one () = 
    let r = init () in
    C.fmpz_one r;
    r

  let cmp a b =
    C.fmpz_cmp a b

  let cmp_si a i = 
    C.fmpz_cmp_si a (Signed.Long.of_int i)

  let equal a b = 
    if C.fmpz_equal a b = 1 then true else false

  let equal_si a i = 
    if C.fmpz_equal_si a (Signed.Long.of_int i) = 1 then true else false

  let neg a =
    let r = init () in
    C.fmpz_neg r a;
    r

  let add a b = 
    let r = init () in 
    C.fmpz_add r a b;
    r

  let add_si a i = 
    let r = init () in
    C.fmpz_add_si r a (Signed.Long.of_int i);
    r

  let sub a b = 
    let r = init () in
    C.fmpz_sub r a b;
    r

  let sub_si a i = 
    let r = init () in
    C.fmpz_sub_si r a (Signed.Long.of_int i);
    r

  let mul a b = 
    let r = init () in
    C.fmpz_mul r a b;
    r

  let mul_si a i = 
    let r = init () in
    C.fmpz_mul_si r a (Signed.Long.of_int i);
    r

  let divexact a b = 
    let r = init () in 
    C.fmpz_divexact r a b;
    r

  let divexact_si a i = 
    let r = init () in 
    C.fmpz_divexact_si r a (Signed.Long.of_int i);
    r

  let pow_ui a i = 
    if i < 0 then raise (Invalid_argument "Negative Exponent");
    let r = init () in
    C.fmpz_pow_ui r a (Unsigned.ULong.of_int i);
    r

  let gcd a b = 
    let r = init () in
    C.fmpz_gcd r a b;
    r

  let lcm a b = 
    let r = init () in
    C.fmpz_lcm r a b;
    r


end


module Acb = struct
  type t = C.acb_ptr

  let init () = 
    let acb_struct = Ctypes.make C.acb_struct in
    let c = Ctypes.addr acb_struct in
    C.acb_init c;
    c

  let clear a = C.acb_clear a

  let set_fmpz = C.acb_set_fmpz

  let init_set_fmpz real = 
    let c = init () in
    set_fmpz c real;
    c

  let set_fmpz_fmpz c real imag = 
    C.acb_set_fmpz_fmpz c real imag

  let init_set_fmpz_fmpz real imag = 
    let c = init () in
    set_fmpz_fmpz c real imag;
    c
    

  let log a prec = 
    let r = init () in
    C.acb_log r a (Signed.Long.of_int prec);
    r

  let get_real_mid_fmpz a = 
    let real = Ctypes.CArray.start (Ctypes.CArray.make C.arb_struct 1) in
    (*C.arb_init mid;*)
    C.acb_get_real real a;
    let m = Fmpz.init () in
    let e = Fmpz.init () in
    C.arf_get_fmpz_2exp m e (C.arb_mid_ptr real);
    (m, e)
    
  let get_imag_mid_fmpz a = 
    let imag = Ctypes.CArray.start (Ctypes.CArray.make C.arb_struct 1) in
    (*C.arb_init mid;*)
    C.acb_get_imag imag a;
    let m = Fmpz.init () in
    let e = Fmpz.init () in
    C.arf_get_fmpz_2exp m e (C.arb_mid_ptr imag);
    (m, e)

  let get_real_imag_mag_upper a = 
    let real = Ctypes.addr (Ctypes.make C.arb_struct) in
    C.acb_get_real real a;
    let imag = Ctypes.addr (Ctypes.make C.arb_struct) in
    C.acb_get_imag imag a;
    let real_mag = Ctypes.addr (Ctypes.make C.mag_struct) in
    C.arb_get_mag real_mag real;
    let imag_mag = Ctypes.addr (Ctypes.make C.mag_struct) in
    C.arb_get_mag imag_mag imag;
    let real_upper = Fmpz.init () in
    C.mag_get_fmpz real_upper real_mag;
    let imag_upper = Fmpz.init () in
    C.mag_get_fmpz imag_upper imag_mag;
    (real_upper, imag_upper)



  let pi prec = 
    let p = init () in
    C.acb_const_pi p (Signed.Long.of_int prec);
    p
  
  let add a b prec = 
    let r = init () in
    C.acb_add r a b (Signed.Long.of_int prec);
    r
  
  let sub a b prec = 
    let r = init () in
    C.acb_sub r a b (Signed.Long.of_int prec);
    r

  let mul a b prec = 
    let r = init () in
    C.acb_mul r a b (Signed.Long.of_int prec);
    r

  let div a b prec = 
    let r = init () in
    C.acb_div r a b (Signed.Long.of_int prec);
    r

  let mul_si a i prec = 
    let r = init () in
    C.acb_mul_si r a (Signed.Long.of_int i) (Signed.Long.of_int prec);
    r

  let div_si a i prec = 
    let r = init () in
    C.acb_div_si r a (Signed.Long.of_int i) (Signed.Long.of_int prec);
    r

  let bits a = Signed.Long.to_int (C.acb_bits a)

  let trim a = 
    let r = init () in
    C.acb_trim r a;
    r

  let rel_error_bits a = Signed.Long.to_int (C.acb_rel_error_bits a)

  let rel_accuracy_bits a = Signed.Long.to_int (C.acb_rel_accuracy_bits a)

  let rel_one_accuracy_bits a = Signed.Long.to_int (C.acb_rel_one_accuracy_bits a)

  let neg a = 
    let r = init () in
    C.acb_neg r a;
    r

  let pow_si a e prec = 
    let r = init () in
    C.acb_pow_si r a (Signed.Long.of_int e) (Signed.Long.of_int prec);
    r
end


module Fmpz_poly = struct
  type t = C.fmpz_poly_t

  let init () = 
    let poly_struct = Ctypes.make C.fmpz_poly_struct in
    let c = Ctypes.addr poly_struct in
    C.fmpz_poly_init c;
    c
  
  let clear a = 
    C.fmpz_poly_clear a

  let set_coef p deg coef = 
    C.fmpz_poly_set_coeff_fmpz p (Signed.Long.of_int deg) coef

  let get_coef p deg = 
    let coef = Fmpz.init () in
    C.fmpz_poly_get_coeff_fmpz coef p (Signed.Long.of_int deg);
    coef

  let get_complex_roots p prec = 
    let d = C.fmpz_poly_degree p in
    let di = Signed.Long.to_int d in
    let roots = C._acb_vec_init d in
    C.arb_fmpz_poly_complex_roots roots p 0 (Signed.Long.of_int prec);
    List.map (
      fun index -> 
        Ctypes.(+@) roots index
    ) (List.init di (fun i -> i))
end

module Fmpz_mat = struct

  type t = C.fmpz_mat_t

  exception Incompatible_Dimensions
  exception Index_Out_of_Bounds

  let init nrows ncols : t = 
    if nrows < 0 || ncols < 0 then raise Index_Out_of_Bounds;
    let mat_struct = Ctypes.make C.fmpz_mat_struct in
    let c = Ctypes.addr mat_struct in
    C.fmpz_mat_init c (Signed.Long.of_int nrows) (Signed.Long.of_int ncols);
    c

  let copy a = 
    let mat_struct = Ctypes.make C.fmpz_mat_struct in
    let c = Ctypes.addr mat_struct in
    C.fmpz_mat_init_set c a;
    c

  let clear a = 
    C.fmpz_mat_clear a
  
  let nb_rows mat = 
    Signed.Long.to_int (C.fmpz_mat_nrows mat)

  let nb_cols mat = 
    Signed.Long.to_int (C.fmpz_mat_ncols mat)

  
  let set_entry mat valu i j = 
    let m, n = nb_rows mat, nb_cols mat in
    if i >= m || j >= n then raise Index_Out_of_Bounds;
    let entry = C.fmpz_mat_entry mat (Signed.Long.of_int i) (Signed.Long.of_int j) in
    Fmpz.set entry valu
  
  let get_entry mat i j = 
    let m, n = nb_rows mat, nb_cols mat in
    if i >= m || j >= n then raise Index_Out_of_Bounds;
    let res = Fmpz.init () in
    let entry = C.fmpz_mat_entry mat (Signed.Long.of_int i) (Signed.Long.of_int j) in
    Fmpz.set res entry;
    res

  let fmpq_init_set_si n d = 
    let fmpq = Ctypes.make C.fmpq in
    let res = Ctypes.addr fmpq in
    C.fmpq_init res;
    C.fmpq_set_si res (Signed.Long.of_int n) (Unsigned.ULong.of_int d);
    res

  let lll_original mat (delta_n, delta_d) (eta_n, eta_d) = 
    let delta = fmpq_init_set_si delta_n delta_d in
    let eta = fmpq_init_set_si eta_n eta_d in
    C.fmpz_mat_lll_original mat delta eta

  let lll_storjohann mat (delta_n, delta_d) (eta_n, eta_d) = 
    let delta = fmpq_init_set_si delta_n delta_d in
    let eta = fmpq_init_set_si eta_n eta_d in
    C.fmpz_mat_lll_storjohann mat delta eta

  let ident m n = 
    let i = init m n in
    C.fmpz_mat_one i;
    i

  let zero m n = 
    let i = init m n in
    C.fmpz_mat_zero i;
    i

  let window mat r1 c1 r2 c2 = 
    let m, n = nb_rows mat, nb_cols mat in
    if r1 < 0 || r2 < 0 || c1 < 0 || c2 < 0 || r1+r2 > m || c1+c2 > n then raise Index_Out_of_Bounds;
    let mat_struct = Ctypes.make C.fmpz_mat_struct in
    let res = Ctypes.addr mat_struct in
    C.fmpz_mat_window_init res mat (Signed.Long.of_int r1) (Signed.Long.of_int c1) (Signed.Long.of_int r2) (Signed.Long.of_int c2);
    copy res


  let equal a b =
    if (C.fmpz_mat_equal a b) = 0 then false
    else true

  let is_zero mat = 
    if (C.fmpz_mat_is_zero mat) = 0 then false
    else true

  let transpose mat = 
    let m, n = nb_rows mat, nb_cols mat in
    let res = init n m in
    C.fmpz_mat_transpose res mat;
    res
  

  let concat_vertical mat1 mat2 = 
    let n1, n2 = nb_cols mat1, nb_cols mat2 in
    if n1 <> n2 then raise Incompatible_Dimensions;
    let res = init ((nb_rows mat1) + (nb_rows mat2)) n1 in
    C.fmpz_mat_concat_vertical res mat1 mat2;
    res

  let concat_horizontal mat1 mat2 = 
    let m1, m2 = nb_rows mat1, nb_rows mat2 in
    if m1 <> m2 then raise Incompatible_Dimensions;
    let res = init m1 ((nb_cols mat1) + (nb_cols mat2)) in
    C.fmpz_mat_concat_horizontal res mat1 mat2;
    res

  let add mat1 mat2 = 
    let (m1, n1) = nb_rows mat1, nb_cols mat1 in
    let (m2, n2) = nb_rows mat2, nb_cols mat2 in
    if m1 <> m2 || n1 <> n2 then raise Incompatible_Dimensions;
    let res = init m1 n1 in
    C.fmpz_mat_add res mat1 mat2;
    res

  let sub mat1 mat2 = 
    let (m1, n1) = nb_rows mat1, nb_cols mat1 in
    let (m2, n2) = nb_rows mat2, nb_cols mat2 in
    if m1 <> m2 || n1 <> n2 then raise Incompatible_Dimensions;
    let res = init m1 n1 in
    C.fmpz_mat_sub res mat1 mat2;
    res

  let neg a = 
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_neg res a;
    res

  let scalar_mult_si a i = 
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_scalar_mul_si res a (Signed.Long.of_int i);
    res
  
  let scalar_mult a z = 
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_scalar_mul_fmpz res a z;
    res

  let divexact_si a i = 
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_scalar_divexact_si res a (Signed.Long.of_int i);
    res

  let divexact a z = 
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_scalar_divexact_fmpz res a z;
    res

  let scalar_mult_2exp a e = 
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_scalar_mul_2exp res a (Unsigned.ULong.of_int e);
    res

  let scalar_tdiv_q_2exp a e = 
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_scalar_tdiv_q_2exp res a (Unsigned.ULong.of_int e);
    res

  let mul a b = 
    let (m1, n1) = nb_rows a, nb_cols a in
    let (m2, n2) = nb_rows b, nb_cols b in
    if n1 <> m2 then raise Incompatible_Dimensions;
    let res = init m1 n2 in
    C.fmpz_mat_mul res a b;
    res

  let kronecker a b = 
    let m1, n1 = nb_rows a, nb_cols a in
    let m2, n2 = nb_rows b, nb_cols b in
    let res = init (m1 * m2) (n1 * n2) in
    for i = 0 to m1 * m2 - 1 do
      for j = 0 to n1 * n2 - 1 do
        let entry = C.fmpz_mat_entry res (Signed.Long.of_int i) (Signed.Long.of_int j) in
        let aentry = C.fmpz_mat_entry res (Signed.Long.of_int (i / m2)) (Signed.Long.of_int (j / n2)) in
        let bentry = C.fmpz_mat_entry res (Signed.Long.of_int (i mod m2)) (Signed.Long.of_int (j mod n2)) in
        C.fmpz_init entry;
        C.fmpz_mul entry aentry bentry
      done;
    done;
    res

  let pow a e = 
    if e < 0 then raise (Invalid_argument "Negative exponent");
    let res = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_pow res a (Unsigned.ULong.of_int e);
    res

  let content a = 
    let c = Fmpz.init () in
    C.fmpz_mat_content c a;
    c

  let trace a = 
    let c = Fmpz.init () in
    C.fmpz_mat_trace c a;
    c

  let hnf a = 
    let h = init (nb_rows a) (nb_cols a) in
    C.fmpz_mat_hnf h a;
    h

  let hnf_transform a = 
    let m, n = (nb_rows a), (nb_cols a) in
    let h = init m n in
    let u = init m m in
    C.fmpz_mat_hnf_transform h u a;
    (h, u)

  let is_in_hnf a = 
    if C.fmpz_mat_is_in_hnf a = 1 then true else false
  

end