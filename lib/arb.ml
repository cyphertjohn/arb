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

  let init nrows ncols : t = 
    let mat_struct = Ctypes.make C.fmpz_mat_struct in
    let c = Ctypes.addr mat_struct in
    (*let c = Ctypes.CArray.start (Ctypes.CArray.make ~finalise:(fun _ -> print_endline "Freeing matrix") C.fmpz_mat_struct 1) in*)
    C.fmpz_mat_init c (Signed.Long.of_int nrows) (Signed.Long.of_int ncols);
    (*for i = 0 to nrows - 1 do
      for j = 0 to ncols - 1 do
        C.fmpz_init (C.fmpz_mat_entry c (Signed.Long.of_int i) (Signed.Long.of_int j))
      done
    done;*)
    c

  let clear a = 
    C.fmpz_mat_clear a

  let set_entry mat valu i j = 
    let entry = C.fmpz_mat_entry mat (Signed.Long.of_int i) (Signed.Long.of_int j) in
    Fmpz.set entry valu

  let get_entry mat i j = 
    let res = Fmpz.init () in
    let entry = C.fmpz_mat_entry mat (Signed.Long.of_int i) (Signed.Long.of_int j) in
    Fmpz.set res entry;
    res
  
  let nb_rows mat = 
    Signed.Long.to_int (C.fmpz_mat_nrows mat)

  let nb_cols mat = 
    Signed.Long.to_int (C.fmpz_mat_ncols mat)

  let fmpq_init_set_si n d = 
    let fmpq = Ctypes.make ~finalise:(fun _ -> print_endline "Freeing fmpq struct"; flush stdout) C.fmpq in
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

end