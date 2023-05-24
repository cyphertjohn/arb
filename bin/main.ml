open Arb
open Arb_zarith


let fmpzify m rows cols = 
  let mat = Fmpz_mat.init rows cols in
  List.iteri (
    fun i row ->
      List.iteri (
        fun j ent ->
          Fmpz_mat.set_entry mat (Fmpzz.zarith_to_fmpz (Z.of_int ent)) i j
      ) row
  ) m;
  mat


let zify m = 
  let rows, cols = Fmpz_mat.nb_rows m, Fmpz_mat.nb_cols m in
  List.init rows (
    fun i ->
      List.init cols (
        fun j ->
          Fmpzz.fmpz_to_zarith (Fmpz_mat.get_entry m i j)
      )
  )

let to_string m = 
  let row_to_string row = 
    List.fold_left (
      fun s e ->
        s ^ "; " ^ (Z.to_string e)
    ) (Z.to_string (List.hd row)) (List.tl row) 
  in
  List.fold_left (
    fun str row ->
      str ^ "\n" ^ (row_to_string row)
  ) (row_to_string (List.hd m)) (List.tl m)

let print_acb a = 
  let (realm, reale) = Acb.get_real_mid_fmpz a in
  let (imagm, image) = Acb.get_imag_mid_fmpz a in
  let bits = Acb.bits a in
  let error_bits = Acb.rel_error_bits a in
  let accuracy_bits = Acb.rel_accuracy_bits a in
  let accuracy_one_bits = Acb.rel_one_accuracy_bits a in
  let str = ((Z.to_string (Fmpzz.fmpz_to_zarith realm))  ^ "*2^" ^ (Z.to_string (Fmpzz.fmpz_to_zarith reale))) 
          ^ " + i" ^ ((Z.to_string (Fmpzz.fmpz_to_zarith imagm))  ^ "*2^" ^ (Z.to_string (Fmpzz.fmpz_to_zarith image))) in
  print_endline str;
  print_endline ("Bits: " ^ (string_of_int bits));
  print_endline ("Error bits: " ^ (string_of_int error_bits));
  print_endline ("Accuracy bits: " ^ (string_of_int accuracy_bits));
  print_endline ("Accuracy one bits: " ^ (string_of_int accuracy_one_bits));
  print_endline ""



let () = 
  let one = Fmpzz.zarith_to_fmpz Z.one in
  let minus_one = Fmpzz.zarith_to_fmpz (Z.minus_one) in
  let x2_minus_x_minus_1 = Fmpz_poly.init () in (*x^2 - x - 1*)
  Fmpz_poly.set_coef x2_minus_x_minus_1 0 minus_one;
  Fmpz_poly.set_coef x2_minus_x_minus_1 1 minus_one;
  Fmpz_poly.set_coef x2_minus_x_minus_1 2 one;
  let prec = 12 in
  let roots = Fmpz_poly.get_complex_roots x2_minus_x_minus_1 prec in
  let num_roots = List.length roots in
  print_endline ("Num of roots:" ^ (string_of_int num_roots));
  let r1 = List.hd roots in
  let r2 = List.hd (List.tl roots) in 

  let acb_minus_one = Acb.init_set_fmpz minus_one in
  let log_sigma1_epsilon_1 = Acb.div (Acb.div_si (Acb.log acb_minus_one prec) 2 prec) (Acb.pi prec) prec in
  let log_sigma2_epsilon_1 = Acb.div (Acb.div_si (Acb.log acb_minus_one prec) 2 prec) (Acb.pi prec) prec in
  let log_sigma1_epsilon_2 = Acb.div (Acb.div_si (Acb.log r1 prec) 2 prec) (Acb.pi prec) prec in
  let log_sigma2_epsilon_2 = Acb.div (Acb.div_si (Acb.log r2 prec) 2 prec) (Acb.pi prec) prec in
  let log_sigma1_epsilon_3 = Acb.div (Acb.div_si (Acb.log r2 prec) 2 prec) (Acb.pi prec) prec in
  let log_sigma2_epsilon_3 = Acb.div (Acb.div_si (Acb.log r1 prec) 2 prec) (Acb.pi prec) prec in

  print_endline "log_sigma1_epsilon1";
  print_acb log_sigma1_epsilon_1;
  print_endline "log_sigma1_epsilon1";
  print_acb log_sigma2_epsilon_1;
  print_endline "log_sigma1_epsilon2";
  print_acb log_sigma1_epsilon_2;
  print_endline "log_sigma2_epsilon2";
  print_acb log_sigma2_epsilon_2;
  print_endline "log_sigma1_epsilon3";
  print_acb log_sigma1_epsilon_3;
  print_endline "log_sigma2_epsilon3";
  print_acb log_sigma2_epsilon_3;

  let close_to_2 = Acb.div (Acb.mul_si (Acb.pi prec) 10 prec) (Acb.mul_si (Acb.pi prec) 5 prec) prec in
  print_endline "Close to 2";
  print_acb close_to_2;

  let example = 
    [
      [1; 0; 0; 0; 0; 0; 16384; 0; 16384];
      [0; 1; 0; 0; 0; -2509; 16384; 2509; 0];
      [0; 0; 1; 0; 0; 2509; 0; -2509; 16384];
      [0; 0; 0; 1; 0; 0; 32768; 0; 0];
      [0; 0; 0; 0; 1; 0; 0; 0; 32768]
    ] in
  let e_z = fmpzify example 5 9 in
  Fmpz_mat.lll_storjohann e_z (75, 100) (51, 100);
  let e_zz = zify e_z in
  print_endline "Reduced mat";
  print_endline (to_string e_zz)

