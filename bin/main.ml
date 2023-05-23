open Arb_zarith


let fmpzify m rows cols = 
  let mat = Fmpz_mat.init rows cols in
  List.iteri (
    fun i row ->
      List.iteri (
        fun j ent ->
          Fmpz_mat.set_entry mat (Fmpz.zarith_to_fmpz (Z.of_int ent)) i j
      ) row
  ) m;
  mat


let zify m = 
  let rows, cols = Fmpz_mat.nb_rows m, Fmpz_mat.nb_cols m in
  List.init rows (
    fun i ->
      List.init cols (
        fun j ->
          Fmpz.fmpz_to_zarith (Fmpz_mat.get_entry m i j)
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



let () = 
  let one = Fmpz.zarith_to_fmpz Z.one in
  let x2_x_1 = Fmpz_poly.init () in (*x^2 + x + 1*)
  Fmpz_poly.set_coef x2_x_1 0 one;
  Fmpz_poly.set_coef x2_x_1 1 one;
  Fmpz_poly.set_coef x2_x_1 2 one;
  let roots = Fmpz_poly.get_complex_roots x2_x_1 10 in
  let num_roots = List.length roots in
  print_endline ("Num of roots:" ^ (string_of_int num_roots));
  let r1 = List.hd roots in (*-1/2+i sqrt(3)/2*)
  let r2 = List.hd (List.tl roots) in (*-1/2-i sqrt(3)/2*)
  let (realmr1, realer1) = Acb.get_real_mid_fmpz r1 in
  let (imagmr1, imager1) = Acb.get_imag_mid_fmpz r1 in  
  print_endline "Real midpoint root 1";
  print_endline ((Z.to_string (Fmpz.fmpz_to_zarith realmr1))  ^ "*2^" ^ (Z.to_string (Fmpz.fmpz_to_zarith realer1)));
  print_endline "Complex midpoint root 1";
  print_endline ((Z.to_string (Fmpz.fmpz_to_zarith imagmr1))  ^ "*2^" ^ (Z.to_string (Fmpz.fmpz_to_zarith imager1)));

  let (realmr2, realer2) = Acb.get_real_mid_fmpz r2 in
  let (imagmr2, imager2) = Acb.get_imag_mid_fmpz r2 in  
  print_endline "Real midpoint root 2";
  print_endline ((Z.to_string (Fmpz.fmpz_to_zarith realmr2))  ^ "*2^" ^ (Z.to_string (Fmpz.fmpz_to_zarith realer2)));
  print_endline "Complex midpoint root 2";
  print_endline ((Z.to_string (Fmpz.fmpz_to_zarith imagmr2))  ^ "*2^" ^ (Z.to_string (Fmpz.fmpz_to_zarith imager2)));


  let example = 
    [
      [1; 0; 0; 0; 0; 0; 256; 0; 256];
      [0; 1; 0; 0; 0; 40; 0; -40; 256];
      [0; 0; 1; 0; 0; -40; 256; 40; 0];
      [0; 0; 0; 1; 0; 0; 512; 0; 0];
      [0; 0; 0; 0; 1; 0; 0; 0; 512]
    ] in
  let e_z = fmpzify example 5 9 in
  Fmpz_mat.lll_storjohann e_z (75, 100) (51, 100);
  let e_zz = zify e_z in
  print_endline "Reduced mat";
  print_endline (to_string e_zz)

