let () = 
  let tmp_out = open_out "arb_stub.c" in
  let ml_out = open_out "arb_stub.ml" in
  output_string ml_out "(**Generated stub file. Not meant to be used by users of the arb package. *)";
  Cstubs.write_ml (Format.formatter_of_out_channel ml_out) ~prefix:"arb_stub" (module Bindings.Bindings);
  output_string tmp_out "#include \"flint/fmpz_mat.h\"\n#include \"flint/fmpz_poly.h\"\n#include \"acb.h\"\n";
  Cstubs.write_c (Format.formatter_of_out_channel tmp_out) ~prefix:"arb_stub" (module Bindings.Bindings);
  flush tmp_out;
  flush ml_out;
  close_out tmp_out;
  close_out ml_out
