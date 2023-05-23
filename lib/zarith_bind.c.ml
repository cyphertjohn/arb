(**Copied from https://github.com/fdopen/ctypes-zarith/blob/master/lib/ctypes_zarith.c.ml *)

let%c () = header {|
#include <gmp.h>
#include <zarith.h>
|}

module MPZ = struct
  let%c t = abstract "__mpz_struct"

  type ptr = t Ctypes.abstract Ctypes.ptr

  external clear : t ptr -> void = "mpz_clear" [@@noalloc]

  external init : t ptr -> void = "mpz_init"

  let make () =
    (* allocate_n zero initializes the memory. It's safe to pass
       such a struct to mpz_clear. *)
    let r = Ctypes.allocate_n ~finalise:clear t ~count:1 in
    init r;
    r

  [%%c
  external set : zt_:(Z.t[@ocaml_type]) -> tptr_:t ptr -> void
    = {|
   value z = $zt_; /* not converted. The usual rules for stub code must be
                      obeyed (accessors, memory management (GC), etc.) */
   __mpz_struct * p = $tptr_; /* already converted to a native c type, will not
                                 be garbage collected during the stub code */
   ml_z_mpz_set_z(p, z);
|}]

  let of_z x =
    let r = make () in
    set x r;
    r

  [%%c
  external to_z : ptr_:t ptr -> (Z.t[@ocaml_type])
    = {|
  return (ml_z_from_mpz($ptr_));
|}]

end
