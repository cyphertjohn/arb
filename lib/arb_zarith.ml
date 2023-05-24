module Fmpzz = struct
  open Arb.Fmpz

  let zarith_to_fmpz z = 
    let ptr = Zarith_bind.MPZ.of_z z in
    let ptr_mpz = Ctypes.from_voidp (Arb.C.mpz_struct) (Ctypes.to_voidp ptr) in
    fmpz_from_mpz ptr_mpz

  let fmpz_to_zarith fmpz = 
    let mpz = mpz_from_fmpz fmpz in
    Zarith_bind.MPZ.to_z (Ctypes.from_voidp Zarith_bind.MPZ.t (Ctypes.to_voidp mpz))
end

(*module Acb = Arb.Acb

module Fmpz_poly = Arb.Fmpz_poly

module Fmpz_mat = Arb.Fmpz_mat*)