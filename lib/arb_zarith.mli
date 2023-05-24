(** Conversion from zarith to Fmpz.t*)
module Fmpzz : sig
  
  (** Convert a zarith integer to an fmpz.*)
  val zarith_to_fmpz : Z.t -> Arb.Fmpz.t

  (** Convert an fmpz value to a zarith integer.*)
  val fmpz_to_zarith : Arb.Fmpz.t -> Z.t
  
end
