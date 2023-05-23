(** The zbind library provides a conversion from Zarith.Z.t to a Ctypes pointer (essentially a *void). At the level of C this resulting pointer can be cast as an mpz_t type.
    For the arb package this pointer is passed to the arb code as an mpz_t value. The actual conversion from Zarith.Z.t to *void is done via a function
    distributed with the zarith package found in zarith.h.
    
    This code is not meant to be used by users of the faugere package.*)

module MPZ : sig

  type t

  type ptr = t Ctypes.abstract Ctypes.ptr

  val t : t Ctypes.abstract Ctypes.typ

  val of_z : Z.t -> ptr

  val to_z : ptr -> Z.t

end
