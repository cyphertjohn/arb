(** Interface to some functions of the arb and flint library. For more information on the
    functionality of these functions see the arb documentation.*)

(** The type of mpz's so outside libraries can connect.*)
module C : sig
  type mpz_struct
  val mpz_struct : mpz_struct Ctypes.structure Ctypes.typ
  type mpz_t = mpz_struct Ctypes.structure Ctypes.ptr
  val mpz_t : mpz_t Ctypes.typ
end

(** Flint's version of mpz's. *)
module Fmpz : sig
  type t

  (** Get a fresh fmpz ready for use*)
  val init : unit -> t

  (** Clear and free the memory of an fmpz. Probably not needed due to ocamls GC.*)
  val clear : t -> unit

  (** Get a string representation of the fmpz.*)
  val get_str : t -> string

  (** Set a fmpz to the given string with repsect to the given radix.*)
  val set_str : t -> string -> int -> unit

  (** Initialize and set an fmpz to the given string with respect to the given radix.*)
  val init_set_str : string -> int -> t

  (** Get an fmpz value from an mpz value*)
  val fmpz_from_mpz : C.mpz_t -> t

  (** Get an mpz value from an fmpz value*)
  val mpz_from_fmpz : t -> C.mpz_t

  (** [set a b] sets [a] to the value [b]*)
  val set : t -> t -> unit

end

(** Arb complex ball's. *)
module Acb : sig
  type t

  (** Get a fresh acb ready for use.*)
  val init : unit -> t

  (** Clear and free the memory of an acb. Probably not needed due to ocamls GC.*)
  val clear : t -> unit

  (** Set an acb to the given real and imaginary integers.*)
  val set_fmpz_fmpz : t -> Fmpz.t -> Fmpz.t -> unit

  (** Initialize and set an acb to the give real and imaginary integers. *)
  val init_set_fmpz_fmpz : Fmpz.t -> Fmpz.t -> t

  (** Calculates the complex logarithm of the given input to the given precision.*)
  val log : t -> int -> t

  (** Each component of an acb has a midpoint and a width value. This returns the midpoint
      value of the real component as integers (m, e). The midpoint is m * 2^e.*)
  val get_real_mid_fmpz : t -> Fmpz.t * Fmpz.t

  (** Each component of an acb has a midpoint and a width value. This returns the midpoint
    value of the imag component as integers (m, e). The midpoint is m * 2^e.*)
  val get_imag_mid_fmpz : t -> Fmpz.t * Fmpz.t

  (** Gives the value pi to a given precision.*)
  val pi : int -> t

  (** Adds two acbs to the given precision.*)
  val add : t -> t -> int -> t

  (** Subtracts two acbs to the given precision.*)
  val sub : t -> t -> int -> t

  (** Multiplies two acbs to the given precision.*)
  val mul : t -> t -> int -> t

  (** Divides two acbs to the given precision.*)
  val div : t -> t -> int -> t

  (** [mul_si a i prec]. Multiplies [a] by the integer [i] to the given precision.*)
  val mul_si : t -> int -> int -> t

  (** [div_si a i prec]. Divides [a] by the integer [i] to the given precision.*)
  val div_si : t -> int -> int -> t
end

(** Flint's integer polynomials.*)
module Fmpz_poly : sig
  type t

  (** Initializes an fmpz poly ready for use.*)
  val init : unit -> t

  val clear : t -> unit

  (** From the flint documentation. [set_coef poly n x], sets coefficient [n] of [poly] to the fmpz value [x]. 
      Coefficient numbering starts from zero and if [n] is beyond the current length of [poly] then the polynomial 
      is extended and zero coefficients inserted if necessary.*)
  val set_coef : t -> int -> Fmpz.t -> unit

  (** [get_coef poly n], gets the [n]-th coefficient of [poly]. Coefficient numbering is from zero and if 
  [n] is set to a value beyond the end of the polynomial, zero is returned.*)
  val get_coef : t -> int -> Fmpz.t


  (** From the arb documentation. 
      Writes to [roots] all the real and complex roots of the polynomial [poly], computed to at least [prec] accurate bits. 
      The root enclosures are guaranteed to be disjoint, so that all roots are isolated.
      The real roots are written first in ascending order (with the imaginary parts set exactly to zero). 
      The following nonreal roots are written in arbitrary order, but with conjugate pairs grouped together 
      (the root in the upper plane leading the root in the lower plane).
      The input polynomial must be squarefree. For a general polynomial, compute the squarefree part 
      or do a full squarefree factorization to obtain the multiplicities of the roots.
      All roots are refined to a relative accuracy of at least prec bits. The output values will generally have higher actual 
      precision, depending on the precision needed for isolation and the precision used internally by the algorithm.
      This implementation should be adequate for general use, but it is not currently competitive with state-of-the-art isolation 
      methods for finding real roots alone.*)
  val get_complex_roots : t -> int -> Acb.t list
end 

(** Flint's integer matrices.*)
module Fmpz_mat : sig
  type t

  (** Initialize an integer matrix of the given size.*)
  val init : int -> int -> t

  val clear : t -> unit

  (** Sets the i, j th entry to the given fmpz value.*)
  val set_entry : t -> Fmpz.t -> int -> int -> unit

  (** Returns a copy of the i, j th entry of the given matrix.*)
  val get_entry : t -> int -> int -> Fmpz.t

  (** The number of rows of the matrix. *)
  val nb_rows : t -> int

  (** The number of cols of the matrix. *)
  val nb_cols : t -> int

  (** [lll_original mat (delta_n, delta_d) (eta_n, eta_d)] performs in place
      lll reduction to the given matrix. Delta and eta are parameters to the
      lll algorithm. The algorithm gives gaurentees for delta in the range (0.25, 1\].
      Most often delta = 0.75. Eta needs to be in the range (0.5, sqrt(delta)).
      Most often eta = 0.51*)
  val lll_original : t -> int * int -> int * int -> unit

  (** [lll_storjohann mat (delta_n, delta_d) (eta_n, eta_d)] performs in place
      lll reduction to the given matrix. Delta and eta are parameters to the
      lll algorithm. The algorithm gives gaurentees for delta in the range (0.25, 1\].
      Most often delta = 0.75. Eta needs to be in the range (0.5, sqrt(delta)).
      Most often eta = 0.51. According to the flint documentation this algorithm has
      better complexity in the lattice dimension compared to the original algorithm.*)
  val lll_storjohann : t -> int * int -> int * int -> unit
end