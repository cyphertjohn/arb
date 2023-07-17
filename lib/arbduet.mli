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

  (** [set a b] sets [a] to the value [b]*)
  val set : t -> t -> unit

  (** Get an fmpz value from an mpz value*)
  val fmpz_from_mpz : C.mpz_t -> t

  (** Get an mpz value from an fmpz value*)
  val mpz_from_fmpz : t -> C.mpz_t

  val mul_2exp : t -> int -> t

  (** Returns a copy. *)
  val init_set : t -> t

  (** Return an fmpz of an int *)
  val of_int : int -> t

  (** Returns an integer of an fmpz. If the input can't fit into an integer the function has undefined behavior.*)
  val to_int : t -> int

  val zero : unit -> t

  val one : unit -> t
  
  val cmp : t -> t -> int

  val cmp_si : t -> int -> int

  val equal : t -> t -> bool

  val equal_si : t -> int -> bool

  val neg : t -> t

  val add : t -> t -> t

  val add_si : t -> int -> t

  val sub : t -> t -> t

  val sub_si : t -> int -> t

  val mul : t -> t -> t
  
  val mul_si : t -> int -> t

  val divexact : t -> t -> t

  val divexact_si : t -> int -> t

  val pow_ui : t -> int -> t

  val gcd : t -> t -> t

  val lcm : t -> t -> t

end


(** Arb complex ball's. *)
module Acb : sig
  type t

  (** Get a fresh acb ready for use.*)
  val init : unit -> t

  (** Clear and free the memory of an acb. Probably not needed due to ocamls GC.*)
  val clear : t -> unit

  val set_fmpz : t -> Fmpz.t -> unit

  val init_set_fmpz : Fmpz.t -> t

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

  (** If x = a + bi, then [(c, d) = get_real_imag_mag_upper x] such that [c] and [d] are integers
      and |a| <= [c] and |b| <= [d]. Note the following warning from the arb docs: These functions 
      are unsafe: the user must check in advance that x is of reasonable magnitude. If x is infinite 
      or has a bignum exponent, an abort will be raised. If the exponent otherwise is too large or 
      too small, the available memory could be exhausted resulting in undefined behavior.*)
  val get_real_imag_mag_upper : t -> Fmpz.t * Fmpz.t

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

  val bits : t -> int

  val trim : t -> t

  val rel_error_bits : t -> int

  val rel_accuracy_bits : t -> int
  
  val rel_one_accuracy_bits : t -> int

  val neg : t -> t

  val pow_si : t -> int -> int -> t

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

  exception Incompatible_Dimensions
  exception Index_Out_of_Bounds

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

  (** Create a truncated identity matrix of the given size.
      @raise Index_Out_of_Bounds if either of the inputs are negative.*)
  val ident : int -> int -> t

  (** Create a zero matrix of the given size
      @raise Index_Out_of_Bounds if either of the inputs are negative. *)
  val zero : int -> int -> t

  (** [window m r1 c1 r2 c2] returns a copy of the [r2-r1] by [c2 - c1] submatrix whose [(0, 0)]
      entry is the [(r1, c1)] entry of the input.
      @raise Index_Out_of_Bounds if the indices are incompatible. *)
  val window : t -> int -> int -> int -> int -> t

  val equal : t -> t -> bool

  val is_zero : t -> bool
  
  val transpose : t -> t

  (** Concats two matrices vertically.
      @raise Incompatible_Dimensions if the inputs are incompatible.*)
  val concat_vertical : t -> t -> t
  
  (** Concats two matrices horizontall.
      @raise Incompatible_Dimensions if the inputs are incompatible.*)
  val concat_horizontal : t -> t -> t

  (** @raise Incompatible_Dimensions if the inputs are not the same size.*)
  val add : t -> t -> t

  (** @raise Incompatible_Dimensions if the inputs are not the same size.*)
  val sub : t -> t -> t
  
  val neg : t -> t
  
  val scalar_mult_si : t -> int -> t

  val scalar_mult : t -> Fmpz.t -> t

  (** Set A = B / c, where B is an fmpz_mat_t and c is a scalar respectively 
      which is assumed to divide all elements of B exactly.*)
  val divexact_si : t -> int -> t

  val divexact : t -> Fmpz.t -> t

  (** Multiplies the entries of the matrix by 2^e*)
  val scalar_mult_2exp : t -> int -> t

  (** Divides the matrix by 2^e rounded down towards 0.*)
  val scalar_tdiv_q_2exp : t -> int -> t

  (** @raise Incompatible_Dimensions if the inputs are not compatible.*)
  val mul : t -> t -> t

  (** Computes the Kronecker product of the given matrices.*)
  val kronecker : t -> t -> t

  (** @raise Invalid_argument "Negative Exponent"*)
  val pow : t -> int -> t

  (** The gcd of the elements of the matrix*)
  val content : t -> Fmpz.t

  (** The trace of the given matrix.*)
  val trace : t -> Fmpz.t

  (** The hnf of the given matrix.*)
  val hnf : t -> t

  (** [hnf_transform a = (h, u)] where [h] is the hnf of [a]
      and [u] is such that [ua = h]*)
  val hnf_transform : t -> t * t

  val is_in_hnf : t -> bool

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