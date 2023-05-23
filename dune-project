(lang dune 2.9)

(name arb)

(generate_opam_files false)

(source
 (github cyphertjohn/arb))

(authors "John Cyphert")

(maintainers "John Cyphert")

(license MIT)


(package
 (name arb)
 (synopsis "Provides a limited interface to the arb library")
 (depends
   (ctypes :build)
   (ctypes-foreign :build)
   (ppx_cstubs :build)
 )
 (depopts
   (zarith :build)
 ))

