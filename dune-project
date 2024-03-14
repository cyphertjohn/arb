(lang dune 2.9)

(name arbduet)

(generate_opam_files false)

(source
 (github cyphertjohn/arb))

(authors "John Cyphert")

(maintainers "John Cyphert")

(license MIT)


(package
 (name arbduet)
 (synopsis "Provides a limited interface to the arb library. Required for the duet program analyzer")
 (depends
   (ctypes (and :build (>= 0.22.0)))
   (ctypes-foreign (and :build (>= 0.22.0)))
   (ppx_cstubs :build)
 )
 (depopts
   (zarith :build)
 ))

