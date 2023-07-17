(lang dune 2.9)

(name arbduet)

(generate_opam_files true)

(source
 (github cyphertjohn/arb))

(authors "John Cyphert")

(maintainers "John Cyphert")

(license MIT)


(package
 (name arbduet)
 (synopsis "Provides a limited interface to the arb library. Required for the duet program analyzer")
 (depends
   (ctypes :build)
   (ctypes-foreign :build)
   (ppx_cstubs :build)
 )
 (depopts
   (zarith :build)
 ))

