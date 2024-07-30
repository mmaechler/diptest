- The vignette `inst/doc/diptest-issues.Rnw`
  has mentioned the "new" (May 2011)  dip.test() function;
  Still have not finalized *analyzing* the simulations in ./stuff/

- Consider an analogue qnormDiptab  which is constructed using
   `rnorm(.)` instead of `runif(.)` simulations.
    - This idea is old; Werner Stuetzle's student, Jeremy Tantrum, 
	  did things in this direction in 2003; see, `stuff/jeremy-unimodality.R`

- Visualize the l.c.m. and g.c.m. and the modal interval !
	 g.c.m = greatest convex minorant =: G(x)
     l.c.m = least   concave majorant =: H(x)
