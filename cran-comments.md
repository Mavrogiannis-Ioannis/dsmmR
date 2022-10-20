# Test environments

* R CMD check results

There were no ERRORs or WARNINGs. There are 2 NOTEs:

* ONLY on windows...
checking for detritus in the temp directory
     'lastMiKTeXException'

As seen in the following ticket, it seems like this is a bug on MikTex's end.
<https://github.com/r-hub/rhub/issues/503#issuecomment-1014483806>


Running the command below eliminates the issues of spelling and HTML validation descibed in 2 and 3.
``` r
rhub::check_for_cran(env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "true", `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "false", `_R_CHECK_RD_VALIDATE_RD2HTML_` = "false"))
```
However, the issue at 1 remains.


