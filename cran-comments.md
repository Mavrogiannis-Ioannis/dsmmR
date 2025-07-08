# version 1.0.7 --------------------------------------------------------------------

# R Version: 4.5.1

# R CMD check results.

No ERRORs, WARNINGs or NOTEs.


# R CMD check results through `rhub::rhub_check()`.

There are no WARNINGs or NOTEs, only ERRORs.
   
   Similar to previous accepted version.

### Regarding a local Rstudio ERROR, after all CMD check passes.

Solution: `devtools::check(document = FALSE)` 

Info: `https://github.com/rstudio/rstudio/issues/15805`



# version 1.0.5 --------------------------------------------------------------------

# R Version: 4.4.1


# R CMD check results.

No ERRORs, WARNINGs or NOTEs.



# R CMD check results through `rhub::rhub_check(branch = 'master')`.

There are no WARNINGs or NOTEs, only ERRORs.
   
   Errors:
   All errors are associated with github actions.
   (1) * `rchk` (input number 19):
      Github action fails to finish, with the error message:
      `No files were found with the provided path: check. No artifacts will be uploaded.`
   (2-6) * `c23`, `clang16`, `clang17`, `clang18`, `clang19` (input numbers 6, 8, 9, 10, 11):
      Github action fails to finish, in part with the error message:
      `undefined symbol: __asan_option_detect_stack_use_after_return`

### Regarding the ERRORs:

   Regarding (1) : 
      As seen in the issue below, this could be a bug in `rchk`'s side:
      <https://github.com/actions/upload-artifact/issues/232#issuecomment-964235360>
      
   Regarding (2-6) : 
      Similar problems were encountered before from other users:
      <https://gcc.gnu.org/pipermail/gcc-help/2019-March/136892.html>



# version 1.0.2 --------------------------------------------------------------------

# R Version: 4.3.1

# R CMD check results.

No ERRORs, WARNINGs or NOTEs.


# R CMD check results through `rhub::check_for_cran()`.

There are no ERRORs or WARNINGs, only NOTEs.

   Notes:
   (1) * Solutions similar to what was described in v0.0.96.



# version 1.0.1 --------------------------------------------------------------------

# R Version: 4.2.2

# R CMD check results.

No ERRORs, WARNINGs or NOTEs.


# R CMD check results through `rhub::check_for_cran()`.

There are no ERRORs or WARNINGs, only NOTEs.

   Notes:
   (1) * Solutions similar to what was described in v0.0.96.


# version 0.0.96 --------------------------------------------------------------------

# R Version : 4.2.2

# R CMD check results.

No ERRORs, WARNINGs or NOTEs.


# R CMD check results through `rhub::check_for_cran()`.

There are no ERRORs or WARNINGs, only NOTEs.



## Windows Server 2022, R-devel, 64 bit

   Notes:
   (1) * Possibly misspelled words in DESCRIPTION:
           Barbu (27:5, 31:5)
           Limnios (27:18)
           Vergne (29:5, 31:19)
   
   (2) * checking HTML version of manual ... NOTE
            Skipping checking math rendering: package 'V8' unavailable
   (3) * checking for detritus in the temp directory ... NOTE
          Found the following files/directories:
            'lastMiKTeXException'

### Regarding the NOTEs
   Regarding (1) : These are real names.
   
   Regarding (2) : The solution regarding this NOTE seems to be a bug regarding 
   version 4.2.0 of R: <https://groups.google.com/g/r-sig-mac/c/7u_ivEj4zhM>.
   
   Regarding (3) : As seen in the following ticket, it seems like this is a bug 
   on MikTex's end:
   <https://github.com/r-hub/rhub/issues/503#issuecomment-1014483806>.


   By using the command :
   
   ```
   rhub::check_for_cran(env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "true", 
                                     `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true", 
                                     `_R_CHECK_RD_VALIDATE_RD2HTML_` = "false"))
   ```

   These NOTEs do not persist.



## One more _unique_ NOTE may appear in other releases except WINDOWS.

   For example, in:
   
## Ubuntu Linux 20.04.1 LTS, R-release, GCC
   
   
   (4) * checking HTML version of manual ... NOTE
            Skipping checking HTML validation: no command 'tidy' found


### Regarding the NOTE:

   Regarding (4) : This appears to be similar to previous NOTE (2).


