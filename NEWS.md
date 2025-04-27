# R2sample 4.0.0
                It is now possible to use a routine to generate new simulated data to
                find p values, the parametric bootstrap approach. This is needed for the
                goodness-of-fit/twosample hybrid problem.
                Also some changes to the hidden interior routines   
# R2sample 3.1.1
                Fixed a minor bug
# R2sample 3.1.0
                Improved the routines that do power calculations for better speed. 
                Included a timing routine to see whether using a single core is 
                faster than using multiple processing.
                Some minor changes to other routines.
# R2sample 3.0.0
                Added routine to allow benchmarking of new user
                supplied tests. Some minor changes to other routines.
                   
# R2sample 2.2.0
                some minor bug fixes, additions to vignette
                
# R2sample 2.1.0
                Added routines for calculating p values adjusted for simultaneous testing
                
# R2sample 1.1.0
                fixed a bug in calculation of chi square test, made cpp routines 
                invisible, fixed issue with help titles
                  
# R2sample 1.0.0
 
                fixed a code problem in TS_disc_cpp on line 146
                made some changes to the arguments of twosample_power   

# R2sample 0.0.4

                changed line 66 in TS_cont_cpp.cpp from 
                   while ( (x[j]<=sxy[i]) && (j<nx))
                to
                   while ( (j<nx) && (x[j]<=sxy[i]) )   
                to avoid heap-buffer-overflow error.   
                   

# R2sample 0.0.3

                Changed & to && in a TS_cont_cpp.cpp

# R2sample 0.0.2

                Changed | to || in a number of the C++ routines per request from the CRAN Team  

# R2sample 0.0.1

* Added a `NEWS.md` file to track changes to the package.

* Oct 10, 2022: Added \value to run_shiny(), added () to function names
*               eliminated \dontrun(), added toy examples
*               changed cat to message
*               changed Maintainer to Authors@R  
* Oct. 11, 2022: Eliminated all \dontrun and \donttest, added toy examples
                 Searched and eliminated all empty spaces in DESCRIPTION that I could find.  


              
