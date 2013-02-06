Rachel
======

Authors:
  A.B. Hayes and D. Cline, University of Rochester

GUI for Gosia

Rachel is an interface to the Coulomb excitation fit code Gosia.  Gosia itself
is maintained by a collaboration between Warsaw University, Koeln and
Rochester.  More information is on the Gosia wiki:

    Gosia
        http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php/What_is_Gosia
    Rachel itself
        http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php/Rachel,_a_GUI_for_Gosia

Borrowed Code:
    Modified the filename completion class from user "samplebias."
        http://stackoverflow.com/questions/5637124/tab-completion-in-pythons-raw-input
        The original author is http://stackoverflow.com/users/538718/samplebias

    Code for differential Rutherford cross
    sections was translated from ruthx.for
    by C.Y. Wu.

    The Clebsch-Gordan coefficients function was
    translated from the fortran function "NED" by
    A. Quirantes, Dept.  of Applied Physics,
    University of Granada.
    NED version 20 May 2.003
    See http://www.ugr.es/~aquiran/codigos.htm

    The code elast.c (Oak Ridge) is used
    to calculate stopping powers using
    Zeigler's formulas.  The version
    distributed with Rachel has been
    modified to operate on the command
    line by J.M. Allmond.

    The moving Gaussian smoothing routine is by
    S. Harden, U. Florida.

    Rachel reads the Radware .ags ascii level
    scheme files of GLS Version 3.0.
    D. C. Radford, Sept 1999.
    Later versions should be compatible.


