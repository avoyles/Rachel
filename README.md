Rachel, GUI for Gosia
======

[Mind-blowing demo video](http://www.pas.rochester.edu/~hayes/beta_rachel/main_ad.html)

![Rachel snapshot](http://www-user.pas.rochester.edu/~gosia/mediawiki/images/4/41/Guisnapshot.png)


Authors
------

  Code and Physics:  A.B. Hayes, University of Rochester

  Physics:           D. Cline, University of Rochester

About
------

Rachel is an interface to the Coulomb excitation fit code Gosia.  Gosia itself
is maintained by a collaboration between Warsaw University, Koeln and
Rochester.  More information is on the [Gosia wiki](http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php/Rachel,_a_GUI_for_Gosia)

[Gosia Wiki: What is Gosia?](http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php/What_is_Gosia)

Borrowed Code
------

Filename completion by user "samplebias."
http://stackoverflow.com/questions/5637124/tab-completion-in-pythons-raw-input
The original author is http://stackoverflow.com/users/538718/samplebias
This code was modified to understand "..", "~" and to correctly complete through subdirectories.

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

