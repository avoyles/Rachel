Rachel, GUI for Gosia
======

Rachel is currently being cleaned up and made PEP8 compliant in preparation to
turn over maintenance to users.  Please fork if you are interested in doing
maintenance or improvements.

This is a DEVELOPMENT version.  You may find bugs.  If you want a
completely-tested version, download [the release version](http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php/Main_Page#Downloads)

:boom: [Mind-blowing demo video](http://youtu.be/moVVC-GODzQ)

![Rachel snapshot](http://www-user.pas.rochester.edu/~gosia/mediawiki/images/4/41/Guisnapshot.png)


Authors
------

  Code and Physics:  A.B. Hayes, University of Rochester

  Physics:           D. Cline, University of Rochester

About
------

Rachel is an interface to the Coulomb excitation fit code Gosia.  More
information is on the [Gosia wiki](http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php?title=Rachel_GUI&oldid=830).

Gosia itself is maintained by a collaboration between Warsaw University, Koeln
and Rochester.  [Gosia Wiki: What is Gosia?](http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php?title=What_is_Gosia%3F&oldid=297).


Fast Setup
------

(A complete installer is coming soon.)

1. Un-tar the archive, or, better yet, clone the [github repo](http://github.com/adamhayes/Rachel), so you can contribute.

2. Move into the directory containing this README file (the directory containing rachel.py).

3. Run the compile-all.sh script.  This is usually done by typing
   './compile-all.sh' on Linux systems.  compile-all will compile the Gosia
   version distributed with Rachel and the Elast version modified for
   command-line use by J.M. Allmond.

   The compile-all script will also create a '.rachel_setup' file in your home
   directory, which can be copied to another account and/or a working
   directory.

4. './compile-all.sh' clean will remove the compiled executables.

5. Make a working directory and cd into it.  This directory should not be
   shared by different calculations, i.e., you can calculate whatever you can
   put in one Rachel session at a time, but you shouldn't switch back and forth
   between other 'pickle.jar' files for unrelated calculations.  (This is the
   reason that the session file name is fixed.  The working directory will be
   used in a "stateful" way, analogous to the statefulness of any Gosia working
   directory.)

6. Start Rachel, by typing 'python [path]rachel.py' at the Linux prompt, e.g.
   'python /home/hayes/Rachel/rachel.py'  Do not run Rachel in the background.
   Rachel will tell you if any python libraries need to be installed

7. Begin by reading a level scheme file.  At this point, you should refer to
   the Rachel section of the [Gosia Manual](http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php/Main_Page#Downloads)
   and new versions of the tutorial videos.  (The beta-version tutorial videos
   have been taken down.  New videos will be posted soon.)

Notes
------

* The current version of Rachel will check the Rochester server for version
  updates or other important information at startup.

* Rachel keeps at least 5 backup versions of your work, in case of user error,
  however, you can recover from almost any exception using the Reactivate
  button and Undo/Redo.  The Undo/Redo buttons are very robust.

* If you somehow managed to crash the GUI (see above), you can also recover
  with 'python [path]/rachel.py -r', but you have to do this the first time you
  restart Rachel after the crash!  You might want to use the -r option if you
  did a lot of work without saving before the crash.

* You could make your working directory into a git repo, if you know the basics
  of git.  This is very useful for complicated work where you might want to
  revert to some old calculation, model, set of matrix elements, etc.  git
  restores the proper state of the working directory!

* Check out the new "Logs" button!

Borrowed Code and Credits
------

Thanks to Mitch Allmond for the compiling script.

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

