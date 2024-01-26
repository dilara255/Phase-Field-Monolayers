Copy the files to ../controlScripts first.

Multiple copies of runner can be created and set to different parameters, and runBatch can call a sequence of such files.

Multiple instances of runBatch can also be ran in parallel. The simulation itself is single-threaded, but it may occupy some time from a second core for reasons. In some cases, this usage may be high (TODO: look into this).

makeGif is called from runBatch, but it can also be called directly using "python makeGif.py *folder*", where the folder is a folder containing eithr PGMs (created by the CL simulation runner) or jpgs. It may fail in  case there are several thousand image files in the folder called.

