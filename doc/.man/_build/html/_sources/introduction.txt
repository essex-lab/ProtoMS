*************
Introduction
*************

ProtoMS is short for “Prototype Molecular Simulation”, and is a software package that was originally designed by Christopher Woods to perform protein-ligand binding free energy calculations during his PhD. Julien Michel latter added numerous features and used the program extensively during his PhD. The program is now routinely used by several members of the research group of Jonathan Essex. This document has been written to try and explain how to use ProtoMS.

This document is presented as a hyperlinked PDF file, and this is the format that I recommend you use when you read it. This is because this document may change between versions of ProtoMS, and you do not want to have to keep printing it out, and viewing it electronically means that you can take advantage of the search facilities of your PDF viewer, and you can click on highlighted links to jump around between sections. You should also see bookmarks for all of the chapters, sections and subsections in the left pane, and may also see thumbnails of all of the pages.

The user manual has been written as a reference manual, with extensive hyperlinking to allow you to quickly dip in and out to find the information you need. While you could read it from start to end, it would be a boring and repetitive read and you probably wouldn’t learn much! I recommend that you play with the examples that come with ProtoMS. You can then use the links in those descriptions to dip in and out of the user manual, thus obtaining a more detailed knowledge of how the examples, and thus ProtoMS, work.


===========
Formatting
===========

The following formats are used throughout this document. Program commands or contents of files will be written in shaded boxes, e.g. ::

    temperature float

where float is a floating point option to the command. If this option is given a value (e.g. 25.0, then it is written like this ::

    temperature 25.0

The following options are standard to many commands

* **float** A floating point number.

* **integer** An integer. Most integer options given to ProtoMS are positive integers, greater than 0. This will always be made clear with the command.

* **logical** A logical, true or false option. Possible values for this option are true or false, yes or no or on or off, depending on your personal preference.

* **filename** This is the name of a file. Note that while ProtoMS is mostly case insensitive, file handling is dependent on the operating system you are using, so the filenames may be case sensitive. UNIX/Linux are examples of operating systems where case is important, while case is not important for Windows.


