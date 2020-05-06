*************
Introduction
*************

ProtoMS is short for “Prototype Molecular Simulation”, and is a software package that was originally designed by Dr. Christopher Woods to perform protein-ligand binding free energy calculations during his PhD. Dr. Julien Michel and Dr. Michael Bodnarchuk latter added numerous features and used the program extensively during their PhDs. Dr. Samuel Genheden, Dr. Richard Bradshaw, Dr. Gregory Ross, Dr. Chris Cave-Ayland, Dr. Ana Cabedo Martinez, Hannah Bruce-Macdonald and James Graham have since then made numerous additions to the code among them a complete revision of the tools used to setup and analyse the simulation results.

The program is routinely used by several members of the research group of Professor Jonathan Essex for the development of new techniques to perform free energy calculations. This document has been written to try and explain how to use ProtoMS.

The user manual has been written as a reference manual, with extensive hyperlinking to allow you to quickly dip in and out to find the information you need. While you could read it from start to end, it would be a boring and repetitive read and you probably wouldn’t learn much! We recommend that you engage with the tutorials that come with ProtoMS. You can then use the links in those descriptions to dip in and out of the user manual, thus obtaining a more detailed knowledge of how the examples, and thus ProtoMS, work.


===========
Formatting
===========

The following formats are used throughout this document. Program commands or contents of files will be written in monotype, e.g. ::

    temperature float

where float is a floating point option to the command. If this option is given a value (e.g. 25.0), then it is written like this ::

    temperature 25.0

The following options are standard to many commands

* **float** A floating point number.

* **integer** An integer. Most integer options given to ProtoMS are positive integers, greater than 0. This will always be made clear with the command.

* **logical** A logical, true or false option. Possible values for this option are true or false, yes or no or on or off, depending on your personal preference.

* **filename** This is the name of a file. Note that while ProtoMS is mostly case insensitive, file handling is dependent on the operating system you are using, so the filenames may be case sensitive. UNIX/Linux are examples of operating systems where case is important, while case is not important for Windows.


