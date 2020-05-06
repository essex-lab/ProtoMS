
c################################
c
c  Variables that are used for
c  this 'split' type functions
c
c  (C) Christopher Woods
c
c################################

c      the maximum number of words from a split line
       integer SPLITTER_MAXSPLIT
c      the maximum number of characters in each word
       integer SPLITTER_MAXLEN
       
       parameter(SPLITTER_MAXSPLIT=50)
       parameter(SPLITTER_MAXLEN=128)
       
c      the number of words returned by the split
       integer nwords
c      the string array that holds all of the words - the
c      length of each string must match SPLITTER_MAXLEN
       character*128 words(SPLITTER_MAXSPLIT)
       
c      the number of words returned by split2
       integer nwords2
c      the string array that holds all of the words
       character*128 words2(SPLITTER_MAXSPLIT)
       common /SPLITTERVARS/ words,nwords,words2,nwords2
       save /SPLITTERVARS/
