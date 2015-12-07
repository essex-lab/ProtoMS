diff -B <(grep -vE '^\s*(#|parfile|protein|solute|solvent|$)' run_free.cmd) <(grep -vE '^\s*(#|parfile|protein|solute|solvent|$)' $PROTOMSHOME/tests/RETI_dbl/run_free.cmd)
