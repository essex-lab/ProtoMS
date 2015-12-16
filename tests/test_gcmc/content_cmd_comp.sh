diff -B <(grep -vE '^\s*(#|parfile|protein|solute|solvent|$)' run_bnd.cmd) <(grep -vE '^\s*(#|parfile|protein|solute|solvent|$)' $PROTOMSHOME/tests/gcmc/run_bnd.cmd)
