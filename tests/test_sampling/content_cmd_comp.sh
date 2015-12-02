diff -B <(grep -vE '^\s*(#|parfile|protein|solute|solvent|$)' run_bnd.cmd) <(grep -vE '^\s*(#|parfile|protein|solute|solvent|$)' $PROTOMSHOME/tests/sampling/run_bnd.cmd)
