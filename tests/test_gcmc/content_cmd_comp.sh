diff -B <(grep -vE '^\s*(#|parfile|protein1|solute1|solvent1|$)' run_bnd.cmd) <(grep -vE '^\s*(#|parfile|protein1|solute1|solvent1|$)' $PROTOMSHOME/tests/gcmc/run_bnd.cmd)
