diff -B <(grep -vE '^\s*(#|parfile|protein1|solute1|$)' run_bnd.cmd) <(grep -vE '^\s*(#|parfile|protein1|solute1|$)' $PROTOMSHOME/tests/sampling/run_bnd.cmd)
