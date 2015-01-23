#! /bin/sh

DBFILE="$HOME/public_html/synergy/browser/database.sqlite"

rm -f "$DBFILE"

#    python insert_supplementary_data.py -t 0 -f "$f" -c 2 --tf="$reg" database.sqlite "ischip"

SUPPFILE="results/pol2max_and_delays_final.txt"
SUPPFILE2="../R/gene_structures.txt"

FIGNAME_POL2="http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots/\${probe_name}_hmc_final.png"

ORIGRESULTFILE="results/hmc_results_to_browser_final.txt"
GENEFILE="../R/analysed_genes.txt"
RESULTFILE="results/hmc_results_to_browser_final.filtered.txt"
FIGNAME=$FIGNAME_POL
REGNAME="pol2"

grep -f $GENEFILE $ORIGRESULTFILE > $RESULTFILE

insert_results.py -f $RESULTFILE --log-likelihood-column=2 "$DBFILE" "pol2rnaseq" "2.5% quantile" "$REGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=3 "$DBFILE" "pol2rnaseq" "25% quantile" "$REGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=4 "$DBFILE" "pol2rnaseq" "50% quantile" "$REGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=5 "$DBFILE" "pol2rnaseq" "75% quantile" "$REGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=6 "$DBFILE" "pol2rnaseq" "97.5% quantile" "$REGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=7 "$DBFILE" "pol2rnaseq" "begdev" "$REGNAME"
insert_results.py -f "$SUPPFILE" --log-likelihood-column=5 "$DBFILE" "pol2rnaseq" "corr" "$REGNAME"
insert_results.py -f "$SUPPFILE" --log-likelihood-column=3 "$DBFILE" "pol2rnaseq" "pol2tmax" "$REGNAME"
insert_results.py -f "$SUPPFILE2" --log-likelihood-column=2 "$DBFILE" "pol2rnaseq" "maxTrLength" "$REGNAME"
insert_results.py -f "$SUPPFILE2" --log-likelihood-column=3 "$DBFILE" "pol2rnaseq" "lastProportion" "$REGNAME"
insert_figures.py -d ',' "$DBFILE" "pol2tmax" "$REGNAME" "pol2rnaseq" <<EOF
http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots/\${probe_name}_hmc_final.png, Pol II model, , 2
http://www.cs.helsinki.fi/u/ahonkela/synergy/premrna_plots/\${probe_name}_premrna.png, pre-mRNA data, , 1
EOF

insert_aliases.py -f $HOME/projects/pol2rnaseq/data/ensembl69_aliases.txt "$DBFILE" "HGNC" "pol2rnaseq"

chmod 444 "$DBFILE"

