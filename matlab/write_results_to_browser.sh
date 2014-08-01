#! /bin/sh

DBFILE="$HOME/public_html/synergy/browser/database.sqlite"

rm -f "$DBFILE"

#    python insert_supplementary_data.py -t 0 -f "$f" -c 2 --tf="$reg" database.sqlite "ischip"

SUPPFILE="results/pol2max_and_delays_2013-08-30.txt"
SUPPFILE2="../R/gene_structures.txt"

FIGNAME_POL2="http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots/\${probe_name}_hmc_2013-08-30.png"
FIGNAME_PREMRNA="http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots_2013-11-05/\${probe_name}_hmc_2013-11-05.png"

RESULTFILE="results/hmc_results_to_browser_2013-08-30.txt"
FIGNAME=$FIGNAME_POL
REGNAME="pol2"

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
http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots/\${probe_name}_hmc_2013-08-30.png, Pol II model, , 3
http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots_2013-11-05/\${probe_name}_hmc_2013-11-05.png, pre-mRNA model, , 2
http://www.cs.helsinki.fi/u/ahonkela/synergy/premrna_plots/\${probe_name}_premrna.png, pre-mRNA data, , 1
EOF


RESULTFILE="results/hmc_results_to_browser_2013-09-04.txt"
FIGNAME="http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots_2013-09-04/\${probe_name}_hmc_2013-09-04.png"
REGNAME="pol2-norm"

insert_results.py -f $RESULTFILE --log-likelihood-column=2 "$DBFILE" "pol2rnaseq" "2.5% quantile" "$REGNAME" "$FIGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=3 "$DBFILE" "pol2rnaseq" "25% quantile" "$REGNAME" "$FIGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=4 "$DBFILE" "pol2rnaseq" "50% quantile" "$REGNAME" "$FIGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=5 "$DBFILE" "pol2rnaseq" "75% quantile" "$REGNAME" "$FIGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=6 "$DBFILE" "pol2rnaseq" "97.5% quantile" "$REGNAME" "$FIGNAME"
insert_results.py -f $RESULTFILE --log-likelihood-column=7 "$DBFILE" "pol2rnaseq" "begdev" "$REGNAME" "$FIGNAME"
insert_results.py -f "$SUPPFILE" --log-likelihood-column=5 "$DBFILE" "pol2rnaseq" "corr" "$REGNAME" "$FIGNAME"



RESULTFILE="results/hmc_results_to_browser_2013-11-05.txt"
FIGNAME=$FIGNAME_PREMRNA
REGNAME="pre-mRNA"

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
http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots/\${probe_name}_hmc_2013-08-30.png, Pol II model, , 2
http://www.cs.helsinki.fi/u/ahonkela/synergy/hmc_plots_2013-11-05/\${probe_name}_hmc_2013-11-05.png, pre-mRNA model, , 3
http://www.cs.helsinki.fi/u/ahonkela/synergy/premrna_plots/\${probe_name}_premrna.png, pre-mRNA data, , 1
EOF




#    python insert_results.py -f "$f" --log-likelihood-column=6 --experiment-set="2010-11-19" --figure-name="Models" --figure-filename='http://users.ics.tkk.fi/ahonkela/dros/expro_pics/${probe_name}_exp.png' database.sqlite "Drosophila Affy series" "posterior2/4" "$reg" 'http://users.ics.tkk.fi/ahonkela/multitf/model_pics_2010-11-19/multitf8aTestMCMC24-Nov-2010GeneExp${probe_name}Rep3.png'
#    python insert_results.py -f "$f" --log-likelihood-column=7 --experiment-set="2010-11-19" --figure-name="Models" --figure-filename='http://users.ics.tkk.fi/ahonkela/dros/expro_pics/${probe_name}_exp.png' database.sqlite "Drosophila Affy series" "posterior32" "$reg" 'http://users.ics.tkk.fi/ahonkela/multitf/model_pics_2010-11-19/multitf8aTestMCMC24-Nov-2010GeneExp${probe_name}Rep3.png'
#    python insert_figures.py -d ',' -f figures.txt --experiment-set="2010-11-19" database.sqlite "marll" "$reg" "Drosophila Affy series"
#    python insert_figures.py -d ',' -f figures.txt --experiment-set="2010-11-19" database.sqlite "posterior2/4" "$reg" "Drosophila Affy series"
#    python insert_figures.py -d ',' -f figures.txt --experiment-set="2010-11-19" database.sqlite "posterior32" "$reg" "Drosophila Affy series"
#done
#echo "done"

insert_aliases.py -f $HOME/projects/pol2rnaseq/data/ensembl69_aliases.txt "$DBFILE" "HGNC" "pol2rnaseq"

#python insert_supplementary_data.py -t 0 -f dros_twi_results.txt -c 3 database.sqlite "hasinsitu"
#python insert_supplementary_data.py -t 0 -f dros_twi_results.txt -c 4 database.sqlite "isinsitu"

chmod 444 "$DBFILE"

