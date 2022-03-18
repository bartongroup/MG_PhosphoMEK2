#ssh www-prod mkdir -p /var/www/html/dag.compbio.dundee.ac.uk/share/marek/phosphomek2
#scp .htaccess doc/analysis.html www-prod:/var/www/html/dag.compbio.dundee.ac.uk/share/marek/phosphomek2

rsync -rvm --include='*/' --include='.htaccess' --include='doc/analysis.html' --include='tab/*' --exclude='*' . cluster:/cluster/gjb_lab/mgierlinski/public_html/phosphomek2