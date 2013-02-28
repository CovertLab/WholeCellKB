#Updates whole-cell knowledge base haystack indices
#
#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 2012-07-17

export DJANGO_SETTINGS_MODULE=settings
export PYTHONPATH=.:..
python2.7 manage.py update_index --age=2 > log/haystack.index.log
chown -R WholeCell:WholeCell xapian_index
chmod ugo+w xapian_index
