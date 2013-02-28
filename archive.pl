#!/usr/bin/perl

#archive.pl
#  Generates zip archive of WholeCellKB code
#
#Author: Jonathan Karr, jkarr@stanford.edu
#Affiliation: Covert Lab, Department of Bioengineering, Stanford University
#Last updated: 2012-08-10

$url=trim(`svn info | grep URL | cut -f2 -d" "`);

#clean previous archive
`rm -rf WholeCellKB-*.zip`;

#archive code
`svn export $url WholeCellKB`;
`rm -rf WholeCellKB/public/fixtures/`;
`rm -rf WholeCellKB/.htaccess`;
`rm -rf WholeCellKB/archive.pl`;
`rm -rf WholeCellKB/DevelopersGuide.tex`;

`export PYTHONPATH=$PYTHONPATH:.`;
`export DJANGO_SETTINGS_MODULE=settings`;
`rm -rf static/public/doc`;
`svn up static/public/doc`;
`epydoc -o static/public/doc --name=WholeCellKB --url=http://wholecellkb.stanford.edu --html --graph all --parse-only --docformat plaintext public`;
`cp -R static/public/doc WholeCellKB/static/public/`;

`python2.7 manage.py graph_models public | dot -Tsvg -o static/public/img/data_model.svg`;
`python2.7 manage.py graph_models public | dot -Tpng -o static/public/img/data_model.png`;
`cp static/public/img/data_model.* WholeCellKB/static/public/img/`;

`zip WholeCellKB-code.zip -r WholeCellKB`;
`rm -rf WholeCellKB`;

#archive data
`mkdir WholeCellKB`;
`svn export $url/public/fixtures/data.fna WholeCellKB/data.fna`;
`svn export $url/public/fixtures/data.sql WholeCellKB/data.sql`;
`svn export $url/public/fixtures/data.xlsx WholeCellKB/data.xlsx`;
`svn export $url/README.txt WholeCellKB/README.txt`;
`zip WholeCellKB-data.zip -r WholeCellKB`;
`rm -rf WholeCellKB`;

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
