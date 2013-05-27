'''
Whole-cell knowledge base URL patterns

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from django.conf.urls.defaults import patterns, include, url
from public.views import *

# enable the admin:
from django.contrib import admin
admin.autodiscover()

# admin interface
urlpatterns = patterns('',
	url(r'^admin/doc/*', include('django.contrib.admindocs.urls')),
	url(r'^admin/*', include(admin.site.urls)),
)

# authentication
urlpatterns += patterns('public.views',
	url(r'^login/*$', 'login'),
	url(r'^login/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'login'),
	
	url(r'^logout/*$', 'logout'),
	url(r'^logout/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'logout'),
)

# public interface
urlpatterns += patterns('public.views',	
	url(r'^about/*$', 'about'),
	url(r'^about/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'about'),
	
	url(r'^tutorial/*$', 'tutorial'),
	url(r'^tutorial/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'tutorial'),
	
	url(r'^contributors/*$', 'contributors'),
	url(r'^contributors/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'contributors'),
	
	url(r'^contributor/(?P<username>[\w\d]+)/*$', 'contributor'),
	url(r'^contributor/(?P<username>[\w\d]+)/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'contributor'),
	
	url(r'^search/*', 'search'),
	url(r'^search/(?P<species_wid>[a-zA-Z0-9_\-]+)/*', 'search'),
	url(r'^list/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<model_type>\w+)/*$', 'list'),	
	url(r'^detail/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<wid>[a-zA-Z0-9_\-]+)/*$', 'detail'),
	url(r'^viewPropertyInSimulation/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<class_name>[a-zA-Z0-9_]+)/(?P<property_name>[a-zA-Z0-9_]+)/*$', 'viewPropertyInSimulation'),
	url(r'^edit/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<wid>[a-zA-Z0-9_\-]+)/*$', 'edit'),
	url(r'^delete/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<wid>[a-zA-Z0-9_\-]+)/*$', 'delete'),
	url(r'^add/(?P<species_wid>[a-zA-Z0-9_\-]+)/(?P<model_type>[a-zA-Z0-9_]+)/*$', 'add'),
	url(r'^add/(?P<model_type>[a-zA-Z0-9_]+)/*$', 'add'),
	
	url(r'^export/*$', 'exportData'),
	url(r'^export/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'exportData'),
	
	url(r'^import/*$', 'importData'),
	url(r'^import/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'importData'),
	
	url(r'^validate/(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'validate'),
	
	url(r'^sitemap.xml$', 'sitemap'),
	url(r'^sitemap_toplevel.xml$', 'sitemap_toplevel'),
	url(r'^sitemap_(?P<species_wid>[a-zA-Z0-9_\-]+).xml$', 'sitemap_species'),
	
	url(r'^$', 'index'),
	url(r'^(?P<species_wid>[a-zA-Z0-9_\-]+)/*$', 'index'),
)
