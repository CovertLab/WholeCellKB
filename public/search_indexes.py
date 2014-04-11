'''
Whole-cell knowledge base haystack indices

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from django.contrib.auth.models import User
from django.db.models.fields import BooleanField, NullBooleanField
from django.db.models.fields.related import ForeignKey, ManyToManyField
from haystack import site
from haystack.indexes import CharField, IntegerField, SearchIndex, ModelSearchIndex
from public.helpers import getModels, getModelDataFields
from public.models import Entry, SpeciesComponent
import datetime
import settings

def truncate_fields(obj):
	for key, val in obj.iteritems():
		if isinstance(val, dict):
			obj[key] = truncate_fields(val)
		elif isinstance(val, (list, tuple, )):
			idx = -1
			for subval in val:
				idx += 1
				obj[key][idx] = truncate_fields(subval)
		elif isinstance(val, (str, unicode)) and key not in ['content', 'text']:
			obj[key] = val[:225]
	return obj

class EntryIndex(ModelSearchIndex):
	def get_updated_field(self):
		return 'last_updated_date'
		
	# Hack to avoid error: "xapian.InvalidArgumentError: Term too long (> 245)"
	# See: https://groups.google.com/forum/?fromgroups#!topic/django-haystack/hRJKcPNPXqw
	def prepare(self, object):
		self.prepared_data = truncate_fields(super(EntryIndex, self).prepare(object))		
		return self.prepared_data
		
	class Meta:
		pass
		
class SpeciesComponentIndex(EntryIndex):
	species_id = IntegerField(model_attr='species__id')
	species_wid = CharField(model_attr='species__wid')
			
	class Meta:
		pass
		
def format_field_for_indexing(field=None, related=None, depth=0):
	name_obj = 'sub_'*depth + 'object'
	name_child = 'sub_' + name_obj
	if (field is None and related is None) or (field is not None and related is not None):
		return ''
		
	if related is not None:
		field_name = related.get_accessor_name()
		field_model = related.model
		field = related.field
	else:
		field_name = field.name
		if field.rel is not None:
			field_model = field.rel.to
			
	if field_name in ['id', 'created_user', 'created_date', 'last_updated_user', 'last_updated_date']:
		return ''
		
	if isinstance(field, ManyToManyField) or related is not None:
		if issubclass(field_model, Entry):			
			return ''
		else:
			result  = '\t'*depth + '{%% for %s in %s.%s.all %%}\n' % (name_child, name_obj, field_name, )
			for subfield in field_model._meta.fields + field_model._meta.many_to_many:
				result += format_field_for_indexing(field=subfield, depth=depth+1)
			result += '\t'*depth + '{% endfor %}\n'
			return result
	elif isinstance(field, ForeignKey):
		if issubclass(field_model, Entry):
			return ''
		else:
			result  = '\t'*depth + '{%% with %s=%s.%s %%}\n' % (name_child, name_obj, field_name, )
			for subfield in field_model._meta.fields + field_model._meta.many_to_many:
				result += format_field_for_indexing(field=subfield, depth=depth+1)
			result += '\t'*depth + '{% endwith %}\n'
			return result
	elif (field.choices is not None) and (len(field.choices) > 0) and (not isinstance(field, (BooleanField, NullBooleanField))):
		return '\t'*depth + '{%% get_choice_verbose_name "%s" "%s" "%s" %s.%s %%}\n' % (field.model._meta.app_label, field.model.__name__, field_name, name_obj, field_name, )
	else:
		if field.name != 'sequence':
			return '\t'*depth + '{{ %s.%s }}\n' % (name_obj, field_name, )
		
	return ''
		
for modelname, model in getModels().iteritems():
	#make full text template
	filename = '%s/templates/search/indexes/%s/%s_%s.txt' % (
		settings.ROOT_DIR, model._meta.app_label, modelname.lower(), 'text')
	f = open(filename, 'w')
	f.write('{% load templatetags %}\n')
	for field in getModelDataFields(model):
		f.write(str(format_field_for_indexing(field=field)))
	for related in model._meta.get_all_related_objects() + model._meta.get_all_related_many_to_many_objects():
		if related.get_accessor_name() != '+':
			f.write(str(format_field_for_indexing(related=related)))
	f.close()
	
	#create index class
	if issubclass(model, SpeciesComponent):
		index = type(model.__class__.__name__ + 'Index', (SpeciesComponentIndex, ), {})
	else:
		index = type(model.__class__.__name__ + 'Index', (EntryIndex, ), {})
	
	#register index	
	site.register(model, index)
