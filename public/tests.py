'''
Whole-cell knowledge base tests

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from django.db import models as dj_models
from django.test import TestCase
from public import helpers, models

'''
TODO
- import/excel works correctly
'''

class Test(TestCase):
	def test_data_model(self):		
		for name, model in helpers.getModels().iteritems():
			validate_fields(model)
			validate_model(model)
		for tmp in helpers.getEntryDatas():
			validate_model(tmp[1])
				
def validate_fields(model):
	for field in model._meta.fields + model._meta.many_to_many:
		if not field.editable or field.auto_created:
			continue
		if issubclass(model, models.Entry) and isinstance(field, dj_models.OneToOneField) and field.rel.parent_link:
			continue
		if field.unique and not (issubclass(model, models.Entry) and field.name == 'wid'):
			raise TypeError('%s model field %s cannot be unique' % (model.__name__, field.name, ))
		if isinstance(field, (
			dj_models.BooleanField, dj_models.CharField,
			dj_models.FloatField, dj_models.NullBooleanField, dj_models.IntegerField, 
			dj_models.PositiveIntegerField, dj_models.SlugField, dj_models.TextField,
			)):
			continue
		if isinstance(field, dj_models.OneToOneField):
			raise TypeError('Invalid type for %s field %s' % (model.__name__, field.name, ))
		if isinstance(field, (dj_models.ForeignKey, dj_models.ManyToManyField)):
			if issubclass(field.rel.to, models.Entry) or field.rel.to == models.Evidence:
				pass
			elif issubclass(model, models.Entry):
				validate_fields(field.rel.to)
			else:
				raise TypeError('Invalid type for %s field %s' % (model.__name__, field.name, ))

def validate_model(model):	
	if len(model._meta.unique_together) > 0:
		raise TypeError('%s model unique_together not supported' % model.__name__)
		
	if model.full_clean.im_func.__module__ != 'django.db.models.base':
		raise TypeError('Cannot override validation in %s model' % model.__name__)
	if model.clean_fields.im_func.__module__ != 'django.db.models.base':
		raise TypeError('Cannot override validation in %s model' % model.__name__)
	if model.clean.im_func.__module__ != 'django.db.models.base':
		raise TypeError('Cannot override validation in %s model' % model.__name__)
	if model.validate_unique.im_func.__module__ != 'django.db.models.base':
		raise TypeError('Cannot override validation in %s model' % model.__name__)