'''
Whole-cell knowledge base forms

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from django import forms
from public.helpers import getObjectTypes, getModels
from public.models import Species, SpeciesComponent

class SearchForm(forms.Form):
	q = forms.CharField(
		required = True,
		widget = forms.TextInput,
		label = 'query',
		help_text = 'Enter search term(s)',
		initial = ''
		)

class ExportDataForm(forms.Form):
	FORMAT_CHOICES = (
		('bib', 'BibTex'),
		('xlsx', 'Excel'),
		('html', 'HTML'),
		('json', 'JSON'),		
		('pdf', 'PDF'),
		('xml', 'XML'),
	)
	
	species = forms.ChoiceField(
		required = False,
		widget = forms.RadioSelect, 
		label = 'species',
		help_text = 'Select species to export'
		)
	model_type = forms.MultipleChoiceField(
		required = False,
		widget = forms.CheckboxSelectMultiple, 
		label = 'type',
		help_text = 'Select entry types to export'
		)
	all_model_types = forms.ChoiceField(
		choices = (
			('True', 'All'),
		),
		required = False,
		widget = forms.RadioSelect, 
		label = 'All',
		help_text = 'Select all entry types'
	)
	format = forms.ChoiceField(
		choices = FORMAT_CHOICES, 
		initial = FORMAT_CHOICES[0][0],
		required = True,
		widget = forms.RadioSelect, 		
		label = 'format',
		help_text = 'Select an output format'
		)
	
	def __init__(self, *args, **kwargs):
		super(ExportDataForm, self).__init__(*args, **kwargs)
		
		choices = []
		for species in Species.objects.values('wid', 'name').all():
			choices.append((species['wid'], species['name'], ))
		self.fields['species'].choices = choices
		
		model_types = getObjectTypes()
		models = getModels()
		choices = []
		for model_type in model_types:
			choices.append((model_type, models[model_type]._meta.verbose_name_plural, ))
		self.fields['model_type'].choices = choices
		self.fields['model_type'].initial = model_types

class ImportDataForm(forms.Form):
	species = forms.ChoiceField(
		required = False,
		widget = forms.Select, 
		label = 'species',
		help_text = 'Select species to export'
		)
	new_species = forms.CharField(
		required = False,
		widget = forms.TextInput,
		label = 'species',
		help_text = 'Choose a name for a new PGDB'
		)	
	file = forms.FileField(
		required = True,
		widget = forms.ClearableFileInput, 
		label = '(Excel or FASTA) file',
		help_text = 'Select a file (Excel or FASTA) to import'
		)
		
	def __init__(self, *args, **kwargs):
		super(ImportDataForm, self).__init__(*args, **kwargs)
		
		choices = []
		for species in Species.objects.values('wid', 'name').all():
			choices.append((species['wid'], species['name'], ))
		self.fields['species'].choices = choices
