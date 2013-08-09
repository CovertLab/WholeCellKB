'''
Whole-cell knowledge base views

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from copy import deepcopy
from django.contrib.auth import login as auth_login, logout as auth_logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from django.core.urlresolvers import reverse
from django.db.models import Count, Sum, Avg
from django.db.models.fields import BooleanField, NullBooleanField, AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField
from django.db.models.fields.related import OneToOneField, RelatedObject, ManyToManyField, ForeignKey
from django.db.models.query import EmptyQuerySet
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.utils.text import capfirst
from django.views.decorators.debug import sensitive_post_parameters
from django.views.decorators.cache import never_cache
from django.views.decorators.csrf import csrf_protect
from haystack.query import SearchQuerySet
from itertools import chain
from public import models
from public.forms import ExportDataForm, ImportDataForm
from public.helpers import getEntry, format_field_detail_view, objectToQuerySet, render_queryset_to_response, getObjectTypes, getModel, get_invalid_objects, get_edit_form_fields, get_edit_form_data
from public.helpers import validate_object_fields, validate_model_objects, validate_model_unique, save_object_data, batch_import_from_excel, readFasta, writeExcel
from public.models import Entry, Parameter, Reference, Species, SpeciesComponent, ModelProperty, SimulationProperty, Type
from public.helpers import PropertyDefinitionFilter, PygmentsStyle, PygmentsFormatter
import pygments
from pygments.lexers import MatlabLexer
from pygments.token import Token
from StringIO import StringIO
from urlparse import urlparse
import numpy
import os
import settings
import tempfile

MODEL_CODE_BASE_DIR = '/home/projects/WholeCell/simulation'

def index(request, species_wid=None):
	if species_wid is not None and species_wid != '':
		species = Species.objects.get(wid=species_wid)
	else:
		species = Species.objects.all()
		if len(species) > 0:
			species = species[0]
		else:
			species = None
	content = []
	if species is not None:		
		content.append([
			[0, 'Compartments', models.Compartment.objects.filter(species__id = species.id).count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Compartment'})],
		])
		
		chrs = models.Chromosome.objects.filter(species__id = species.id)
		chrcontent = chrs.aggregate(length=Sum('length'));
		try:
			gc_content = sum([chr.get_gc_content() * chr.length for chr in chrs]) / chrcontent['length']		
		except TypeError, e:
			gc_content = 0
		content.append([
			[0, 'Chromosomes', chrs.count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Chromosome'})],
			[1, 'Length', chrcontent['length'], 'nt'],
			[1, 'GC-content', ('%0.1f' % (gc_content * 100)), '%'],
		])
				
		tus = models.TranscriptionUnit.objects.filter(species__id = species.id).annotate(num_genes = Count('genes'))
		mons = tus.filter(num_genes__lte = 1)
		nPolys = tus.filter(num_genes__gt = 1).count()
		content.append([
			[0, 'Transcription units', tus.count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionUnit'})],			
			[1, 'Monocistrons', tus.count() - nPolys],
			[1, 'Polycistrons', nPolys],
		])
		
		content.append([
			[0, 'Genes', models.Gene.objects.filter(species__id = species.id).count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'})],
			[1, 'mRNA', models.Gene.objects.filter(species__id = species.id, type__wid='mRNA').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=mRNA'],
			[1, 'rRNA', models.Gene.objects.filter(species__id = species.id, type__wid='rRNA').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=rRNA'],
			[1, 'sRNA', models.Gene.objects.filter(species__id = species.id, type__wid='sRNA').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=sRNA'],
			[1, 'tRNA', models.Gene.objects.filter(species__id = species.id, type__wid='tRNA').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Gene'}) + '?type=tRNA'],
		])
		
		content.append([
			[0, 'Chromosome features', 
				models.ChromosomeFeature.objects.filter(species__id = species.id).count(),
				None,
				reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ChromosomeFeature'}),
				],
			[1, 'DnaA boxes', 
				models.ChromosomeFeature.objects.filter(species__id = species.id, type__parent__wid='ChromosomeFeature-DnaA_box').count(),
				],				
			[1, 'Short tandem repeats', 
				models.ChromosomeFeature.objects.filter(species__id = species.id, type__parent__wid='ChromosomeFeature-Short_Tandem_Repeat').count(),
				],
			[1, 'Other', 
				models.ChromosomeFeature.objects.filter(species__id = species.id).exclude(type__parent__wid='ChromosomeFeature-DnaA_box').exclude(type__parent__wid='ChromosomeFeature-Short_Tandem_Repeat').count()],
		])
		
		content.append([
			[0, 'Metabolites', models.Metabolite.objects.filter(species__id = species.id).count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Metabolite'})],
			[1, 'Amino acids', 
				models.Metabolite.objects.filter(species__id = species.id, type__wid='amino_acid').count() + 
				models.Metabolite.objects.filter(species__id = species.id, type__wid='modified_amino_acid').count() +
				models.Metabolite.objects.filter(species__id = species.id, type__wid='non-standard_amino_acid').count() +
				models.Metabolite.objects.filter(species__id = species.id, type__wid='vitamin_non-standard_amino_acid').count()				
				],
			[1, 'Antibiotic', 
				models.Metabolite.objects.filter(species__id = species.id, type__wid='antibiotic').count() + 
				models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='antibiotic').count()
				],
			[1, 'Gases', 
				models.Metabolite.objects.filter(species__id = species.id, type__wid='gas').count() + 
				models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='gas').count()
				],
			[1, 'Ions', 
				models.Metabolite.objects.filter(species__id = species.id, type__wid='ion').count() + 
				models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='ion').count()
				],
			[1, 'Lipids', 
				models.Metabolite.objects.filter(species__id = species.id, type__wid='lipid').count() + 
				models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='lipid').count() +
				models.Metabolite.objects.filter(species__id = species.id, type__parent__parent__wid='lipid').count()
				],
			[1, 'Vitamins', 
				models.Metabolite.objects.filter(species__id = species.id, type__wid='vitamin').count() + 
				models.Metabolite.objects.filter(species__id = species.id, type__parent__wid='vitamin').count()
				],
		])
				
		mons = models.ProteinMonomer.objects.filter(species__id = species.id)
		cpxs = models.ProteinComplex.objects.filter(species__id = species.id)
		monDNABind = mons.filter(dna_footprint__length__gt=0).count()
		monIntMem = mons.filter(localization__compartment__wid = 'm').exclude(localization__signal_sequence_type = 'lipoprotein').count() + mons.filter(localization__compartment__wid = 'tm').exclude(localization__signal_sequence_type = 'lipoprotein').count()
		monLipo = mons.filter(localization__signal_sequence_type = 'lipoprotein').count()
		monSecreted = mons.filter(localization__signal_sequence_type = 'secretory').count()
		monTermOrg = mons.filter(localization__compartment__wid = 'tc').count() + mons.filter(localization__compartment__wid = 'tm').count()
		cpxDNABind = cpxs.filter(dna_footprint__length__gt=0).count()
		content.append([
			[0, 'Proteins', mons.count() + cpxs.count()],
				[1, 'Monomers', mons.count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'})],			
					[2, 'DNA-binding', monDNABind],
					[2, 'Integral membrane', monIntMem],
					[2, 'Lipoprotein', monLipo, None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?localization__signal_sequence_type=lipoprotein'],			
					[2, 'Secreted', monSecreted, None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinMonomer'}) + '?localization__signal_sequence_type=secretory'],
					[2, 'Terminal organelle', monTermOrg],					
				[1, 'Complexes', cpxs.count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'ProteinComplex'})],
					[2, 'DNA-binding', cpxDNABind],
		])

		rxns = models.Reaction.objects.filter(species__id = species.id)
		content.append([
			[0, 'Reactions', rxns.count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'})],
			[1, 'DNA damage', rxns.filter(processes__wid='Process_DNADamage').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNADamage'],
			[1, 'DNA repair', rxns.filter(processes__wid='Process_DNARepair').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_DNARepair'],
			[1, 'Metabolic', rxns.filter(processes__wid='Process_Metabolism').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Metabolism'],			
			[1, 'Protein decay', rxns.filter(processes__wid='Process_ProteinDecay').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinDecay'],
			[1, 'Protein modification', rxns.filter(processes__wid='Process_ProteinModification').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ProteinModification'],			
			[1, 'Replication Initiation', rxns.filter(processes__wid='Process_ReplicationInitiation').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_ReplicationInitiation'],
			[1, 'RNA decay', rxns.filter(processes__wid='Process_RNADecay').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNADecay'],
			[1, 'RNA modification', rxns.filter(processes__wid='Process_RNAModification').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAModification'],			
			[1, 'RNA processing', rxns.filter(processes__wid='Process_RNAProcessing').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_RNAProcessing'],			
			[1, 'Transcription', rxns.filter(processes__wid='Process_Transcription').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Transcription'],			
			[1, 'Translation', rxns.filter(processes__wid='Process_Translation').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_Translation'],			
			[1, 'tRNA aminoacylation', rxns.filter(processes__wid='Process_tRNAAminoacylation').count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Reaction'}) + '?processes=Process_tRNAAminoacylation'],
			[1, 'Other', rxns.exclude(processes__wid='Process_DNADamage')
				.exclude(processes__wid='Process_DNARepair')
				.exclude(processes__wid='Process_Metabolism')
				.exclude(processes__wid='Process_ProteinDecay')
				.exclude(processes__wid='Process_ProteinModification')
				.exclude(processes__wid='Process_ReplicationInitiation')				
				.exclude(processes__wid='Process_RNADecay')
				.exclude(processes__wid='Process_RNAModification')
				.exclude(processes__wid='Process_RNAProcessing')
				.exclude(processes__wid='Process_Transcription')
				.exclude(processes__wid='Process_Translation')
				.exclude(processes__wid='Process_tRNAAminoacylation')
				.count()],
		])		
		
		tr = models.TranscriptionalRegulation.objects.filter(species__id = species.id)
		nTus = len(set([x[0] for x in tr.values_list('transcription_unit')]))
		nTfs = len(set([x[0] for x in tr.values_list('transcription_factor')]))
		content.append([
			[0, 'Transcriptional regulation'],
			[1, 'Interactions', tr.count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'TranscriptionalRegulation'})],
			[1, 'Transcriptional regulators', nTfs],
			[1, 'Regulated promoters', nTus],
		])		

		content.append([
			[0, 'Pathways', models.Pathway.objects.filter(species__id = species.id).count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Pathway'})],
		])
		content.append([
			[0, 'Stimuli', models.Stimulus.objects.filter(species__id = species.id).count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Stimulus'})],
		])
		
		
		nCellConc = models.Metabolite.objects.filter(species__id = species.id, intracellular_concentrations__isnull=False).count()
		nMediaConc = models.Metabolite.objects.filter(species__id = species.id, media_concentration__isnull=False).count()		
		nKineticsKeq = models.Reaction.objects.filter(species__id = species.id, keq__isnull=False).count()
		nKineticsKm = \
			models.Reaction.objects.filter(species__id = species.id, kinetics_forward__km__isnull=False).count() + \
			models.Reaction.objects.filter(species__id = species.id, kinetics_backward__km__isnull=False).count()
		nKineticsVmax = \
			models.Reaction.objects.filter(species__id = species.id, kinetics_forward__vmax__isnull=False).count() + \
			models.Reaction.objects.filter(species__id = species.id, kinetics_backward__vmax__isnull=False).count()
		nTranscriptionRates = models.TranscriptionUnit.objects.filter(species__id = species.id, transcription_rate__isnull=False).count()
		nTranslationRates = models.ProteinMonomer.objects.filter(species__id = species.id, translation_rate__isnull=False).count()
		nRnaHl = models.Rna.objects.filter(species__id = species.id, half_life__isnull=False).count()
		nProtHl = models.ProteinMonomer.objects.filter(species__id = species.id, half_life__isnull=False).count()
		nStimuli = models.Stimulus.objects.filter(species__id = species.id, value__isnull=False).count()		
		nTrActivity = models.TranscriptionalRegulation.objects.filter(species__id = species.id, activity__isnull=False).count()
		nTrAffinity = models.TranscriptionalRegulation.objects.filter(species__id = species.id, affinity__isnull=False).count()
		nOther = models.Parameter.objects.filter(species__id = species.id).count()
		nTotParameters = nCellConc + nMediaConc + nKineticsVmax + nTranscriptionRates + nTranslationRates + nRnaHl + nProtHl + nStimuli + nTrAffinity + nTrActivity + nOther

		content.append([
			[0, 'Quantitative parameters', nTotParameters, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid})],
			[1, 'Intracellular concentrations', nCellConc, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'intracellular_concentration'})],
			[1, 'Media concentrations', nMediaConc, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'media_concentration'})],
			[1, 'Protein half-lives', nProtHl, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'prot_half_life'})],			
			[1, 'Reaction K<sub>eq</sub>', nKineticsKeq, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'reaction_keq'})],
			[1, 'Reaction K<sub>m</sub>', nKineticsKm, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'reaction_km'})],
			[1, 'Reaction V<sub>max</sub>', nKineticsVmax, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'reaction_vmax'})],			
			[1, 'RNA half-lives', nRnaHl, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'rna_half_life'})],
			[1, 'Stimulus values', nStimuli, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'stimulus_value'})],
			[1, 'Transcription rates', nTranscriptionRates, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'transcription_rates'})],			
			[1, 'Transcr. reg. activity', nTrActivity, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'treg_activity'})],
			[1, 'Transcr. reg. affinity', nTrAffinity, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'treg_affinity'})],
			[1, 'Translation rates', nTranslationRates, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'translation_rates'})],
			[1, 'Other', nOther, None, reverse('public.views.list_parameters', kwargs={'species_wid': species.wid, 'type': 'other'})],
		])
		
		content.append([
			[0, 'Processes', models.Process.objects.filter(species__id = species.id).count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'Process'})],
		])
		
		content.append([
			[0, 'States', models.State.objects.filter(species__id = species.id).count(), None, reverse('public.views.list', kwargs={'species_wid': species.wid, 'model_type': 'State'})],
		])
		
	nContent = [len(x) for x in content]
	totContent = sum(nContent)
	cum = 0
	idx = 0
	breakIdxs = [0, 0]
	for x in nContent:
		cum += x
		idx += 1
		if cum > totContent * 1/ 3 and breakIdxs[0] == 0:
			breakIdxs[0] = idx
		if cum > totContent * 2 / 3 and breakIdxs[1] == 0:
			breakIdxs[1] = idx		
			
	contentCol1 = []
	contentCol2 = []
	contentCol3 = []
	i = 0
	for x in content[:breakIdxs[0]]:
		i += 1
		for y in x:
			contentCol1.append([i] + y)	
	i = 0
	for x in content[breakIdxs[0]:breakIdxs[1]]:
		i += 1
		for y in x:
			contentCol2.append([i] + y)
	i = 0
	for x in content[breakIdxs[1]:]:
		i += 1
		for y in x:
			contentCol3.append([i] + y)
		
	sources = {
		'total': 0,
		'types': [],
		'dates': [],
		'evidence_parameters': [],
		'evidence_species': [],
		'evidence_media': [],
		'evidence_pH': [],
		'evidence_temperature': [],
	}
	if species is not None:
		refs = models.Reference.objects.filter(species__id = species.id)
		sources['total'] = refs.count()
		sources['types'] = [
				{'type': 'Articles', 'count': refs.filter(species__id = species.id, type__wid='article').count()},
				{'type': 'Books', 'count': refs.filter(species__id = species.id, type__wid='book').count()},
				{'type': 'Thesis', 'count': refs.filter(species__id = species.id, type__wid='thesis').count()},
				{'type': 'Other', 'count': refs.filter(species__id = species.id, type__wid='misc').count()},
			]
		sources['dates'] = refs.filter(year__isnull=False).order_by('year').values('year').annotate(count=Count('year'))
		
			
		nEstimated = models.Parameter.objects.filter(species__id = species.id, value__evidence__is_experimentally_constrained=False).count()
		nExpConstrained = nTotParameters - nEstimated
		sources['evidence_parameters'] = [
			{'type': 'Experimentally constrained', 'count': nExpConstrained},
			{'type': 'Computationally estimated', 'count': nEstimated},
			]
			
		
		sources['evidence_species'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('species').annotate(count = Count('id'))
		sources['evidence_media'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('media').annotate(count = Count('id'))
		sources['evidence_pH'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('pH').annotate(count = Count('id'))
		sources['evidence_temperature'] = models.Evidence.objects.filter(species_component__species__id = species.id).values('temperature').annotate(count = Count('id'))
			
	return render_queryset_to_response(
		request = request, 
		species_wid = species_wid,
		data = {
			'content': [contentCol1, contentCol2, contentCol3],
			'contentRows': range(max(len(contentCol1), len(contentCol2), len(contentCol3))),
			'sources': sources,			
			},		
		templateFile = 'public/index.html')
	
def about(request, species_wid=None):
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		templateFile = 'public/about.html', 
		data = {
			'ROOT_URL': settings.ROOT_URL,
		}
	)	
		
def tutorial(request, species_wid=None):
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		templateFile = 'public/tutorial.html')	

@login_required	
def contributors(request, species_wid=None):
	queryset = User.objects.all().filter(is_active = True)
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		models = [User],
		queryset = queryset,
		templateFile = 'public/contributors.html')	
		
@login_required	
def contributor(request, username, species_wid=None):
	queryset = objectToQuerySet(get_object_or_404(User, username = username), model = User)
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request,
		models = [User],
		queryset = queryset,
		templateFile = 'public/contributor.html')	

def search(request, species_wid = None):
	query = request.GET.get('q', '')
	engine = request.GET.get('engine', 'haystack')
	
	if engine == 'haystack':
		return search_haystack(request, species_wid, query)
	else:
		return search_google(request, species_wid, query)
				
def search_haystack(request, species_wid, query):
	#search
	if species_wid is None:
		species_wid = Species.objects.all()[0].wid
	results = SearchQuerySet().filter(species_wid=species_wid).filter(content=query)
	
	#calculate facets		
	facets = results.facet('model_type')
	tmp = facets.facet_counts()['fields']['model_type']
	modelNameFacet = []
	objectTypes = getObjectTypes(is_public=request.user.is_anonymous())
	models = []
	for tmp2 in tmp:
		modelName = objectTypes[objectTypes.index(tmp2[0])]
		modelNameFacet.append({
			'name':modelName, 
			'verbose_name': getModel(modelName)._meta.verbose_name,
			'count':tmp2[1],
			})
		models.append(getModel(modelName))
	modelNameFacet.sort(lambda x, y:cmp(x['verbose_name'], y['verbose_name']))
	
	#narrow search by facets
	model_type = request.GET.get('model_type', '')
	if model_type:
		results = results.models(getModel(model_type))
		
	#order results
	results = results.order_by('wid')
	
	#convert results to query set
	queryset = EmptyQuerySet()
	for object in results:
		tmp = object.model.objects.none()
		tmp._result_cache.append(object.object)
		queryset = chain(queryset, tmp)
	
	#form response
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		models = models, 
		queryset = queryset, 
		templateFile = 'public/search.html', 
		data = {
			'query': query,
			'engine': 'haystack',
			'model_type': model_type,
			'modelNameFacet': modelNameFacet,
			})

def search_google(request, species_wid, query):
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		templateFile = 'public/googleSearch.html', 
		data = {
			'query': query,
			'engine': 'google',
			})

def list(request, species_wid, model_type):
	try:
		getObjectTypes(is_public=request.user.is_anonymous()).index(model_type)
	except ValueError:
		raise Http404
		
	species = Species.objects.get(wid=species_wid)
	model = getModel(model_type)
	objects = model.objects.all().filter(species__id=species.id)
	
	facet_fields = []	
	for field_full_name in model._meta.facet_fields:
		#facet
		field_names = str(field_full_name).split('__')
		tmp_model = model
		field_verbose_name = []
		for field_name in field_names:
			field = tmp_model._meta.get_field_by_name(field_name)[0]
			field_verbose_name.append(field.verbose_name)
			if isinstance(field, (ForeignKey, ManyToManyField)):
				tmp_model = field.rel.to
		field_verbose_name = ' &#8250; '.join(field_verbose_name)
				
		if isinstance(field, (ForeignKey, ManyToManyField)) and not issubclass(field.rel.to, Entry):
			continue
		
		if isinstance(field, (ForeignKey, ManyToManyField)):
			facet_data = model.objects.filter(species__id=species.id).order_by(field_full_name + '__name').values(field_full_name, field_full_name + '__wid', field_full_name + '__name').annotate(count = Count(field_full_name))
		else:
			facet_data = model.objects.filter(species__id=species.id).order_by(field_full_name).values(field_full_name).annotate(count=Count(field_full_name))
		
		if isinstance(field, (ForeignKey, ManyToManyField)) and issubclass(field.rel.to, Type):
			facet_data_hash = {}			
			for facet in facet_data:
				if facet[field_full_name] is None:
					continue
				facet_data_hash[facet[field_full_name]] = {
					field_full_name: facet[field_full_name],
					field_full_name + '__wid': facet[field_full_name + '__wid'],
					field_full_name + '__name': facet[field_full_name + '__name'],
					'count': None,
					}
			
				for parent in Type.objects.get(id=facet[field_full_name]).get_all_parents():
					facet_data_hash[parent.id] = {
						field_full_name: parent.id,
						field_full_name + '__wid': parent.wid,
						field_full_name + '__name': parent.name,
						'count': None,
						}
			facet_data = [val for x, val in facet_data_hash.iteritems()]			
			facet_data.sort(lambda x, y:cmp(x[field_full_name + '__name'], y[field_full_name + '__name']))
			
			for idx in range(len(facet_data)):
				facet_data[idx]['count'] = Type.objects.get(id=facet_data[idx][field_full_name]).get_all_members().filter(model_type=model_type).count()
			
		facets = []
		for facet in facet_data:
			value = facet[field_full_name]			
			if value is None or unicode(value) == '':
				continue
			
			if isinstance(field, (ForeignKey, ManyToManyField)):
				id = facet[field_full_name + '__wid']
				name = capfirst(facet[field_full_name + '__name'])
			elif (field.choices is not None) and (len(field.choices) > 0) and (not isinstance(field, (BooleanField, NullBooleanField))):	
				id = value
				choices = [x[0] for x in field.choices]
				if id in choices:
					name = field.choices[choices.index(id)][1]
				else:
					name = capfirst(value)
			else:
				id = value
				name = capfirst(value)
			if value is not None and unicode(value) != '':
				facets.append({
					'id': unicode(id), 
					'name': unicode(name),
					'count': facet['count']})
		if len(facets) > 1:
			facet_fields.append({ 
				'name': field_full_name,
				'verbose_name': field_verbose_name, 
				'facets': facets,
				})
	
		#filter
		val = request.GET.get(field_full_name)		
		if val:
			if isinstance(field, (ForeignKey, ManyToManyField)):
				kwargs = {field_full_name + '__wid': val}
			elif isinstance(field, (BooleanField, NullBooleanField)):
				kwargs = {field_full_name: val == 'True'}
			elif isinstance(field, (AutoField, BigIntegerField, DecimalField, FloatField, IntegerField, PositiveIntegerField, PositiveSmallIntegerField, SmallIntegerField)):
				kwargs = {field_full_name: float(val)}
			else:
				kwargs = {field_full_name: val}				
			objects = objects.filter(**kwargs)
			
			if isinstance(field, (ForeignKey, ManyToManyField)) and issubclass(field.rel.to, Type):
				objects = Type.objects.get(wid=val).get_all_members().filter(model_type=model_type).order_by('wid')

	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		models = [model], 
		queryset = objects, 
		templateFile = 'public/list.html', 
		data = {
			'model_type': model_type,
			'model_verbose_name': model._meta.verbose_name,
			'model_verbose_name_plural': model._meta.verbose_name_plural,
			'facet_fields': facet_fields,
			})
	
def list_parameters(request, species_wid, type = None):
	species = Species.objects.get(wid=species_wid)
	
	cellConc = models.Metabolite.objects.filter(species__id = species.id, intracellular_concentrations__isnull=False)
	mediaConc = models.Metabolite.objects.filter(species__id = species.id, media_concentration__isnull=False)
	kineticsKeq = models.Reaction.objects.filter(species__id = species.id, keq__isnull=False)
	kineticsKmFor = models.Reaction.objects.filter(species__id = species.id, kinetics_forward__km__isnull=False)
	kineticsKmRev = models.Reaction.objects.filter(species__id = species.id, kinetics_backward__km__isnull=False)
	kineticsVmaxFor = models.Reaction.objects.filter(species__id = species.id, kinetics_forward__vmax__isnull=False)
	kineticsVmaxRev = models.Reaction.objects.filter(species__id = species.id, kinetics_backward__vmax__isnull=False)
	transcriptionRates = models.TranscriptionUnit.objects.filter(species__id = species.id, transcription_rate__isnull=False)
	translationRates = models.ProteinMonomer.objects.filter(species__id = species.id, translation_rate__isnull=False)
	rnaHl = models.Rna.objects.filter(species__id = species.id, half_life__isnull=False)
	protHl = models.ProteinMonomer.objects.filter(species__id = species.id, half_life__isnull=False)
	stimuli = models.Stimulus.objects.filter(species__id = species.id, value__isnull=False)	
	trActivity = models.TranscriptionalRegulation.objects.filter(species__id = species.id, activity__isnull=False)
	trAffinity = models.TranscriptionalRegulation.objects.filter(species__id = species.id, affinity__isnull=False)
	other = models.Parameter.objects.filter(species__id = species.id)
	
	facets = [		
		{'name': 'Intracellular concentrations', 'count': cellConc.count(), 'id': 'intracellular_concentration'},
		{'name': 'Media concentrations', 'count': mediaConc.count(), 'id': 'media_concentration'},
		{'name': 'Protein half-lives', 'count': protHl.count(), 'id': 'prot_half_life'},			
		{'name': 'Reaction K<sub>eq</sub>', 'count': kineticsKeq.count(), 'id': 'reaction_keq'},
		{'name': 'Reaction K<sub>m</sub>', 'count': kineticsKmFor.count() + kineticsKmRev.count(), 'id': 'reaction_km'},
		{'name': 'Reaction V<sub>max</sub>', 'count': kineticsVmaxFor.count() + kineticsVmaxRev.count(), 'id': 'reaction_vmax'},
		{'name': 'RNA half-lives', 'count': rnaHl.count(), 'id': 'rna_half_life'},
		{'name': 'Stimulus values', 'count': stimuli.count(), 'id': 'stimulus_value'},
		{'name': 'Transcription rates', 'count': transcriptionRates.count(), 'id': 'transcription_rates'},
		{'name': 'Transcr. regulation activity', 'count': trActivity.count(), 'id': 'treg_activity'},
		{'name': 'Transcr. regulation affinity', 'count': trAffinity.count(), 'id': 'treg_affinity'},
		{'name': 'Translation rates', 'count': translationRates.count(), 'id': 'translation_rates'},
		{'name': 'Other', 'count': other.count(), 'id': 'other'},
	]
		
	qs_models = []
	qs = EmptyQuerySet()
	params = []	

	if type is None or type == 'intracellular_concentration':
		qs_models.append(models.Metabolite)
		qs = chain(qs, cellConc)
		for met in cellConc:
			for conc in met.intracellular_concentrations.all():
				constrained = conc.evidence.all().filter(is_experimentally_constrained=True).count() > 0
				refs = EmptyQuerySet()
				for evidence in conc.evidence.all():
					refs = refs | evidence.references.all()
				refs = refs.count()
				params.append({'type': 'intracellular concentration', 'obj': met, 'value': '%.4f' % conc.concentration, 'units': 'mM', 'constrained': constrained, 'n_references': refs})

	if type is None or type == 'media_concentration': 
		qs_models.append(models.Metabolite)
		qs = chain(qs, mediaConc)
		for met in mediaConc:
			conc = met.media_concentration
			constrained = conc.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in conc.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'media concentration', 'obj': met, 'value': '%.4f' % conc.concentration, 'units': 'mM', 'constrained': constrained, 'n_references': refs})
		
	if type is None or type == 'prot_half_life': 
		qs_models.append(models.ProteinMonomer)
		qs = chain(qs, protHl)
		for mon in protHl:
			hl = mon.half_life
			constrained = hl.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in hl.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'half-life', 'obj': mon, 'value': '%.2f' % hl.value, 'units': hl.units, 'constrained': constrained, 'n_references': refs})
			
	if type is None or type == 'reaction_keq':
		qs_models.append(models.Reaction)
		qs = chain(qs, kineticsKeq)
		for rxn in kineticsKeq:
			keq = rxn.keq
			constrained = keq.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in keq.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'forward K<sub>m</sub>', 'obj': rxn, 'value': '%s' % keq.value, 'units': keq.units, 'constrained': constrained, 'n_references': refs})
			
	if type is None or type == 'reaction_km':
		qs_models.append(models.Reaction)
		qs = chain(qs, kineticsKmFor, kineticsKmRev)
		for rxn in kineticsVmaxFor:
			kinetics = rxn.kinetics_forward
			constrained = kinetics.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in kinetics.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'forward K<sub>m</sub>', 'obj': rxn, 'value': '%s' % kinetics.km, 'units': '&mu;M', 'constrained': constrained, 'n_references': refs})
			
		for rxn in kineticsVmaxRev:
			kinetics = rxn.kinetics_backward
			constrained = kinetics.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in kinetics.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'reverse K<sub>m</sub>', 'obj': rxn, 'value': '%s' % kinetics.km, 'units': '&mu;M', 'constrained': constrained, 'n_references': refs})
		
	if type is None or type == 'reaction_vmax':
		qs_models.append(models.Reaction)
		qs = chain(qs, kineticsVmaxFor, kineticsVmaxRev)
		for rxn in kineticsVmaxFor:
			kinetics = rxn.kinetics_forward
			constrained = kinetics.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in kinetics.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'forward v<sub>max</sub>', 'obj': rxn, 'value': '%.2f' % kinetics.vmax, 'units': kinetics.vmax_unit, 'constrained': constrained, 'n_references': refs})
		for rxn in kineticsVmaxRev:
			kinetics = rxn.kinetics_backward
			constrained = kinetics.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in kinetics.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'reverse v<sub>max</sub>', 'obj': rxn, 'value': '%.2f' % kinetics.vmax, 'units': kinetics.vmax_unit, 'constrained': constrained, 'n_references': refs})
		
	if type is None or type == 'rna_half_life':
		qs_models.append(models.Rna)
		qs = chain(qs, rnaHl)
		for rna in rnaHl:
			hl = rna.half_life
			constrained = hl.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in hl.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'half-life', 'obj': rna, 'value': '%.2f' % hl.value, 'units': hl.units, 'constrained': constrained, 'n_references': refs})
			
	if type is None or type == 'stimulus_value':
		qs_models.append(models.Stimulus)
		qs = chain(qs, stimuli)
		for stimulus in stimuli:
			s = stimulus.value
			constrained = s.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in s.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'value', 'obj': stimulus, 'value': '%.2f' % s.value, 'units': s.units, 'constrained': constrained, 'n_references': refs})
			
	if type is None or type == 'transcription_rates':
		qs_models.append(models.TranscriptionUnit)
		qs = chain(qs, transcriptionRates)
		for rna in transcriptionRates:
			tr = rna.translation_rate
			constrained = tr.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in tr.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'transcription rate', 'obj': rna, 'value': '%.2f' % tr.value, 'units': tr.units, 'constrained': constrained, 'n_references': refs})
			
	if type is None or type == 'treg_activity':
		qs_models.append(models.TranscriptionalRegulation)
		qs = chain(qs, trActivity)
		for tr in trActivity:
			activity = tr.activity
			constrained = activity.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in activity.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'fold change activity', 'obj': tr, 'value': '%.2f' % activity.value, 'units': activity.units, 'constrained': constrained, 'n_references': refs})
			
	if type is None or type == 'treg_affinity':
		qs_models.append(models.TranscriptionalRegulation)
		qs = chain(qs, trAffinity)
		for tr in trAffinity:
			affinity = tr.affinity
			constrained = affinity.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in affinity.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'affinity', 'obj': tr, 'value': '%.2f' % affinity.value, 'units': affinity.units, 'constrained': constrained, 'n_references': refs})			
	
	if type is None or type == 'translation_rates':
		qs_models.append(models.ProteinMonomer)
		qs = chain(qs, translationRates)
		for mon in translationRates:
			tr = mon.translation_rate
			constrained = tr.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in tr.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'translation rate', 'obj': mon, 'value': '%.2f' % tr.value, 'units': tr.units, 'constrained': constrained, 'n_references': refs})
		
	if type is None or type == 'other':
		qs_models.append(models.Parameter)
		qs = chain(qs, other)
		for param in other:
			val = param.value
			constrained = val.evidence.all().filter(is_experimentally_constrained=True).count() > 0
			refs = EmptyQuerySet()
			for evidence in val.evidence.all():
				refs = refs | evidence.references.all()
			refs = refs.count()
			params.append({'type': 'value', 'obj': param, 'value': '%s' % val.value, 'units': val.units, 'constrained': constrained, 'n_references': refs})
		
	qs_models = set(qs_models)
	
	params.sort(lambda x, y:cmp(x['type']+x['obj'].wid, y['type']+y['obj'].wid))
	
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		models = qs_models, 
		queryset = qs, 
		templateFile = 'public/list_parameters.html', 
		data = {
			'params': params,
			'type': type,
			'facets': facets,
			})
	
def detail(request, species_wid, wid):
	obj = getEntry(species_wid = species_wid, wid = wid)
	if obj is None:
		raise Http404
		
	model = obj.__class__
	model_type = model.__name__
	fieldsets = deepcopy(model._meta.fieldsets)
	
	#filter out type, metadata
	fieldset_names = [x[0] for x in fieldsets]
	if 'Type' in fieldset_names:
		idx = fieldset_names.index('Type')
		del fieldsets[idx]
		
	#filter out empty fields
	rmfieldsets = []
	for idx in range(len(fieldsets)):
		rmfields = []
		for idx2 in range(len(fieldsets[idx][1]['fields'])):
			if isinstance(fieldsets[idx][1]['fields'][idx2], dict):
				field = fieldsets[idx][1]['fields'][idx2]
				field_name = field['name']
				verbose_name = field['verbose_name']
			else:
				field_name = fieldsets[idx][1]['fields'][idx2]
				field = model._meta.get_field_by_name(field_name)[0]
				if isinstance(field, RelatedObject):
					verbose_name = capfirst(field.get_accessor_name())
				else:
					verbose_name = field.verbose_name
				
			data = format_field_detail_view(obj, field_name, request.user.is_anonymous())
			if (data is None) or (data == ''):
				rmfields = [idx2] + rmfields
				
			if issubclass(model, Parameter) and field_name == 'value':
				is_modeled_url = reverse('public.views.viewParameterInSimulation', kwargs = {
					'species_wid':species_wid, 
					'wid': wid,
					})
			elif (isinstance(field, (ForeignKey, ManyToManyField)) and issubclass(field.rel.to, Parameter)) or \
			   (isinstance(field, (RelatedObject)) and issubclass(field.model, Parameter)):
				is_modeled_url = reverse('public.views.viewParametersInSimulation', kwargs = {
					'species_wid':species_wid, 
					'wid': wid,
					})
			else:
				try:
					if ModelProperty.objects.get(species__wid = species_wid, class_name = model_type, property_name = field_name).simulation_properties.exists():
						is_modeled_url = reverse('public.views.viewPropertyInSimulation', kwargs = {
							'species_wid':species_wid, 
							'class_name': model_type, 
							'property_name': field_name,
							})
				except ObjectDoesNotExist:
					is_modeled_url = ''
			
			fieldsets[idx][1]['fields'][idx2] = {
				'name': field_name,
				'verbose_name': verbose_name.replace(" ", '&nbsp;').replace("-", "&#8209;"), 
				'data': data,
				'is_modeled_url': is_modeled_url,
				}
		for idx2 in rmfields:
			del fieldsets[idx][1]['fields'][idx2]
		if len(fieldsets[idx][1]['fields']) == 0:
			rmfieldsets = [idx] + rmfieldsets
	for idx in rmfieldsets:
		del fieldsets[idx]
	
	#form query set
	qs = objectToQuerySet(obj, model = model)

	#render response
	return render_queryset_to_response(
		species_wid = species_wid,		
		request = request, 
		models = [model],
		queryset = qs,
		templateFile = 'public/detail.html', 
		data = {
			'model_type': model_type,
			'model': model,
			'fieldsets': fieldsets,
			'message': request.GET.get('message', ''),
			})

def viewPropertyInSimulation(request, species_wid, class_name, property_name):
	#get model
	model = getModel(class_name)
	
	#get verbose class name
	verbose_class_name = model._meta.verbose_name
	
	#get verbose property name
	verbose_property_name = property_name
	for fieldset in model._meta.fieldsets:
		for field in fieldset[1]['fields']:
			if isinstance(field, (str, unicode)) and field == property_name:
				tmp = model._meta.get_field_by_name(property_name)
				if len(tmp) > 0:
					verbose_property_name = tmp[0].verbose_name
			if isinstance(field, dict) and field['name'] == property_name:
				verbose_property_name = field['verbose_name']

	#get associated simulation code properties
	qs = ModelProperty.objects.get(
		species__wid = species_wid,
		class_name = class_name,
		property_name = property_name
		).simulation_properties.all().order_by('class_name', 'property_name')
		
	#organize simulation code properties by class
	classes = {}
	for object in qs:
		if not classes.has_key(object.class_name):
			classes[object.class_name] = []
		classes[object.class_name].append(object.property_name)
		
	#highlight code for each simulation class
	object_list = []
	for sim_class_name in classes:
		sim_property_names = classes[sim_class_name]
		sim_property_names.sort()
		pathParts = sim_class_name.split('.')
		codePath = "%s/src/+%s/%s.m" % (MODEL_CODE_BASE_DIR, '/+'.join(pathParts[0:-1]), pathParts[-1])
		if not os.path.isfile(codePath):
			codePath = "%s/src/+%s/@%s/%s.m" % (MODEL_CODE_BASE_DIR, '/+'.join(pathParts[0:-1]), pathParts[-1], pathParts[-1])
			if not os.path.isfile(codePath):
				continue
		
		with open (codePath, "r") as codeFile:
			code = codeFile.read()

		lexer = MatlabLexer()
		lexer.add_filter(PropertyDefinitionFilter(property_names = sim_property_names, tokentype=Token.Name.Variable)) 
		
		tokens = lexer.get_tokens(code)
			
		object_list.append({
			'class_name': sim_class_name,
			'property_names': sim_property_names,
			'code': pygments.format(tokens, PygmentsFormatter(linenos='inline', linenostep=1, style=PygmentsStyle, noclasses=True)),
			})
	
	#render response
	return render_queryset_to_response(
		species_wid = species_wid,		
		request = request, 
		models = [SimulationProperty],
		queryset = qs,
		templateFile = 'public/viewPropertyInSimulation.html', 
		data = {
			'object_list': object_list,
			'verbose_class_name': verbose_class_name,
			'verbose_property_name': verbose_property_name,
			})
			
def viewParameterInSimulation(request, species_wid, wid):
	#get associated simulation property
	qs = Parameter.objects.filter(species__wid = species_wid, wid=wid)
	if not qs[0].state is None:
		sim_class_name = 'edu.stanford.covert.cell.sim.state.%s' % qs[0].state.wid.replace('State_', '')
		verbose_class_name = '%s: %s' % (wid, qs[0].state.name)
	else:
		sim_class_name = 'edu.stanford.covert.cell.sim.process.%s' % qs[0].process.wid.replace('Process_', '')
		verbose_class_name = '%s: %s' % (wid, qs[0].process.name)
	sim_property_name = qs[0].name
	verbose_property_name = qs[0].name

	#highlight code for simulation class
	pathParts = sim_class_name.split('.')
	codePath = "%s/src/+%s/%s.m" % (MODEL_CODE_BASE_DIR, '/+'.join(pathParts[0:-1]), pathParts[-1])
	if not os.path.isfile(codePath):
		codePath = "%s/src/+%s/@%s/%s.m" % (MODEL_CODE_BASE_DIR, '/+'.join(pathParts[0:-1]), pathParts[-1], pathParts[-1])
	
	if os.path.isfile(codePath):
		with open (codePath, "r") as codeFile:
			code = codeFile.read()

		lexer = MatlabLexer()
		lexer.add_filter(PropertyDefinitionFilter(property_names = [sim_property_name], tokentype=Token.Name.Variable)) 
		
		tokens = lexer.get_tokens(code)
			
		object = {
			'class_name': sim_class_name,
			'property_names': [sim_property_name],
			'code': pygments.format(tokens, PygmentsFormatter(linenos='inline', linenostep=1, style=PygmentsStyle, noclasses=True)),
			}
	else:
		raise Http404
	
	#render response
	return render_queryset_to_response(
		species_wid = species_wid,		
		request = request, 
		models = [Parameter],
		queryset = qs,
		templateFile = 'public/viewPropertyInSimulation.html', 
		data = {
			'object_list': [object],
			'verbose_class_name': verbose_class_name,
			'verbose_property_name': verbose_property_name,
			})
	
def viewParametersInSimulation(request, species_wid, wid):
	obj = getEntry(species_wid = species_wid, wid = wid)
	if obj is None:
		raise Http404
		
	model = obj.__class__
	verbose_class_name = model._meta.verbose_name
	verbose_property_name = '%s parameters' % obj.name
	
	#get associated simulation property
	qs = EmptyQuerySet()
	for field in model._meta.get_all_related_objects() + model._meta.get_all_related_many_to_many_objects():
		accessor_name = field.get_accessor_name()
		if accessor_name == '+':
			continue

		if isinstance(field, ForeignKey) and issubclass(field.rel.to, Parameter):
			qs._result_cache.append(getattr(obj, accessor_name))
		if (isinstance(field, ManyToManyField) and issubclass(field.rel.to, Parameter)) or \
		   (isinstance(field, RelatedObject) and issubclass(field.model, Parameter)):
			qs2 = getattr(obj, accessor_name).all()
			if len(qs2) > 0:
				qs._result_cache.extend(qs2._result_cache)
	
	classes = {}	
	for o in qs:
		if not o.state is None:
			sim_class_name = 'edu.stanford.covert.cell.sim.state.%s' % o.state.wid.replace('State_', '')
		else:
			sim_class_name = 'edu.stanford.covert.cell.sim.process.%s' % o.process.wid.replace('Process_', '')
		if not classes.has_key(sim_class_name):
			classes[sim_class_name] = []			
		classes[sim_class_name].append(o.name)

	#highlight code for simulation class
	objects = []
	for sim_class_name in classes:
		sim_property_names = classes[sim_class_name]
		sim_property_names.sort()
		pathParts = sim_class_name.split('.')
		codePath = "%s/src/+%s/%s.m" % (MODEL_CODE_BASE_DIR, '/+'.join(pathParts[0:-1]), pathParts[-1])
		if not os.path.isfile(codePath):
			codePath = "%s/src/+%s/@%s/%s.m" % (MODEL_CODE_BASE_DIR, '/+'.join(pathParts[0:-1]), pathParts[-1], pathParts[-1])
			if not os.path.isfile(codePath):
				continue
		with open (codePath, "r") as codeFile:
			code = codeFile.read()

		lexer = MatlabLexer()
		lexer.add_filter(PropertyDefinitionFilter(property_names = sim_property_names, tokentype=Token.Name.Variable)) 
		
		tokens = lexer.get_tokens(code)
			
		objects.append({
			'class_name': sim_class_name,
			'property_names': sim_property_names,
			'code': pygments.format(tokens, PygmentsFormatter(linenos='inline', linenostep=1, style=PygmentsStyle, noclasses=True)),
			})
	
	#render response
	return render_queryset_to_response(
		species_wid = species_wid,		
		request = request, 
		models = [Parameter],
		queryset = qs,
		templateFile = 'public/viewPropertyInSimulation.html', 
		data = {
			'object_list': objects,
			'verbose_class_name': verbose_class_name,
			'verbose_property_name': verbose_property_name,
			})
			
@login_required		
def add(request, model_type, species_wid=None):
	return edit(request, model_type=model_type, species_wid=species_wid, action='add')
		
@login_required
def edit(request, wid=None, model_type=None, species_wid=None, action='edit'):
	#retrieve object
	if action == 'edit':
		obj = getEntry(species_wid = species_wid, wid = wid)
		if obj is None:
			raise Http404
		model = obj.__class__ 
	else:		
		model = getModel(model_type)
		obj = model()
	
	#save object
	error_messages = {}
	if request.method == 'POST':
		submitted_data = get_edit_form_data(model, request.POST)
		
		data = submitted_data
		data['id'] = obj.id
		data['species'] = species_wid
		data['model_type'] = model.__name__
		
		try:
			#validate is WID unique
			if issubclass(model, SpeciesComponent):
				qs = SpeciesComponent.objects.values('wid', 'model_type').filter(species__wid=species_wid)
			else:
				qs = model.objects.values('wid', 'model_type').all()
				
			if action == 'edit':
				qs = qs.exclude(id=obj.id)
				
			wids = {}
			for x in qs:
				wids[x['wid']] = x['model_type']
			
			if data['wid'] in wids.keys():
				raise ValidationError({'wid': 'Value must be unique'})
				
			wids[data['wid']] = model.__name__
		
			#validate
			data = validate_object_fields(model, data, wids, species_wid, data['wid'])
			validate_model_objects(model, data)
			validate_model_unique(model, [data])
			
			#save
			obj = save_object_data(species_wid, obj, data, {}, request.user, save=False, save_m2m=False)
			obj = save_object_data(species_wid, obj, data, {data['wid']: obj}, request.user, save=True, save_m2m=False)
			obj = save_object_data(species_wid, obj, data, {data['wid']: obj}, request.user, save=True, save_m2m=True)
			
			#redirect to details page
			return HttpResponseRedirect(obj.get_absolute_url())
		except ValidationError as error:
			error_messages = error.message_dict
	
	#form query set
	if action == 'edit':
		obj = getEntry(species_wid = species_wid, wid = wid)
		if obj is None:
			raise Http404
		qs = objectToQuerySet(obj, model = model)
	else:
		obj = None
		qs = model.objects.none()
		
	#display form
	fields, initial_values = get_edit_form_fields(species_wid, model, obj=obj)
	
	if request.method == 'POST':
		initial_values = submitted_data
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		models = [model],
		queryset = qs,
		templateFile = 'public/edit.html', 
		data = {
			'model_verbose_name': model._meta.verbose_name,
			'action': action,
			'fields': fields,
			'references_choices': Reference.objects.filter(species__wid = species_wid).values_list('wid'),
			'initial_values': initial_values,
			'error_messages': error_messages,
			}
		)
			
@login_required		
def delete(request, species_wid, wid):
	#retrieve object
	obj = getEntry(species_wid = species_wid, wid = wid)
	if obj is None:
		raise Http404	
	model = obj.__class__ 
	qs = objectToQuerySet(obj, model = model)
	
	#delete
	if request.method == 'POST':
		obj.delete()
		return HttpResponseRedirect(reverse('public.views.list', kwargs={'species_wid':species_wid, 'model_type': model.__name__}))
		
	#confirmation message
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		models = [model],
		queryset = qs,
		templateFile = 'public/delete.html', 
		data = {
			'model_verbose_name': model._meta.verbose_name
			}
		)

def exportData(request, species_wid=None):	
	getDict = request.GET.copy()
	if getDict.get('format', ''):
		getDict.__setitem__('species', getDict.get('species', species_wid))	
	form = ExportDataForm(request.user.is_anonymous(), getDict or None)
	if not form.is_valid():		
		return render_queryset_to_response(
			species_wid=species_wid,
			request = request,
			templateFile = 'public/exportDataForm.html', 
			data = {
				'form': form
				}
			)
	else:		
		species = Species.objects.get(wid = form.cleaned_data['species'])
		queryset = EmptyQuerySet()
		models = []
		if form.cleaned_data['all_model_types'] == 'True':			
			model_types = getObjectTypes(superclass=Entry, is_public = request.user.is_anonymous())
		else:
			model_types = form.cleaned_data['model_type']
		
		for model_type in model_types:
			model = getModel(model_type)
			if issubclass(model, SpeciesComponent):
				queryset = chain(queryset, model.objects.filter(species__id=species.id).select_related(depth=2).all())
			else:
				queryset = chain(queryset, model.objects.select_related(depth=2).filter(id=species.id))
			models.append(getModel(model_type))
		
		return render_queryset_to_response(
			species_wid = species_wid,
			request = request, 
			queryset = queryset, 
			templateFile = 'public/exportDataResult.html', 
			models = models)
			
def exportDataTemplate(request):
	queryset = EmptyQuerySet()
	models = []
	
	for model_type in getObjectTypes():
		models.append(getModel(model_type))
		
	#write work book
	wb = writeExcel(
		Species(wid = '__NEW__', name = 'New species'),
		queryset, models, request.user.is_anonymous())

	#save to string
	result = StringIO()
	wb.save(filename = result)

	#generate HttpResponse
	response = HttpResponse(
		result.getvalue(),
		mimetype = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
		content_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
	response['Content-Disposition'] = "attachment; filename=data.xlsx"
	return response

@login_required
def importData(request, species_wid=None):
	if request.method == 'POST':
		form = ImportDataForm(request.POST or None, request.FILES)
		
		if form.is_valid():		
			selected_species_wid = form.cleaned_data['species'] or form.cleaned_data['new_species']
			if selected_species_wid == '' or selected_species_wid is None:
				form.species.errors += ['Please select a specices']	
			else:
				#save to temporary file
				originalFileName, originalFileExtension = os.path.splitext(request.FILES['file'].name)
				fid = tempfile.NamedTemporaryFile(suffix = originalFileExtension, delete = False)
				filename = fid.name
				for chunk in request.FILES['file'].chunks():
					fid.write(chunk)
				fid.close()
				
				#read file
				if originalFileExtension == '.xlsx':
					try:
						batch_import_from_excel(selected_species_wid, filename, request.user)
						success = True
						message = 'Data successfully saved!'
					except ValidationError as error:
						success = False
						message = 'Unable to import data: ' + ' '.join(error.messages)
				elif originalFileExtension == '.fna':
					success, message = readFasta(selected_species_wid, filename, request.user)
				else:
					raise Http404
				
				#delete file
				os.remove(filename)
				
				#render response
				return render_queryset_to_response(
					species_wid = species_wid,
					request = request,
					templateFile = 'public/importDataResult.html', 
					data = {
						'success': success,
						'message': message,
						})
	else:
		form = ImportDataForm(None)

	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		templateFile = 'public/importDataForm.html', 
		data = {
			'form': form},
		)
	
def validate(request, species_wid):
	errors = get_invalid_objects(Species.objects.values('id').get(wid=species_wid)['id'])
	
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		templateFile = 'public/validate.html', 
		data = {			
			'errors': errors
			},
		)
		
@sensitive_post_parameters()
@csrf_protect
@never_cache
def login(request, species_wid=None):
	next = request.REQUEST.get('next', '')
	
	if request.method == "POST":
		form = AuthenticationForm(data=request.POST)
		if form.is_valid():
			auth_login(request, form.get_user())
			
			if request.session.test_cookie_worked():
				request.session.delete_test_cookie()

			return HttpResponseRedirect(next)
	else:
		form = AuthenticationForm(request)

	request.session.set_test_cookie()

	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		templateFile = 'public/login.html', 
		data = {
			'form': form,
			'next': next,
		})
		
def logout(request, species_wid=None):
	auth_logout(request)	
	return render_queryset_to_response(
		species_wid = species_wid,
		request = request, 
		templateFile = 'public/logout.html', 
		)
	
def sitemap(request):
	return render_queryset_to_response(
		request = request, 
		templateFile = 'public/sitemap.xml', 
		data = {
			'ROOT_URL': settings.ROOT_URL,
			'qs_species': Species.objects.all(),
		}
	)
	
def sitemap_toplevel(request):
	return render_queryset_to_response(
		request = request, 
		templateFile = 'public/sitemap_toplevel.xml', 
		data = {
			'ROOT_URL': settings.ROOT_URL,
		}
	)
	
def sitemap_species(request, species_wid):
	species = Species.objects.get(wid=species_wid)
	return render_queryset_to_response(
		request = request, 
		templateFile = 'public/sitemap_species.xml', 
		data = {
			'ROOT_URL': settings.ROOT_URL,
			'species': species,
			'entries': SpeciesComponent.objects.filter(species__id = species.id),
		}
	)
