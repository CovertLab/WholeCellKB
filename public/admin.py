'''
Whole-cell knowledge base admin interface

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User
from public.models import UserProfile

''' User profile admin '''
class UserProfileInline(admin.StackedInline):
	model = UserProfile

class UserProfileAdmin(UserAdmin):
	inlines = [ UserProfileInline, ]

''' Register admins '''
admin.site.unregister(User)
admin.site.register(User, UserProfileAdmin)