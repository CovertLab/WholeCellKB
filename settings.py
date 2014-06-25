'''
Whole-cell project settings

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

import os

BASE_URL = 'http://www.wholecell.org'
ROOT_URL = '/mpn/kb'
FULL_URL = BASE_URL + ROOT_URL
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
MODEL_CODE_BASE_DIR = '/home/wholecell/wholecell-mpn/sim'
GOOGLE_SEARCH_ENABLED = False

DEBUG = False
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    #('Jonathan Karr', 'jkarr@stanford.edu'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'wholecellkbmpn',
        'USER': 'wholecellkb',
        'PASSWORD': 'wholecellkb',
        'HOST': 'localhost',
        'PORT': '',
    }
}

ALLOWED_HOSTS = ['wholecelldb.org', 'www.wholecelldb.org']

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# On Unix systems, a value of None will cause Django to use the same
# timezone as the operating system.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'America/Los_Angeles'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/media/"
MEDIA_ROOT = ''

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://media.lawrence.com/media/", "http://example.com/media/"
MEDIA_URL = ''

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = ''

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = ROOT_URL + '/static/'

# Additional locations of static files, use absolute paths
STATICFILES_DIRS = (
    ROOT_DIR + "/static",
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = '!!l%%3iao+m!@q51vzr3a0+(nk_w7d*_c1u-7re!3tlbmh*mgs'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
)

ROOT_URLCONF = os.path.dirname(os.path.realpath(__file__)).split(os.sep)[-1] + '.urls'

#use absolute directories
TEMPLATE_DIRS = (
	ROOT_DIR + "/templates",
)

import os
HAYSTACK_SITECONF = 'public.search_indexes'
HAYSTACK_SEARCH_ENGINE  = 'xapian'
HAYSTACK_XAPIAN_PATH = os.path.join(os.path.dirname(__file__), 'xapian_index')

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',	
    'django.contrib.admin',
    'django.contrib.admindocs',
	'django_extensions',
	
	#apps
	'public',
	
	#helpers
	'haystack',
)

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'class': 'django.utils.log.AdminEmailHandler'
        }
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
    }
}

AUTH_PROFILE_MODULE = 'public.UserProfile'
LOGIN_URL = ROOT_URL + '/login/'
LOGIN_REDIRECT_URL = ROOT_URL + '/'

ABSOLUTE_URL_OVERRIDES = {
    'auth.user': lambda o: ROOT_URL + '/contributor/' +  o.username,
}
