import django.core.handlers.wsgi
import os
import sys

paths = [
	os.path.dirname(os.path.realpath(__file__)) + '/../..',
	os.path.dirname(os.path.realpath(__file__)) + '/..',
	]
for path in paths:
	if path not in sys.path:
		sys.path.append(path)

os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

application = django.core.handlers.wsgi.WSGIHandler()

#change monitor
import monitor
monitor.start(interval = 1.0)
