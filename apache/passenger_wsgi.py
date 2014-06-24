import sys, os
sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), '..', '..'))
sys.path.append(os.path.join(os.getcwd(), '..'))
os.environ['DJANGO_SETTINGS_MODULE'] = "kb.settings"
import django.core.handlers.wsgi
_application = django.core.handlers.wsgi.WSGIHandler()

from kb import monitor
monitor.start(interval = 1.0)

from paste.exceptions.errormiddleware import ErrorMiddleware
application = ErrorMiddleware(_application, debug=True)