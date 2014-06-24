from django.core.management.base import BaseCommand, CommandError
from kb import settings
import os
import subprocess

class Command(BaseCommand):
    args = ''
    help = 'Compiles documentation'

    def handle(self, *args, **options):
        #generate html documentation
        subprocess.call('epydoc -o %s --name=WholeCellKB --url=%s --html --graph all --parse-only --docformat plaintext public' % (
            os.path.join(settings.ROOT_DIR, 'static', 'public', 'doc'), settings.ROOT_URL,
            ), shell=True)
        
        #make images of data model
        subprocess.call('python manage.py graph_models public | dot -Tsvg -o %s' % os.path.join(settings.ROOT_DIR, 'static', 'public', 'img', 'data_model.svg'), shell=True)
        subprocess.call('python manage.py graph_models public | dot -Tpng -o %s' % os.path.join(settings.ROOT_DIR, 'static', 'public', 'img', 'data_model.png'), shell=True)
        
        #status message
        self.stdout.write('Successfully compiled documentation.\n')