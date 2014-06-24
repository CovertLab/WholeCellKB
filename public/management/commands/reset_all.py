from django.core.management import call_command
from django.core.management.base import BaseCommand, CommandError

class Command(BaseCommand):
    args = ''
    help = 'Resets the SQL database, rebuilds search index'

    def handle(self, *args, **options):
        interactive = options['interactive'] if 'interactive' in options else True
        
        #reset SQL db
        call_command('reset', 'public', interactive=interactive)        
        
        #rebuild search index
        call_command('rebuild_index', interactive=interactive)
            
        #status message
        self.stdout.write('Successfully reset SQL and HDF databases and search indices\n')