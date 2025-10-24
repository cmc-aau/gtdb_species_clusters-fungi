import os
import logging

from rich.logging import RichHandler

# Set the module version.
with open(os.path.join(__path__[0], 'VERSION'), 'r') as f:
    __version__ = f.readline().strip()
    
FORMAT = "%(message)s"
logging.basicConfig(
    level=logging.INFO, format=FORMAT, datefmt="[%Y-%m-%d %H:%M:%S]", handlers=[RichHandler()]
)

logging.getLogger("matplotlib").setLevel(logging.WARNING)
