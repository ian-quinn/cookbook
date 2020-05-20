from flask import Flask
from flask import request, current_app
from config import Config
from flask_migrate import Migrate

import os, sys

# application factory functions:
app = Flask(__name__)
app.config.from_object(Config)


from app import routes
